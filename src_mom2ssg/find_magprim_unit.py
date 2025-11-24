#!/usr/bin/env python3
"""
Magnetic Primitive Unit Finding Module

This module provides functions for finding and processing magnetic primitive units
in crystal structures. It handles magnetic moment clustering, species relabeling,
and primitive cell identification for magnetic materials.

Main functionality:
- Relabel atoms by species and magnetic moments
- Find magnetic primitive units
- Reorder crystal structures
- Handle magnetic moment clustering and mapping

Author: MOM2SSG Team
"""

import numpy as np
import copy
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from .getcell import prim_cell

def relabel_by_species_and_magmoms(cell, tol=1e-3, prefix="r", per_component=False):
    """
    Relabel atoms by clustering species and magnetic moments.
    
    This function groups atoms with similar magnetic moments within the same species
    and creates new labels for these clusters. It's essential for identifying
    magnetic primitive units in complex magnetic structures.
    
    Args:
        cell (tuple): Crystal structure containing (lattice, positions, numbers, elements, mag)
            - lattice: 3x3 lattice vectors matrix
            - positions: Nx3 fractional coordinates
            - numbers: Nx1 atom type indices
            - elements: List of element symbols
            - mag: Nx3 magnetic moment vectors
        tol (float): Tolerance for magnetic moment clustering (default: 1e-3)
        prefix (str): Prefix for new labels (default: "r")
        per_component (bool): If True, use max component difference; if False, use L2 norm
    
    Returns:
        tuple: (cell_new, mapping) containing:
            - cell_new: Relabeled crystal structure
            - mapping: Dictionary mapping new labels to species, magnetic moments, and counts
    
    Note:
        The function clusters atoms with similar magnetic moments within each species
        and assigns new labels based on the clustering results.
    """
    # Extract crystal structure components
    magmoms = np.asarray(cell[4], float)
    type_ids = np.asarray(cell[2], int)
    species_list = np.asarray(cell[3], str)

    # Get species for each site
    site_species = [species_list[i-1] for i in type_ids]

    # Ensure magnetic moments are 2D array
    m_arr = np.asarray(magmoms, dtype=float)
    if m_arr.ndim == 1:
        m_arr = m_arr[:, None]
    N, D = m_arr.shape

    def distance(a, b):
        """
        Calculate distance between two magnetic moment vectors.
        
        Args:
            a, b: Magnetic moment vectors
            
        Returns:
            float: Distance between vectors
        """
        if per_component:
            return np.max(np.abs(a - b))
        else:
            return np.linalg.norm(a - b)

    # Initialize clustering data structures
    species_clusters = {}  # Dictionary to store clusters by species
    creation_order = []    # Order of cluster creation
    
    # Cluster atoms by species and magnetic moments
    for idx, (sp, v) in enumerate(zip(site_species, m_arr)):
        v = v.ravel()
        clusters = species_clusters.setdefault(sp, [])

        # Find best matching cluster
        best_j = None
        best_d = np.inf
        for j, c in enumerate(clusters):
            d = distance(v, c["centroid"])
            if d <= tol and d < best_d:
                best_d, best_j = d, j

        if best_j is None:
            # Create new cluster
            c = {
                "sum": v.copy(),
                "count": 1,
                "centroid": v.copy(),
                "members": [idx],
                "label": None, 
                "species": sp,
            }
            clusters.append(c)
            creation_order.append(c)
        else:
            # Add to existing cluster
            c = clusters[best_j]
            c["sum"] += v
            c["count"] += 1
            c["centroid"] = c["sum"] / c["count"]
            c["members"].append(idx)

    # Assign labels to clusters
    for k, c in enumerate(creation_order, start=1):
        c["label"] = f"{prefix}{k}"

    # Create new labels for all atoms
    new_labels = [None] * N
    for sp, clusters in species_clusters.items():
        for c in clusters:
            for i in c["members"]:
                new_labels[i] = c["label"]
    
    # Get ordered list of labels
    order = [c["label"] for c in creation_order]

    # Create mapping dictionary
    mapping = {}
    for c in creation_order:
        mapping[c["label"]] = {
            "species": c["species"],
            "magmom": c["centroid"].copy(),
            "count": c["count"],
        }

    # Create new cell with relabeled atoms
    cell_new = list(copy.deepcopy(cell))
    cell_new[3] = order 

    # Convert string labels to integer indices
    for i in range(len(order)):
        new_labels = replace_in_list(new_labels, prefix + str(i+1), np.int32(i+1))
    cell_new[2] = new_labels

    return cell_new, mapping


def reorder_cell(cell):
    """
    Reorder atoms in the cell by their type IDs.
    
    This function sorts all atoms in the crystal structure by their type IDs,
    ensuring that atoms of the same type are grouped together. This is useful
    for consistent processing and analysis of crystal structures.
    
    Args:
        cell (tuple): Crystal structure containing (lattice, positions, numbers, elements, mag)
            - lattice: 3x3 lattice vectors matrix
            - positions: Nx3 fractional coordinates
            - numbers: Nx1 atom type indices
            - elements: List of element symbols
            - mag: Nx3 magnetic moment vectors
    
    Returns:
        tuple: Reordered crystal structure with atoms sorted by type IDs
    """
    cell_new = list(copy.deepcopy(cell))
    
    # Extract crystal structure components
    frac_coords = np.asarray(cell[1], float)
    magmoms = np.asarray(cell[4], float)
    type_ids = np.asarray(cell[2], int)
    species_list = np.asarray(cell[3], str)
    
    # Sort by type IDs
    idx = np.argsort(type_ids)
    frac_coords = frac_coords[idx]
    magmoms = magmoms[idx]
    type_ids = type_ids[idx]
    
    # Update cell with reordered data
    cell_new[1] = frac_coords
    cell_new[2] = type_ids
    cell_new[4] = magmoms
    cell_new[3] = species_list
    
    return cell_new
    

def relabel_by_species_and_magmoms_back(cell_, mapping, prefix="r"):
    """
    Restore original species labels and magnetic moments from mapping.
    
    This function reverses the relabeling process by restoring the original
    species labels and magnetic moments using the mapping dictionary created
    during the clustering process.
    
    Args:
        cell_ (tuple): Crystal structure with relabeled atoms
        mapping (dict): Mapping dictionary from relabel_by_species_and_magmoms
        prefix (str): Prefix used in the original relabeling (default: "r")
    
    Returns:
        tuple: Crystal structure with restored species labels and magnetic moments
    """
    cell = reorder_cell(cell_)
    cell_new = list(copy.deepcopy(cell))
    
    # Extract crystal structure components
    magmoms = np.asarray(cell[4], float)
    type_ids = np.asarray(cell[2], int)
    species_list = np.asarray(cell[3], str)

    # Create mapping from new labels to original species
    id_i = 0
    id_j = 0
    iddict = {}
    atoms = []
    
    for i in species_list:
        id_j += 1
        atom = mapping[i]['species']
        if atom not in atoms:
            id_i += 1
            atoms.append(atom)
        iddict[id_j] = id_i
    
    # Restore original type IDs and magnetic moments
    cell_new[2] = []
    id_i = 0
    
    for i in type_ids:
        cell_new[2].append(iddict[i])
        cell_new[-1][id_i] = mapping[prefix + str(i)]['magmom']
        id_i = id_i + 1
    
    cell_new[3] = atoms
    
    return cell_new

def find_magprim_unit(cell):
    """
    Find the magnetic primitive unit cell.
    
    This is the main function that identifies the magnetic primitive unit cell
    by clustering atoms with similar magnetic moments, finding the primitive
    cell, and applying the appropriate transformations to restore the original
    magnetic structure.
    
    Args:
        cell (tuple): Crystal structure containing (lattice, positions, numbers, elements, mag)
            - lattice: 3x3 lattice vectors matrix
            - positions: Nx3 fractional coordinates
            - numbers: Nx1 atom type indices
            - elements: List of element symbols
            - mag: Nx3 magnetic moment vectors
    
    Returns:
        tuple: Magnetic primitive unit cell with transformed magnetic moments
    
    Note:
        The function uses Phonopy to find the primitive cell and applies
        the standardization rotation matrix to the magnetic moments.
        If the input cell only contains (lattice, positions, numbers, mag),
        placeholder element labels are generated to keep numeric type IDs aligned.
    """
    # Step 0: Normalize inputs that lack an explicit element list
    if len(cell) == 4:
        lattice, positions, numbers, mag = cell
        type_ids = np.asarray(numbers, int)
        max_type_id = int(type_ids.max()) if type_ids.size else 0
        # Placeholder element labels aligned with numeric type IDs (1..max_type_id)
        elements = [f"type{t}" for t in range(1, max_type_id + 1)]
        working_cell = (lattice.copy(), positions.copy(), list(numbers), elements, mag.copy())
    else:
        working_cell = cell

    # Step 1: Relabel atoms by species and magnetic moments
    cell_new, mapping = relabel_by_species_and_magmoms(working_cell)
    
    # Step 2: Find primitive cell
    cell_new2 = prim_cell(cell_new)
    
    # Step 3: Create Phonopy structure for symmetry analysis
    unitcell = PhonopyAtoms(
        cell=cell_new[0],
        scaled_positions=cell_new[1],
        numbers=cell_new[2]
    )
    
    # Step 4: Get standardization rotation matrix from Phonopy
    ph = Phonopy(unitcell, primitive_matrix="auto")
    U = ph.symmetry.dataset["std_rotation_matrix"]
    
    # Step 5: Restore original labels and apply transformation
    cell_out = relabel_by_species_and_magmoms_back(cell_new2, mapping)
    cell_out[-1] = cell_out[-1] @ U.T
    
    P = np.linalg.inv(ph.primitive_matrix)
    
    if len(cell) == 4:
        del cell_out[-2]
    
    return cell_out, P

def replace_in_list(original_list, a, b):
    """
    Replace all occurrences of value 'a' with value 'b' in a list.
    
    This utility function creates a new list where all instances of value 'a'
    are replaced with value 'b', while keeping all other values unchanged.
    
    Args:
        original_list (list): Input list to modify
        a: Value to be replaced
        b: Replacement value
    
    Returns:
        list: New list with replacements made
    """
    modified_list = [
        b if item == a else item 
        for item in original_list
    ]
    return modified_list
