from irssg import findAllOp,identify_SG_lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from scipy.linalg import null_space

from numpy.linalg import norm, inv
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element

from .mag_atom import ELEMENT_TABLE,analyze_magnetic_atoms
from .poscar_io import write_poscar

# Helper function to identify different Wyckoff positions using space group symmetry
def identify_wyckoff_positions(positions, atom_type, cell, tol=1e-3):
    """
    Identify different Wyckoff positions for atoms of the same type using space group symmetry.
    
    Args:
        positions: List of atomic positions for the specific atom type
        atom_type: Atomic type number
        cell: Complete cell tuple (lattice, positions, numbers, mag)
        tol: Tolerance for symmetry operations
        
    Returns:
        dict: {wyckoff_id: [positions]} mapping Wyckoff position ID to list of positions
    """
    lattice, all_positions, numbers, mag = cell
    
    # Create a temporary structure with only the atoms of the specific type
    temp_positions = []
    temp_numbers = []
    for i, pos in enumerate(all_positions):
        if numbers[i] == atom_type:
            temp_positions.append(pos)
            temp_numbers.append(atom_type)
    
    if not temp_positions:
        return {}
    
    
    # Use dummy element for the specific atom type
    structure = Structure(lattice, [Element.from_Z(atom_type)] * len(temp_positions), temp_positions)
    
    # Get space group symmetry information
    sga = SpacegroupAnalyzer(structure, symprec=tol)
    sym_data = sga.get_symmetry_dataset()
    
    # Get Wyckoff positions and equivalent atoms
    wyckoff_letters = sym_data['wyckoffs']
    equivalent_atoms = sym_data['equivalent_atoms']
    
    # Group atoms by their Wyckoff position
    wyckoff_groups = {}
    wyckoff_id = 0
    
    # Create mapping from original position indices to temp structure indices
    pos_to_temp_idx = {}
    temp_idx = 0
    for i, pos in enumerate(all_positions):
        if numbers[i] == atom_type:
            pos_to_temp_idx[tuple(pos)] = temp_idx
            temp_idx += 1
    
    # Group atoms by Wyckoff equivalence
    used_positions = set()
    for i, pos in enumerate(positions):
        if tuple(pos) in used_positions:
            continue
            
        # Find the corresponding index in temp structure
        if tuple(pos) not in pos_to_temp_idx:
            continue
            
        temp_idx = pos_to_temp_idx[tuple(pos)]
        equiv_atom = equivalent_atoms[temp_idx]
        
        # Find all atoms equivalent to this one
        wyckoff_group = [pos]
        used_positions.add(tuple(pos))
        
        # Find other atoms in the same Wyckoff position
        for j, other_pos in enumerate(positions):
            if j != i and tuple(other_pos) not in used_positions:
                if tuple(other_pos) in pos_to_temp_idx:
                    other_temp_idx = pos_to_temp_idx[tuple(other_pos)]
                    if equivalent_atoms[other_temp_idx] == equiv_atom:
                        wyckoff_group.append(other_pos)
                        used_positions.add(tuple(other_pos))
        
        wyckoff_groups[wyckoff_id] = wyckoff_group
        wyckoff_id += 1
    
    return wyckoff_groups
    
def findDimension(cell): # collinear or coplanar or no-coplanar, return 3, 2, 1
    mag = cell[-1].copy()
    rank = np.linalg.matrix_rank(mag)
    return rank



# ************************ function ssg2magmom **************************
# 
# from the operations of a SSG to generate MAGMOM
# check the new MAGMOM satisfy the original SSG
# SUPPORTS MULTIPLE MAGNETIC ATOM TYPES
#
# ***********************************************************************
def ssg2magmom(ssg_dic, cell, magmom_dict,tol=1e-3):
    lattice, positions, numbers, mag = cell
    dim_mag = ssg_dic['dim']
    mag_info = analyze_magnetic_atoms(cell,magmom_dict)
    
    # Support multiple magnetic atom types
    magnetic_species = mag_info['magnetic_species_moments']  # {type_id: mag_moment}
    magnetic_indices = mag_info['magnetic_index']  # list of magnetic atom type IDs
    
    print(f"Found {len(magnetic_indices)} types of magnetic atoms.")
    
    # Group magnetic atoms by type
    mag_atoms_by_type = {}
    for type_id in magnetic_indices:
        mag_atoms_by_type[type_id] = {
            'positions': [positions[i] for i, num in enumerate(numbers) if num == type_id],
            'indices': [i for i, num in enumerate(numbers) if num == type_id],
            'moment': magnetic_species[type_id]
        }
        print(f"{ELEMENT_TABLE[type_id]}: {len(mag_atoms_by_type[type_id]['positions'])} atoms, moment = {magnetic_species[type_id]}")
    
    # Use the first magnetic atom type as reference for symmetry analysis
    ref_type_id = magnetic_indices[0]
    positions_mag = mag_atoms_by_type[ref_type_id]['positions']
    # print(positions_mag)
    # load SSG
    ssgNum = ssg_dic['ssgNum']
    HRot = ssg_dic['HRotC']
    HTau = ssg_dic['HTauC']
    QRot = ssg_dic['QRotC']
    QTau = ssg_dic['QTauC']
    # UROt can have equivalent SSG , here we choose the 0
    URot = ssg_dic['URot'][0]
    rotlist = []
    taulist = []
    spinlist = []
    # get all the operation under the operation of the translation generators
    gid = int(ssgNum.split('.')[0])
    prim_vec = identify_SG_lattice(gid)[1] # each col is a prim basis vector   
    tau_col = ssg_dic['superCell']
    t1 = np.array([tau_col[0][0], tau_col[1][0], tau_col[2][0]]) 
    t2 = np.array([tau_col[0][1], tau_col[1][1], tau_col[2][1]]) 
    t3 = np.array([tau_col[0][2], tau_col[1][2], tau_col[2][2]])
    pure_t = []
    for t_append in [t1, t2, t3]:
        pure_t.append(prim_vec @ t_append)
    # print(pure_t)
    P = np.array(pure_t).T
    # print(P)
    for i, hrot in enumerate(HRot):
        htau = HTau[i]
        hrot1 = inv(P) @ hrot @ P
        htau1 = inv(P) @ htau
        for j, spin in enumerate(URot):
            qrot = QRot[j]
            qtau = QTau[j]
            qrot1 = inv(P) @ qrot @ P
            qtau1 = inv(P) @ qtau

            # r1 t1 * r2 t2 = r1r2| r1t2 + t1
            # qrot|qtau * hrot|htau = qrot*hrot| qrot*htau + qtau
            rot_new = qrot1 @ hrot1
            tau_new = qrot1 @ htau1 + qtau1
            # rotation matrix under the primitive(supercell) lattice, i.e. lattice_s
            # rot and tau in the dictionary are written under conventional cell
            # rot_P = pure_t @ rot_new @ inv(pure_t)
            # tau_P = tau_new @ inv(pure_t)
            # print(tau_new)
            # print(tau_P)
            # print('======')
            rotlist.append(rot_new)
            taulist.append(tau_new)
            spinlist.append(spin)
    # find the operations that being the onsite symmetry
    spin_onsite = []
    for iop, irot in enumerate(rotlist):
        itau = taulist[iop]
        pos_onsite = positions_mag[0]
        new_pos = irot @ pos_onsite
        s1 = (new_pos[0]+6 + itau[0]) % 1
        s2 = (new_pos[1]+6 + itau[1]) % 1
        s3 = (new_pos[2]+6 + itau[2]) % 1
        if norm(pos_onsite - np.array([s1, s2, s3])) < 1e-3:
            # print('=====================')
            # print(irot)
            # print(itau)
            spin_onsite.append(spinlist[iop])
    # unique the spin onsite
    onsite = [np.eye(3)]
    for mat in spin_onsite:
        if not any(np.linalg.norm(mat - existing) < 1e-3 for existing in onsite):
            onsite.append(mat)
    mat_zeros = []
    for mat in onsite:
        mat_zeros.append(mat - np.eye(3))
    if dim_mag == 2:
        mat_zeros.append(np.diag((1,1,-1)) - np.eye(3))
    if dim_mag == 1:
        mat_zeros.append(np.diag((-1,-1,1)) - np.eye(3))

    # find the null_space
    stack_mat = np.vstack(mat_zeros)
    common_null_space = null_space(stack_mat)
    # print(len(spinlist))
    # print(stack_mat)
    # print(common_null_space)
    degree_freedom = np.size(common_null_space, 1)
    # print('degree of freedom of onsite symmetry:', degree_freedom)
    if degree_freedom == 0:
        print('False because not permitted magnetic moments for magnetic atoms.')
        return False,[],[]
    # Generate magnetic moments for each atom type
    coefficient = [3.1, 4.3, 5.1]
    
    # Generate reference magnetic moment for the first atom type
    magmom_ref = np.array([0.,0.,0.])
    for i in range(degree_freedom):
        al = coefficient[i] * common_null_space[:, i]
        magmom_ref += al
    magmom_ref = mag_atoms_by_type[ref_type_id]['moment'] * magmom_ref * (1/norm(magmom_ref))
    
    # Helper function to check if position is already in list
    def notInpos(s1, s2, s3, mag_pos):
        for pos in mag_pos:
            if norm(np.array([s1, s2, s3]) - pos) < 1e-3:
                return False
        return True
    
    # Generate magnetic moments for all atom types
    # For single magnetic atom type, use Wyckoff-based logic like multi-mag
    if len(magnetic_indices) == 1:
        type_id = magnetic_indices[0]
        # print(f"Using Wyckoff-based single magnetic atom logic for type {type_id}...")
        
        # Get Wyckoff positions for this atom type
        wyckoff_groups = identify_wyckoff_positions(
            mag_atoms_by_type[type_id]['positions'], 
            type_id, 
            cell, 
            tol
        )
        # print(f"Found {len(wyckoff_groups)} Wyckoff positions for type {type_id}")
        
        magmom_dict = {type_id: {}}
        mag_pos_dict = {type_id: {}}
        
        for wyckoff_id, wyckoff_positions in wyckoff_groups.items():
            # print(f"  Processing Wyckoff position {wyckoff_id} with {len(wyckoff_positions)} atoms...")
            
            # Get reference position for this Wyckoff position
            ref_pos = wyckoff_positions[0]
            ref_moment = mag_atoms_by_type[type_id]['moment']
            
            # Scale the reference moment by the ratio of magnetic moments
            moment_scale = ref_moment / mag_atoms_by_type[ref_type_id]['moment']
            magmom1 = moment_scale * magmom_ref
            
            magmom_list = []
            mag_pos_list = []
            
            # Apply symmetry operations to generate magnetic moments
            for iop, irot in enumerate(rotlist):
                itau = taulist[iop]
                new_pos = irot @ ref_pos
                s1 = (new_pos[0]+6 + itau[0]) % 1
                s2 = (new_pos[1]+6 + itau[1]) % 1
                s3 = (new_pos[2]+6 + itau[2]) % 1
                
                mag_pos_list.append(np.array([s1, s2, s3]))
                magmom_list.append(spinlist[iop] @ magmom1)
            
            # Assign magnetic moments to actual atom positions
            for atom_pos in wyckoff_positions:
                min_dist = float('inf')
                best_magmom = magmom1  # Default to reference moment
                
                # Find the best matching symmetry-generated position
                for i, sym_pos in enumerate(mag_pos_list):
                    diff = np.array(atom_pos) - np.array(sym_pos)
                    diff = diff - np.round(diff)  # Apply periodic boundary conditions
                    dist = norm(diff)
                    
                    if dist < min_dist:
                        min_dist = dist
                        best_magmom = magmom_list[i]
                
                # If no good match found, use reference moment
                if min_dist >= tol:
                    best_magmom = magmom1
                
                mag_pos_list.append(np.array(atom_pos))
                magmom_list.append(best_magmom)
            
            magmom_dict[type_id][wyckoff_id] = magmom_list
            mag_pos_dict[type_id][wyckoff_id] = mag_pos_list
            # print(f"  Generated {len(magmom_list)} magnetic moments for Wyckoff position {wyckoff_id}")
        
        # print(f"Final: {sum(len(mom_list) for mom_list in magmom_dict[type_id].values())} magnetic moments for {len(mag_atoms_by_type[type_id]['positions'])} atoms of type {type_id}")
        
    else:
        # Use multi-magnetic atom logic for multiple atom types
        magmom_dict = {}  # {type_id: {wyckoff_id: [magmom_list]}}
        mag_pos_dict = {}  # {type_id: {wyckoff_id: [pos_list]}}
        
        for type_id in magnetic_indices:
            # print(f"Generating magnetic moments for atom type {type_id}...")
            
            # Identify different Wyckoff positions for this atom type using space group symmetry
            wyckoff_groups = identify_wyckoff_positions(
                mag_atoms_by_type[type_id]['positions'], 
                type_id, 
                cell, 
                tol
            )
            # print(f"Found {len(wyckoff_groups)} Wyckoff positions for type {type_id}")
            
            magmom_dict[type_id] = {}
            mag_pos_dict[type_id] = {}
            
            for wyckoff_id, wyckoff_positions in wyckoff_groups.items():
                # print(f"  Processing Wyckoff position {wyckoff_id} with {len(wyckoff_positions)} atoms...")
                
                # Get reference position for this Wyckoff position
                ref_pos = wyckoff_positions[0]
                ref_moment = mag_atoms_by_type[type_id]['moment']
                
                # Scale the reference moment by the ratio of magnetic moments
                moment_scale = ref_moment / mag_atoms_by_type[ref_type_id]['moment']
                magmom1 = moment_scale * magmom_ref
                
                magmom_list = [magmom1]
                mag_pos_list = [ref_pos]
                
                # Apply symmetry operations to generate all magnetic moments for this Wyckoff position
                for iop, irot in enumerate(rotlist):
                    itau = taulist[iop]
                    new_pos = irot @ ref_pos
                    s1 = (new_pos[0]+6 + itau[0]) % 1
                    s2 = (new_pos[1]+6 + itau[1]) % 1
                    s3 = (new_pos[2]+6 + itau[2]) % 1
                    
                    if notInpos(s1, s2, s3, mag_pos_list):
                        mag_pos_list.append(np.array([s1, s2, s3]))
                        magmom_list.append(spinlist[iop] @ magmom1)
                
                magmom_dict[type_id][wyckoff_id] = magmom_list
                mag_pos_dict[type_id][wyckoff_id] = mag_pos_list
                # print(f"  Generated {len(magmom_list)} magnetic moments for Wyckoff position {wyckoff_id}")
    
    # Assign magnetic moments to all atoms in the cell
    for i, pos in enumerate(positions):
        atom_type = numbers[i]
        if atom_type in magnetic_indices:
            # Find the closest position in the generated list across all Wyckoff positions
            min_dist = float('inf')
            best_magmom = np.array([0., 0., 0.])
            
            for wyckoff_id, wyckoff_positions in mag_pos_dict[atom_type].items():
                for j, pos1 in enumerate(wyckoff_positions):
                    # Calculate distance considering periodic boundary conditions
                    diff = np.array(pos) - np.array(pos1)
                    diff = diff - np.round(diff)  # Apply periodic boundary conditions
                    dist = norm(diff)
                    if dist < min_dist:
                        min_dist = dist
                        best_magmom = magmom_dict[atom_type][wyckoff_id][j]
            
            if min_dist < tol:
                mag[i] = best_magmom
            else:
                print(f"Warning: Could not find matching position for atom {i} of type {atom_type}")
                # Use the reference moment as fallback
                mag[i] = magmom_dict[atom_type][0][0]  # Use first available moment
    
    cell_new = (lattice, positions, numbers, mag)
    # Debug output (can be enabled by setting debug=True)
    debug = False  # Set to True for debugging
    if debug:
        print('cell_new:')
        print(cell_new)
        write_poscar(cell_new, file_name = 'POSCAR.test')    # print(cell_new)
    # write_poscar(cell_new, file_name = 'output/'+ ssgNum)
    (lattice, positions, numbers, mag) = cell_new
    # print('numbers:', numbers)
    # numbers = [1,1,1,1,2,2,3,3,3,3,3,3,3,3,3,3,3,3]
    # cell_new = (lattice, positions, numbers, mag)
    # Save original cell_new before findAllOp modifies it
    cell_original = (cell_new[0].copy(), cell_new[1].copy(), cell_new[2].copy(), cell_new[3].copy())
    
    ops_dic = findAllOp(cell_new, tol = tol, tolm=1e-4)
    # ops_dic must maintain all the ops in ssgNum, since the magmom is construct using this ops
    # so we only need to check it has the same number of operation
    # and check it has the same dimension (L, P, NP)
    # print(ops_dic)
    dim_cell = findDimension(cell_new)
    
    if not dim_cell == dim_mag:
        print('Construct the SSG cell failed because the result dimension in not equal to the SSG.')
        return False,[],ops_dic
    if not len(ops_dic['HRotC']) == len(HRot):
        # print('False because H.')
        print('Construct the SSG cell failed.')
        return False,[],ops_dic
    if not len(ops_dic['RotC']) == len(rotlist):
        # print('False because RotC.')
        print('Construct the SSG cell failed.')
        # Even if validation fails, write the POSCAR file
        # write_poscar(cell_original, file_name = 'SPOSCAR')
        return False,[],ops_dic
    Gnum_ssg = int(ssgNum.split('.')[0])
    if not Gnum_ssg == int(ops_dic['Gnum']):
        # print('False because G number.')
        print('Construct the SSG cell failed.')
        # Even if validation fails, write the POSCAR file
        # write_poscar(cell_original, file_name = 'SPOSCAR')
        return False,[],ops_dic
    
    # for write POSCAR file - use original cell to preserve magnetic moments
    
    return True,cell_original,ops_dic