#!/usr/bin/env python3
"""
SSG (Spin Space Group) Standardization Module

This module provides functions for standardizing magnetic crystal structures
according to their Spin Space Group (SSG) symmetry. It handles the conversion
of magnetic structures to standardized forms for consistent analysis.

Main functionality:
- Standardize magnetic crystal structures to canonical SSG forms
- Handle collinear and coplanar magnetic configurations
- Transform spin orientations to standard reference frames
- Generate Wyckoff positions for magnetic atoms

Author: MOM2SSG Team
"""

# Standard library imports
import sys
import copy
import random

# Third-party library imports
from numpy.linalg import norm, inv, det
import numpy as np
from spglib import get_symmetry_dataset
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

# Local module imports
from .pos2abr import pos2abr
from .spintrans import _find_similarity_matrix
from .find_magprim_unit import find_magprim_unit
from .getcell import prim_cell, std_cell
from .poscar_io import read_poscar, write_poscar
from .SG_isomorphism import find_sg_iso_transform
from .SG_utils import identify_SG_lattice
from .load_ssgdata import load_one_ssg, construct_std_operations
from .find_ssg_operation import findAllOp
from .wyckoff import get_swyckoff
from .small_func import *
from .small_func import _rot_about_axis


# Removed duplicate function _rot_between - using _rotation_between instead

def align_collinear_to_z(vecs, force_positive=False, zero_tol=1e-5):
    """
    Align collinear vectors to the z-axis direction.
    
    This function finds the dominant direction among a set of vectors and
    rotates them so that this direction aligns with the z-axis.
    
    Args:
        vecs (array-like): Nx3 array of vectors to align
        force_positive (bool): If True, ensure z-component is positive
        zero_tol (float): Tolerance for considering vectors as zero
    
    Returns:
        tuple: (aligned_vectors, rotation_matrix, reference_vector)
               - aligned_vectors: Rotated vectors aligned to z-axis
               - rotation_matrix: 3x3 rotation matrix used
               - reference_vector: The reference vector used for alignment
    """
    V = np.asarray(vecs, float)
    norms = np.linalg.norm(V, axis=1)

    # Find non-zero vectors
    nz = np.where(norms > zero_tol)[0]
    if nz.size == 0:
        # All vectors are zero, return identity transformation
        return np.eye(3), V.copy(), None
    
    # Use the vector with largest magnitude as reference
    ref_idx = int(nz[np.argmax(norms[nz])])
    ref = V[ref_idx]

    # Find rotation matrix to align reference vector with z-axis
    R = _rotation_between(ref, np.array([0.0, 0.0, 1.0]))

    # Apply rotation to all vectors
    aligned = (R @ V.T).T
    
    # Optionally force positive z-component
    if force_positive:
        neg = (norms > zero_tol) & (aligned[:, 2] < 0)
        aligned[neg] *= -1.0
    
    # Clean up very small values
    aligned[abs(aligned) < 10*zero_tol] = 0.0
    return aligned

def deal_collinear_spin(vec):
    """
    Process collinear spin vectors by aligning them to z-axis.
    
    This is a wrapper function that aligns collinear magnetic moments
    to the z-axis direction for standardization.
    
    Args:
        vec (array-like): Nx3 array of magnetic moment vectors
    
    Returns:
        array: Aligned magnetic moment vectors
    """
    vec_aligned = align_collinear_to_z(vec)
    return vec_aligned



    
def addA(operations):
    """
    Generate auxiliary atom positions for symmetry analysis.
    
    This function creates a set of auxiliary atoms at positions that
    are related by the symmetry operations. These auxiliary atoms help
    in the standardization process by providing additional reference points.
    
    Args:
        operations (dict): Dictionary containing 'RotC' (rotations) and 'TauC' (translations)
    
    Returns:
        array: Nx3 array of auxiliary atom fractional coordinates
    """
    # Use a fixed seed for reproducible results
    random.seed()    
    A_frac = [random.random(), random.random(), random.random()]
    # Use fixed coordinates for consistency
    A_frac = [0.12345, 0.2715333316, 0.417356188543321]
    A_frac_all = []
    
    # Generate positions for auxiliary atoms under each symmetry operation
    for i in range(len(operations['RotC'])):
        R = operations['RotC'][i]
        tau = operations['TauC'][i]
        # Apply symmetry operation: rotation + translation
        A_frac_per = (R @ A_frac + tau) % 1.0
        A_frac_all.append(A_frac_per)
    
    # Remove duplicate positions
    A_frac_all = dedup_vectors(A_frac_all)
    
    return np.array(A_frac_all)
    

def _rotation_between(a, b, tol=1e-5):
    """
    Find rotation matrix that rotates vector a to vector b.
    
    Uses Rodrigues' rotation formula to compute the rotation matrix
    that transforms vector a to align with vector b.
    
    Args:
        a (array-like): Source vector (3D)
        b (array-like): Target vector (3D)
        tol (float): Tolerance for numerical comparisons
    
    Returns:
        array: 3x3 rotation matrix
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    
    # Handle zero vectors
    if na < tol or nb < tol:
        return np.eye(3)
    
    # Normalize vectors
    a, b = a/na, b/nb
    
    # Compute cross product and dot product
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = float(np.dot(a, b))
    
    # Handle parallel vectors
    if s < tol:
        if c > 0: 
            # Vectors are parallel and in same direction
            return np.eye(3)
        # Vectors are parallel but opposite - need 180° rotation
        tmp = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
        axis = np.cross(a, tmp)
        axis /= np.linalg.norm(axis)
        return _rot_about_axis(axis, np.pi)
    
    # Rodrigues' rotation formula
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * ((1 - c) / (s * s))




def rotate_coplanar_to_xy(vectors, ref=None, tol=1e-6):
    V = np.asarray(vectors, float)
    if V.ndim == 1: V = V[None, :]

    n = _plane_normal_from_vectors(V, tol)
    ez = np.array([0,0,1.0])
    if np.dot(n, ez) < 0: 
        n = -n
    R1 = _rotation_between(n, ez)
    V1 = V @ R1.T
    

    if ref is not None:
        v_xy = (np.asarray(ref, float) @ R1.T).copy()
    else:
        idx = np.argmax([np.linalg.norm(np.array([vv[0], vv[1]])) for vv in V1])
        v_xy = V1[idx].copy()
    v_xy[2] = 0.0
    try:
        angle = np.atan2(v_xy[1], v_xy[0]) if np.linalg.norm(v_xy) >= tol else 0.0
    except:
        angle = np.arctan2(v_xy[1], v_xy[0]) if np.linalg.norm(v_xy) >= tol else 0.0
    Rz = _rot_about_axis(ez, -angle)

    R = Rz @ R1
    V_rot = (vectors @ R.T) if np.asarray(vectors).ndim == 2 else (vectors @ R.T).ravel()
    return V_rot

def judge_raw_translation(operations):
    num_sym = len(operations['spin'])
    for i in range(num_sym):
        
        rot = operations['RotC'][i]
        if not np.allclose(rot,np.eye(3)):
            continue
        
        srot = operations['spin'][i]
        
        if not np.allclose(srot,np.eye(3,dtype=np.float128),atol=1e-3):
            continue
        
        tau = operations['TauC'][i] % 1.0
        
        if not np.allclose(tau,np.array([0.0,0.0,0.0]),atol=1e-3):
            return True
    return False
        

def transform_super_cell_(cell,  M, C, tol=1e-4):

    # C = C.T
    
    lattice, positions, numbers, elements, mag = cell
    M = np.asarray(M, float)
    
    lattice_new = M.T @ lattice

    P = np.asarray(positions, float)   # Fractional coordinates (row vectors)
    P_new = wrap01(P @ C.T)

    return [lattice_new, P_new, numbers, elements, mag]



def standardize_ssg_cell(file_name,ssgnum_in, cell, operations, ssg_list, tol=1e-3,tolm=1e-4,symm_flag=False, check_flag=True):
    """
    Standardize a magnetic crystal structure according to its Spin Space Group (SSG).
    
    This is the main function that transforms a magnetic crystal structure into its
    canonical standardized form according to the specified SSG number. The function
    handles various types of magnetic configurations including collinear and coplanar
    magnetic moments.
    
    Args:
        ssgnum_in (str): Target SSG number (e.g., '1.1', '1.2L', '1.3P')
        cell (tuple): Input crystal structure (lattice, positions, numbers, elements, mag)
        operations (dict): SSG operations dictionary
        ssg_list (list): List of available SSG data
        symm_flag (bool): Whether to output symmetric structure
        check_flag (bool): Whether to perform verification after standardization
    
    Returns:
        tuple: Standardized crystal structure
        
    Raises:
        SystemExit: If standardization fails or verification fails
        
    Note:
        The function performs the following steps:
        1. Check for raw translations and find magnetic primitive unit if needed
        2. Add auxiliary atoms for symmetry analysis
        3. Find isomorphism transformations
        4. Apply transformations to standardize the structure
        5. Handle spin standardization (collinear/coplanar)
        6. Generate Wyckoff positions
        7. Verify the result (if check_flag=True)
    """
    
    # Step 1: Check for raw translations and prepare cell
    flag = judge_raw_translation(operations)
    spg_ssg = ssgnum_in.split('.')[0]  # Extract space group number
    cell = read_poscar(file_name=file_name)
    
    # Find magnetic primitive unit if raw translations are present
    if flag:
        cell,_ = find_magprim_unit(cell)
    
    # Prepare cell for symmetry analysis (remove elements, keep magnetic moments)
    cell_ = cell[:-2] + cell[-1:]
    operations = findAllOp(cell_, tol=tol,tolm=tolm)

    # Standardize the crystal structure
    spg, cell_std = std_cell(cell,tol=tol)

    # if str(spg) != spg_ssg:
    if True:
        A_frac_all = addA(operations)
        
        cell_tmp = copy.deepcopy(list(cell))
        
        cell_tmp[1] = np.concatenate([cell_tmp[1], A_frac_all], axis=0)
        cell_tmp[2] += [cell_tmp[2][-1]+1 for i in range(len(A_frac_all))]
        cell_tmp[3].append('A')
        cell_tmp[4] = np.concatenate([cell_tmp[4], np.zeros([len(A_frac_all),3])], axis=0)
        cell = tuple(cell_tmp)

        # write_poscar(cell, 'POSCAR.addA')

        spg, cell_std = std_cell(cell,tol=tol)
    
    cell_prim = prim_cell(cell_std,tol=tol)
    
    cell_pos2abr,shift = pos2abr(cell_prim)

    unitcell = PhonopyAtoms(
        cell=cell[0],
        scaled_positions=cell[1],
        numbers=cell[2]
    )

    ph = Phonopy(unitcell, primitive_matrix="auto",symprec=tol)
    S1 = inv(ph.primitive_matrix)
    U = ph.symmetry.dataset["std_rotation_matrix"]
    # print(abs(cell[0] - S1.T@cell_prim[0]@U).max())   # == 0
    
    cell2 = list(copy.deepcopy(cell))
    origin_shift= ph.symmetry.dataset["origin_shift"]

    P = ph.symmetry.dataset["transformation_matrix"]

    cell2[1] = wrap01(cell2[1] + origin_shift@inv(P.T) + shift@inv(S1.T))
    cell = copy.deepcopy(cell2)
    

    
    W = find_sg_iso_transform(eval(spg_ssg))['iso_rot']
    isotau = find_sg_iso_transform(eval(spg_ssg))['iso_tau']



    for ssg in ssg_list:
        ssgnum = ssg['ssgNum']
        # if 1:
        if ssgnum == ssgnum_in:
            ssgdic = load_one_ssg(ssgnum)
            superCell = ssgdic['superCell']
            
            break
    
    AClist = []
    for i,j in zip(W,isotau):
        C = inv(superCell)@inv(i)@S1
        

        if np.allclose(det(C),1.0, atol=1e-4):
            if is_integer_matrix(C, tol=1e-4):
                A = i
                for k in j:
                    AClist.append([A,k,C])
    
    if AClist == []:
        print("Can't find A!!! Standardization failed!!!")
        sys.exit()
    
    rots,taus,spins_all = construct_std_operations(ssgdic)
    spg = eval(ssgnum_in.split('.')[0])
    trans_mat_BP_SP = identify_SG_lattice(spg)[1] @ superCell
    l1 = []
    
    for i in range(len(rots)):
        # print(inv((trans_mat_BP_SP)) @ rots[i] @ trans_mat_BP_SP,(inv(trans_mat_BP_SP) @ taus[i]+0.5)%1.0-0.5)
        l1.append([inv((trans_mat_BP_SP)) @ rots[i] @ trans_mat_BP_SP,(inv(trans_mat_BP_SP) @ taus[i]+0.5)%1.0-0.5, np.eye(3)])
        
    cell_bak = copy.deepcopy(cell)
    findP = False
    
    for comb in AClist:
        cell = copy.deepcopy(cell_bak)
        A = comb[0]
        At = comb[1]
        C = comb[2]
        # cell_trans = transform_super_cell(cell, inv(A),C)
        # cart_trans_mat = cell_prim[0].T @ A @ inv(cell_prim[0].T)
        cell_trans = transform_super_cell_(cell, inv(S1)@superCell, C)
        cell_trans[1] = (cell_trans[1] - inv(A@superCell) @ At)%1.0
        # cell_trans[1] = (cell_trans[1] - inv(superCell.T) @ At)%1.0      


        cell_sc_output = copy.deepcopy(list(cell_trans))
        cell_sc_output[0] = cell_sc_output[0] @ U.T
        cell_sc_output[-1] = cell_sc_output[-1] @ U.T
        
        if 'A' in cell_sc_output[-2]:
            cell_sc_output[-2].remove('A')
            A_id = cell_sc_output[2][-1]
            A_num = 0
            
            while cell_sc_output[2][-1] == A_id:
                A_num += 1
                del cell_sc_output[2][-1]
                
            cell_sc_output[1] = cell_sc_output[1][:-A_num]
            cell_sc_output[-1] = cell_sc_output[-1][:-A_num]
        
        write_poscar(cell_sc_output, 'SPOSCAR.spinnorot')

        #----------------------------------------------------------
        ###spin
        if ssgnum_in[-1] == 'L':
            flagL = True
        else:
            flagL = False
        
        flagP = False
        
        if ssgnum_in[-1] == 'P':
            try:
                cell_sc_ = copy.deepcopy(list(cell_sc_output))

                cell_sc_[-1] = rotate_coplanar_to_xy(cell_sc_[-1])
                
                cell_sc_output = copy.deepcopy(cell_sc_)
                flagP = True

            except:
                flagL = True
            

        if flagL:
            cell_sc_ = copy.deepcopy(list(cell_sc_output))
            cell_sc_[-1] = deal_collinear_spin(cell_sc_[-1])
                    
            cell_sc_output = copy.deepcopy(cell_sc_)
        
        
        if not flagL:
            cell_sc_ = cell_sc_output[:-2] + cell_sc_output[-1:]
            operations = findAllOp(cell_sc_, tol = tol,tolm=tolm)
            
            if operations['QLabel'] == '1' or operations['QLabel'] == '-1':
                P = np.eye(3)
                findP = True
                
            else:
                
                    
                l2 = []
                
                for j in range(len(operations['RotC'])):
                    # print([operations['RotC'][j],(operations['TauC'][j]+0.5)%1.0-0.5])
                    l2.append([operations['RotC'][j],(operations['TauC'][j]+0.5)%1.0-0.5,np.eye(3)])
                    
                order = l1_order_in_l2(l1, l2, atol=1e-4, rtol=1e-4)
                
                if order == False:
                    # print(compare_lists_with_report(l1,l2, atol=1e-4, rtol=1e-4))
                    continue

                spin_new = [copy.deepcopy(operations['spin'][j]) for j in order]
                spin2 = dedup_real_3x3(spin_new, tol=1e-3)

                
                if flagP and operations['QLabel'] == 'm':
                    spin1 = [np.eye(3), np.diag([-1,1,1])]
                    P = _find_similarity_matrix(operations['QLabel'], spin1, spin2, ordered=True)
                    findP = True
                else:
                    for spin_id,spins in enumerate(spins_all):
                        if spins == []:
                            continue
                        spin1 = dedup_real_3x3(spins, tol=1e-3)

                        # print_matlist(spin1)
                
                        # print(operations['QLabel'])
                        
                        P = _find_similarity_matrix(operations['QLabel'], spin1, spin2, ordered=True)
                        if P is not None:
                            findP = True
                            break

            if not findP:
                continue
            else:
                for i in range(len(cell_sc_output[-1])):
                    cell_sc_output[-1][i] = P.T@cell_sc_output[-1][i]
                
                write_poscar(cell_sc_output, 'SPOSCAR')
                break
        else:
            
            write_poscar(cell_sc_output, 'SPOSCAR')
            break
    
        

    if (not findP) and (not flagL):
        print('Spin standardization failed!!! P not found!')
        sys.exit()

    cell_sc_ = cell_sc_output[:-2] + cell_sc_output[-1:]
    operations = findAllOp(cell_sc_, tol = tol,tolm=tolm)
    cell_sc_output_acc,wyckoff = get_swyckoff(cell_sc_output,operations)
    # output_wyckoff(cell_sc_output_acc, wyckoff, filename='swyckoff.out')
    
    if symm_flag:
        write_poscar(cell_sc_output_acc, 'SPOSCAR.symm')
    
    
    ###check
    
    from .MOM2SSG import get_ssg
    if check_flag:
    
        ssg_list, ssgnum_out,operations = get_ssg('SPOSCAR',ssg_list=ssg_list,tol=tol,tolm=tolm)
        if ssgnum_out != ssgnum_in:
            print('Warning: the identified SSG number is different from the input one!!! Please check carefully!')
            sys.exit()
        
        else:
            try:
                spin_id
            except:
                spin_id = 0
                
            l1 = []
            
            for i in range(len(rots)):
                if ssgnum_out[-1] == 'L':
                    l1.append([inv((trans_mat_BP_SP)) @ rots[i] @ trans_mat_BP_SP,(inv(trans_mat_BP_SP) @ taus[i]+0.5)%1.0-0.5, det(spins_all[0][i])*np.eye(3)])
                elif ssgnum_out[-1] == 'P' and operations['QLabel'] == 'm':
                    if det(spins_all[0][i]) > 0:
                        l1.append([inv((trans_mat_BP_SP)) @ rots[i] @ trans_mat_BP_SP,(inv(trans_mat_BP_SP) @ taus[i]+0.5)%1.0-0.5, np.eye(3)])
                    else:
                        l1.append([inv((trans_mat_BP_SP)) @ rots[i] @ trans_mat_BP_SP,(inv(trans_mat_BP_SP) @ taus[i]+0.5)%1.0-0.5, np.diag([-1,1,1])])
                else:
                    l1.append([inv((trans_mat_BP_SP)) @ rots[i] @ trans_mat_BP_SP,(inv(trans_mat_BP_SP) @ taus[i]+0.5)%1.0-0.5, spins_all[spin_id][i]])


            l2 = []

            for j in range(len(operations['RotC'])):
                if ssgnum_out[-1] == 'L':
                    l2.append([operations['RotC'][j],(operations['TauC'][j]+0.5)%1.0-0.5, det(operations['spin'][j])*np.eye(3)])
                else:
                    l2.append([operations['RotC'][j],(operations['TauC'][j]+0.5)%1.0-0.5, operations['spin'][j]])

            order = l1_order_in_l2(l1, l2, atol=1e-4, rtol=1e-4)
            if order == False:
                print('Spin standardize failed!!!')
    
