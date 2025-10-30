"""
MOM2SSG: Magnetic Order to Spin Space Group Analysis

This module provides comprehensive analysis tools for determining spin space groups
from magnetic crystal structures. It can identify magnetic configurations, analyze
symmetry operations, and determine the corresponding spin space group labels.

Main functionality:
- Read magnetic crystal structures from POSCAR or mcif files
- Identify non-magnetic space groups and magnetic point groups
- Determine spin space group numbers and labels
- Analyze supercell relationships and volume ratios
- Generate symmetry operations and Wyckoff positions

Author: MOM2SSG Development Team
"""

# Standard library imports
import os
import pickle
import warnings

from math import sin, pi
import sys
import argparse

# Third-party imports
import numpy as np
from numpy.linalg import norm, det

import copy
from spglib import get_symmetry_dataset

# Local module imports
from .small_func import *
from .SG_utils import sg_symbol_from_number
from .SSGLabel import get_SSG_label
from .lirssg import generate_irssg_in
from .poscar_io import read_poscar_no_elements, mcif2cell, write_poscar
from .load_ssgdata import load_ssg_list
from .wyckoff import get_swyckoff, output_wyckoff
from .find_ssg_operation import findAllOp
from .get_ssg_number import search4ssg
from .std2ssg import standardize_ssg_cell,addA
from .spintrans import spin_axis
from .eqvpg2label import get_std_pg
from .find_magprim_unit import find_magprim_unit


def format_output(dim_mag,axis_vector,spin_rot_list,operations,lps,pg_op_num,nonmag_sym,ssgnum,format_ssg,cell,tol=1e-3):
    num_operator = len(operations['spin'])
    
    space_international = nonmag_sym['international']
    spg_op_num = nonmag_sym['rotations'].shape[0]
    ncell_pos_ssg = count_identity_pairs(operations['RotC'], operations['spin'], tol=1e-3) #find pure translation number
    ncell_pos_asg = count_identity_pairs(operations['RotC'], [np.eye(3)]*len(operations['RotC']), tol=1e-3)
    assert ncell_pos_asg%ncell_pos_ssg == 0
    ncell_ssg_asg = int(ncell_pos_asg)//int(ncell_pos_ssg)
    assert ncell_ssg_asg%operations['Ik'] == 0
    assert spg_op_num%num_operator == 0
    
    if ssgnum == 'need more loop':
        print('The SSG number can not be identified!')
    else:

        print(f'The SSG number: {ssgnum}')
        print(f'The SSG label: {format_ssg}')
        
    if 'Collinear' in lps:
        print('''I: Collinear SSG along z' direction; So: C\N{INFINITY}={C\N{INFINITY}z',Mx'C\N{INFINITY}z'}''')
    elif 'Coplanar' in lps:
        print('''II: Coplanar SSG in the x'y' plane; So: Cs={E,Mz'}''')
    elif 'Noncoplanar' in lps:
        print('''III: Noncoplanar SSG; So: C1={E}''')
        
    
    print("The SSG G = So x Go")
    print('P (spin part of Go): ' + get_std_pg(operations['QLabel'])[1])
    #print('H (lattice part of Go): '+ sg_symbol_from_number(operations['Gnum']) + f' ({num_operator//ncell_pos_ssg} operations)') #wzj
    print('H (lattice part of Go): '+ sg_symbol_from_number(operations['Gnum']))
    
    # print(f"The volume of POSCAR is {ncell_pos_ssg} times of the SSG cell.")

    # print(f"The volume of SSG cell is {ncell_ssg_asg} times of ASG cell.")
    
    # print(f"The volume of SSG_SPG cell is {ncell_ssg_asg//operations['Ik']} times of ASG cell.")
    # print(f"The volume of SSG cell is {operations['Ik']} times of SSG_SPG cell.")
    
    print('='*40)
    
    print('Spin space group operations: {U||R|v}')
    print("U is given in x'y'z' coordinates, while R is given in the lattice basis of POSCAR.")
    
    if 'Collinear' in lps:
        print('''U= I_2 + \u03BE in type I ''')
    elif 'Coplanar' in lps:
        print('''U= \u03B6_{2x2} + I_1 in type II''')
        
    print(f'# Number: {num_operator}')
    # if dim_mag == 1:
    #     print('{Spin|| Ri  | taui}')
    # elif dim_mag == 2:
    #     print('{   Spin O2 || Ri  | taui}')
    # else:
    #     print('{      Spin O3    || Ri  | taui}')

    if dim_mag == 1:
        print('{  \u03BE ||   Ri   |  taui }')
    elif dim_mag == 2:
        print('{   \u03B6_{2x2}   ||    Ri   |  taui }')
    else:
        print('{          U         ||   Ri    |  taui }')
        
    for i in range(num_operator):
        time_reversal = det(operations['spin'][i])
        if time_reversal > 0:
            time_reversal_str = 'without time-reversal'
        else:
            time_reversal_str = 'with time-reversal'
            
        print(f'# {i+1:>3d}   {time_reversal_str}')
        
        operations["TauC"][i] = (operations["TauC"][i]+0.5)%1-0.5
        if dim_mag == 1:
            print(f"{spin_rot_list[i]:>4d}",end='  ')
            print("".join(f"{v:>2d} " for v in np.asarray(operations["RotC"][i][0, :].reshape(-1), int)),end='  ')
            print("".join(f"{v:>6.3f} " for v in np.asarray(operations["TauC"][i][0].reshape(-1), float)+1e-6))
            
            print("      ",end='')
            print("".join(f"{v:>2d} " for v in np.asarray(operations["RotC"][i][1, :].reshape(-1), int)),end='  ')
            print("".join(f"{v:>6.3f} " for v in np.asarray(operations["TauC"][i][1].reshape(-1), float)+1e-6))
            
            print("      ",end='')
            print("".join(f"{v:>2d} " for v in np.asarray(operations["RotC"][i][2 :].reshape(-1), int)),end='  ')
            print("".join(f"{v:>6.3f} " for v in np.asarray(operations["TauC"][i][2].reshape(-1), float)+1e-6))


        elif dim_mag == 2:
            print("".join(f"{v:>6.3f} " for v in np.asarray(spin_rot_list[i][0, :], float)+1e-6),end='  ')
            print("".join(f"{v:>2d} " for v in np.asarray(operations["RotC"][i][0, :].reshape(-1), int)),end='  ')
            print("".join(f"{v:>6.3f} " for v in np.asarray(operations["TauC"][i][0].reshape(-1), float)+1e-6))
            
            
            print("".join(f"{v:>6.3f} " for v in np.asarray(spin_rot_list[i][1, :], float)+1e-6),end='  ')
            print("".join(f"{v:>2d} " for v in np.asarray(operations["RotC"][i][1, :].reshape(-1), int)),end='  ')
            print("".join(f"{v:>6.3f} " for v in np.asarray(operations["TauC"][i][1].reshape(-1), float)+1e-6))
            
            print(" "*16,end='')
            print("".join(f"{v:>2d} " for v in np.asarray(operations["RotC"][i][2 :].reshape(-1), int)),end='  ')
            print("".join(f"{v:>6.3f} " for v in np.asarray(operations["TauC"][i][2].reshape(-1), float)+1e-6))


        else:
            print("".join(f"{v:>6.3f} " for v in np.asarray(spin_rot_list[i][0, :], float)+1e-6),end='  ')
            print("".join(f"{v:>2d} " for v in np.asarray(operations["RotC"][i][0, :].reshape(-1), int)),end='  ')
            print("".join(f"{v:>6.3f} " for v in np.asarray(operations["TauC"][i][0].reshape(-1), float)+1e-6))
            
            
            print("".join(f"{v:>6.3f} " for v in np.asarray(spin_rot_list[i][1, :], float)+1e-6),end='  ')
            print("".join(f"{v:>2d} " for v in np.asarray(operations["RotC"][i][1, :].reshape(-1), int)),end='  ')
            print("".join(f"{v:>6.3f} " for v in np.asarray(operations["TauC"][i][1].reshape(-1), float)+1e-6))
            
            print("".join(f"{v:>6.3f} " for v in np.asarray(spin_rot_list[i][2, :], float)+1e-6),end='  ')
            print("".join(f"{v:>2d} " for v in np.asarray(operations["RotC"][i][2 :].reshape(-1), int)),end='  ')
            print("".join(f"{v:>6.3f} " for v in np.asarray(operations["TauC"][i][2].reshape(-1), float)+1e-6))
    
    print()
    print("The redefined axises in spin space: (x',y',z')=(x,y,z)D")
    
    # print(f"[{axis_vector[1][0]:>6.3f}, {axis_vector[1][1]:>6.3f}, {axis_vector[1][2]:>6.3f}]:x'")
    # print(f"[{axis_vector[2][0]:>6.3f}, {axis_vector[2][1]:>6.3f}, {axis_vector[2][2]:>6.3f}]:y'")
    # print(f"[{axis_vector[0][0]:>6.3f}, {axis_vector[0][1]:>6.3f}, {axis_vector[0][2]:>6.3f}]:z'")
    
    A = [axis_vector[1], axis_vector[2], axis_vector[0]]
    rows = [f"{r0+1e-6:>6.3f} {r1+1e-6:>6.3f} {r2+1e-6:>6.3f}" for r0, r1, r2 in zip(*A)]
    prefix = "D = "
    print(prefix + rows[0])
    indent = " " * len(prefix)
    for r in rows[1:]:
        print(indent + r)
    
    print("D U D^{-1}  is in Cartesian coordinates, which is in ssg.data.")
    
    print('='*40)

    print(f"Atomic space group: {nonmag_sym['number']}, {space_international}")

    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    
    
    unitcell_pos = PhonopyAtoms(
        cell=cell[0],
        scaled_positions=cell[1],
        numbers=cell[2]
    )
    
    ph = Phonopy(unitcell_pos, primitive_matrix="auto",symprec=tol)
    R33 = inv(ph.primitive_matrix)
    det33 = np.linalg.det(R33)
    
    asg_cell = R33.T @ cell[0]
    print('Atomic primitive cell')
    print(f":a  [{asg_cell[0,0]+1e-6:>6.3f}, {asg_cell[0,1]+1e-6:>6.3f}, {asg_cell[0,2]+1e-6:>6.3f}]:")
    print(f":b  [{asg_cell[1,0]+1e-6:>6.3f}, {asg_cell[1,1]+1e-6:>6.3f}, {asg_cell[1,2]+1e-6:>6.3f}]:")
    print(f":c  [{asg_cell[2,0]+1e-6:>6.3f}, {asg_cell[2,1]+1e-6:>6.3f}, {asg_cell[2,2]+1e-6:>6.3f}]:")
    print()
    
    if ncell_pos_ssg > 1:  # magnetic super cell
        cell_mag_unit, P = find_magprim_unit(cell)
    else:
        cell_mag_unit = copy.deepcopy(cell)
        P = np.eye(3)
    
    unitcell_mag = PhonopyAtoms(
        cell=cell_mag_unit[0],
        scaled_positions=cell_mag_unit[1],
        numbers=cell_mag_unit[2]
    )
    ph = Phonopy(unitcell_mag, primitive_matrix="auto",symprec=tol)
    R32 = inv(ph.primitive_matrix)
    det32 = np.linalg.det(R32)
    
    if operations['Gnum'] == nonmag_sym['number']:
        P2 = np.eye(3)
    else:
        
        A_frac_all = addA(operations)
        cell_tmp = copy.deepcopy(list(cell_mag_unit))
        cell_tmp[1] = np.concatenate([cell_tmp[1], A_frac_all], axis=0)
        cell_tmp[2] += [cell_tmp[2][-1]+1 for i in range(len(A_frac_all))]
        cell_tmp[3] = np.concatenate([cell_tmp[3], np.zeros([len(A_frac_all),3])], axis=0)
        unitcell_ssg_spg_A = PhonopyAtoms(
            cell=cell_tmp[0],
            scaled_positions=cell_tmp[1],
            numbers=cell_tmp[2]
        )
        ph = Phonopy(unitcell_ssg_spg_A, primitive_matrix="auto",symprec=tol)
        P2 = inv(ph.primitive_matrix)
        
    R31 = R32 @ inv(P2)

    det31 = np.linalg.det(R31)
    
    
    # Header line for R1/R2/R3 with det values
    hdr1 = f"H:(a1,b1,c1) = (a,b,c)R1, det(R1) = {round(det31)}"
    hdr2 = f"T0:(a2,b2,c2) = (a,b,c)R2, det(R2) = {round(det32)}"
    hdr3 = f"POSCAR:(a3,b3,c3) = (a,b,c)R3, det(R3) = {round(det33)}"
    sep_hdr = "   "
    print(hdr1 + sep_hdr + hdr2 + sep_hdr + hdr3)
    
    # Print R1, R2, R3 side-by-side in three lines, aligned under their headers
    rows_R1 = [f"{r[0]+1e-6:>6.3f} {r[1]+1e-6:>6.3f} {r[2]+1e-6:>6.3f}" for r in R31]
    rows_R2 = [f"{r[0]+1e-6:>6.3f} {r[1]+1e-6:>6.3f} {r[2]+1e-6:>6.3f}" for r in R32]
    rows_R3 = [f"{r[0]+1e-6:>6.3f} {r[1]+1e-6:>6.3f} {r[2]+1e-6:>6.3f}" for r in R33]
    p1, p2, p3 = "R1 = ", "R2 = ", "R3 = "
    i1, i2, i3 = " " * len(p1), " " * len(p2), " " * len(p3)
    # Compute start columns for R2 and R3 to align with header segments
    pos_R2 = len(hdr1) + len(sep_hdr)
    pos_R3 = pos_R2 + len(hdr2) + len(sep_hdr)
    # First row with prefixes
    line0 = p1 + rows_R1[0]
    line0 = line0.ljust(pos_R2) + p2 + rows_R2[0]
    line0 = line0.ljust(pos_R3) + p3 + rows_R3[0]
    print(line0)
    # Subsequent rows aligned under numbers
    line1 = i1 + rows_R1[1]
    line1 = line1.ljust(pos_R2) + i2 + rows_R2[1]
    line1 = line1.ljust(pos_R3) + i3 + rows_R3[1]
    print(line1)
    line2 = i1 + rows_R1[2]
    line2 = line2.ljust(pos_R2) + i2 + rows_R2[2]
    line2 = line2.ljust(pos_R3) + i3 + rows_R3[2]
    print(line2)
        
    print()
    
    print('N_ASG/N_SSG = ',spg_op_num//num_operator)
    print()
    if ssgnum != 'need more loop':
        print('More information about this SSG can be found at')
        print('https://cmpdc.iphy.ac.cn/ssg/ssgs/'+ssgnum)
    
    
    if abs(det33) > abs(det32)+0.5:
        is_super_cell = True
    else:
        is_super_cell = False
        
    return is_super_cell,cell_mag_unit
    
    
def get_ssg(file_name,ssg_list=None,tol=1e-3,tolm=1e-4):
    
    if ssg_list is None:
        ssg_list = load_ssg_list()
        
    if os.path.exists(file_name):

        cell,elements = read_poscar_no_elements(file_name)

    else:
        # print('Have no input files !!!')
        raise ValueError("Can't find input file!!!")
    
    lattice = cell[0].copy()
    position = cell[1].copy()
    numbers = cell[2].copy()
    mag = cell[-1].copy()
    
    if all(norm(m1) < tolm for m1 in mag):
        print('It is a non-magnetic structure')
        cell_nonmag = (lattice, position, numbers)
        nonmag_sym = get_symmetry_dataset(cell_nonmag)
        print('The SG number:',nonmag_sym['number'],'(',nonmag_sym['international'],')')
        sys.exit()
    else:
        operations = findAllOp(cell, tol = tol,tolm=tolm)
        
        # print(operations)
        ssgnum = search4ssg(cell, ssg_list, tol=tol, tolm=tolm)
        


        return ssg_list,ssgnum,operations

def build_argparser():
    parser = argparse.ArgumentParser(description='MOM2SSG: Magnetic Order to Spin Space Group Analysis')
    parser.add_argument('-c', type=str, default='POSCAR', help='Input file name')
    parser.add_argument('--standardize', action='store_true', help='Standardize the input')
    parser.add_argument('--tolerance', type=float, default=1e-3, help='Tolerance for space symmetry operations')
    parser.add_argument('--magtolerance', type=float, default=1e-4, help='Tolerance for spin operations')
    
    return parser


def main():
    """
    Main function for MOM2SSG analysis.
    
    This function parses command line arguments and performs SSG analysis
    on the input crystal structure file.
    """
    # Parse command line arguments
    parser = build_argparser()
    args = parser.parse_args()
    
    # Extract arguments to variables
    input_file = args.c  # Input file name
    standardize_flag = args.standardize  # Whether to standardize input
    tolerance = args.tolerance  # Tolerance for symmetry operations
    magtolerance = args.magtolerance  # Tolerance for spin operations
    
    warnings.filterwarnings('ignore')
    
    # Check for input file
    if os.path.exists(input_file):
        ssg_list = load_ssg_list()
        ssg_list_ = copy.deepcopy(ssg_list)
        
        if input_file.endswith('.mcif'):
            cell, elements = mcif2cell(input_file)
            print(f'Reading {input_file} as input......')
            vasp_input = False
        else:
            cell, elements = read_poscar_no_elements(input_file)
            print(f'Reading {input_file} as input......')
            vasp_input = True
    else:
        raise ValueError(f'Input file not found: {input_file}')
    
    lattice = cell[0].copy()
    position = cell[1].copy()
    numbers = cell[2].copy()
    mag = cell[3].copy()
    cell_nonmag = (lattice, position, numbers)
    nonmag_sym = get_symmetry_dataset(cell_nonmag,symprec=tolerance)
    
    if all(norm(m1) < magtolerance for m1 in mag):
        print('The SG number:',nonmag_sym['number'],'(',nonmag_sym['international'],')')
        sys.exit()
    else:
        
        ssg_ops = findAllOp(cell, tol=tolerance,tolm=magtolerance)
        np.save('ssgop.npy',ssg_ops,allow_pickle=True)
        dim_mag = findDimension(cell,tolm=magtolerance)
        
        if dim_mag == 1:
            axis_vector = [line_normal_vector(mag,tolm=magtolerance)]
            if abs(axis_vector[0][-1]) > 0.7:
                axis_vector[0] = axis_vector[0] / axis_vector[0][-1]* abs(axis_vector[0][-1])
            axis_vector.append(orthonormal_basis_from_vector(axis_vector[0])[0])
            spin_rot_list = [1 if det(i) > 0 else -1 for i in ssg_ops['spin']]
            
        elif dim_mag == 2:
            axis_vector = [generate_normal_vector(mag,tolm=magtolerance)]
            if abs(axis_vector[0][-1]) > 0.7:
                axis_vector[0] = axis_vector[0] / axis_vector[0][-1]* abs(axis_vector[0][-1])
                
            if ssg_ops['QLabel'] == '1' or ssg_ops['QLabel'] == 'm':
                axis_vector.append(orthonormal_basis_from_vector(axis_vector[0])[0])
            else:
                axis_vector1 = spin_axis(ssg_ops['spin'],ssg_ops['QLabel'])
                axis_vector.append(axis_vector1[1])
                
            
            spin_rot_list = [change_basis_O3(i,axis_vector[1],axis_vector[0],tol=magtolerance)[0:2,0:2] for i in ssg_ops['spin']]
            
            
            
        else:
            spin_rot_list = copy.deepcopy(ssg_ops['spin'])
            axis_vector = spin_axis(spin_rot_list,ssg_ops['QLabel'])
            
        axis_vector.append(np.cross(axis_vector[0],axis_vector[1]))
        
        cell_addelements = (lattice,position,numbers,elements,np.array(mag))
        cell_output_acc,wyckoff = get_swyckoff(cell_addelements, ssg_ops, tol=tolerance)
        output_wyckoff(cell_output_acc,wyckoff)
        write_poscar(cell_output_acc,file_name='POSCAR.symm')
        
        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        pg_path = os.path.join(pkg_dir, 'ssg_data', 'PG_dic_list.npy')
        if not os.path.exists(pg_path):
            raise FileNotFoundError(f"PG_dic_list.npy not found: {pg_path}")
        loaded_list = np.load(pg_path, allow_pickle=True)
        loaded_list = loaded_list.tolist()
        for pg in loaded_list:
            if pg['HM_label'] == ssg_ops['QLabel']:
                pg_op_num = len(pg['ops'])
                break
        else:
            pg_op_num = 0
            

        ssgnum = search4ssg(cell, ssg_list, tol = tolerance,tolm=magtolerance)

        if ssgnum[-1] == 'L':
            lps = 'Collinear magnetic configuration'
        elif ssgnum[-1] == 'P':
            lps = 'Coplanar magnetic configuration'
        elif ssgnum == 'need more loop':
            lps = ''
        else:
            lps = 'Noncoplanar magnetic configuration'
        
        if ssgnum != 'need more loop':
            format_ssg = get_SSG_label(ssgnum,ssg_list_)
        else:
            format_ssg = 'Cannot find SSG number and international symbol!!!'
        
        generate_irssg_in(ssg_ops['Gnum'], format_ssg, cell, mag, ssg_ops, tolm=magtolerance)
        
        is_super_cell, cell_mag_unit=format_output(dim_mag,axis_vector,spin_rot_list,ssg_ops,lps,pg_op_num,nonmag_sym,ssgnum,format_ssg,cell)
        
        if is_super_cell:
            print('Warning: The POSCAR cell is not a SSG primitive cell!!! The SSG primitive cell is output in POSCAR.ssg_primitive.')
            cell_addelements_magunit = (cell_mag_unit[0],cell_mag_unit[1],cell_mag_unit[2],elements,cell_mag_unit[-1])
            write_poscar(cell_addelements_magunit,file_name='POSCAR.ssg_primitive')
            
        
        
            
        if ssgnum != 'need more loop':
            if standardize_flag and vasp_input:
                standardize_ssg_cell(input_file,ssgnum,cell,ssg_ops,ssg_list_,tol=tolerance,tolm=magtolerance,symm_flag=True,check_flag=True)

if __name__ == '__main__':
    main()
