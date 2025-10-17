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

from spglib import get_symmetry_dataset
import copy

# Local module imports
from .small_func import *
from .SG_utils import sg_symbol_from_number
from .SSGLabel import get_SSG_label
from .lirssg import generate_irssg_in
from .poscar_io import read_poscar_no_elements, mcif2cell
from .load_ssgdata import load_ssg_list
from .wyckoff import get_swyckoff, output_wyckoff
from .find_ssg_operation import findAllOp
from .get_ssg_number import search4ssg
from .std2ssg import standardize_ssg_cell
from .spintrans import spin_axis
from .eqvpg2label import get_std_pg


def format_output(dim_mag,axis_vector,spin_rot_list,operations,lps,pg_op_num,nonmag_sym,ssgnum,format_ssg):
    space_international = nonmag_sym['international']
    spg_op_num = nonmag_sym['rotations'].shape[0]
    print(f"Atomic space group: {nonmag_sym['number']}, {space_international}")
    ncell_pos_ssg = count_identity_pairs(operations['RotC'], operations['spin'], tol=1e-3)
    print(f"The volume of POSCAR is {ncell_pos_ssg} times of the SSG cell.")

    ncell_pos_asg = count_identity_pairs(operations['RotC'], [np.eye(3)]*len(operations['RotC']), tol=1e-3)
    assert ncell_pos_asg%ncell_pos_ssg == 0
    ncell_ssg_asg = int(ncell_pos_asg)//int(ncell_pos_ssg)
    print(f"The volume of SSG cell is {ncell_ssg_asg} times of ASG cell.")
    assert ncell_ssg_asg%operations['Ik'] == 0
    print(f"The volume of SSG_SPG cell is {ncell_ssg_asg//operations['Ik']} times of ASG cell.")
    print(f"The volume of SSG cell is {operations['Ik']} times of SSG_SPG cell.")
    
    num_operator = len(operations['spin'])

    assert spg_op_num%num_operator == 0
    print('N_ASG/N_SSG = ',spg_op_num//num_operator)
    
    print(lps)
       
    print(f"z' = [{axis_vector[0][0]:>6.3f}, {axis_vector[0][1]:>6.3f}, {axis_vector[0][2]:>6.3f}]")
    print(f"x' = [{axis_vector[1][0]:>6.3f}, {axis_vector[1][1]:>6.3f}, {axis_vector[1][2]:>6.3f}]")

    if ssgnum == 'need more loop':
        print('The SSG number can not be identified!')
    else:

        print(f'The SSG number: {ssgnum}')
        print(f'The SSG label: {format_ssg}')
        
    
    
    
    # print('Non-magnetic space group: ' + space_international + f' ({spg_op_num})',end = '  ;  ')
    print('CSG H(space part): '+ sg_symbol_from_number(operations['Gnum']),end = ',  ')
    print('NPG P(spin part): ' + get_std_pg(operations['QLabel'])[1])
    print('='*40)
    
    print('Spin space group operations:')

    print(f'# ssg_symm: {num_operator}')
    if dim_mag == 1:
        print('{Spin|| Ri  | taui}')
    elif dim_mag == 2:
        print('{   Spin O2 || Ri  | taui}')
    else:
        print('{      Spin O3    || Ri  | taui}')
    for i in range(num_operator):
        time_reversal = det(operations['spin'][i])
        if time_reversal > 0:
            time_reversal_str = 'without time reversal'
        else:
            time_reversal_str = 'with time reversal'
            
        print(f'# {i+1:>3d}   {time_reversal_str}')
        
        operations["TauC"][i] = (operations["TauC"][i]+0.5)%1-0.5
        if dim_mag == 1:
            print(f"{spin_rot_list[i]:>4d}",end='  ')
            print("".join(f"{v:>2d}" for v in np.asarray(operations["RotC"][i][0, :].reshape(-1), int)),end=' ')
            print("".join(f"{v:>6.3f}" for v in np.asarray(operations["TauC"][i][0].reshape(-1), float)))
            
            print("      ",end='')
            print("".join(f"{v:>2d}" for v in np.asarray(operations["RotC"][i][1, :].reshape(-1), int)),end=' ')
            print("".join(f"{v:>6.3f}" for v in np.asarray(operations["TauC"][i][1].reshape(-1), float)))
            
            print("      ",end='')
            print("".join(f"{v:>2d}" for v in np.asarray(operations["RotC"][i][2 :].reshape(-1), int)),end=' ')
            print("".join(f"{v:>6.3f}" for v in np.asarray(operations["TauC"][i][2].reshape(-1), float)))


        elif dim_mag == 2:
            print("".join(f"{v:>6.3f}" for v in np.asarray(spin_rot_list[i][0, :], float)),end=' ')
            print("".join(f"{v:>2d}" for v in np.asarray(operations["RotC"][i][0, :].reshape(-1), int)),end=' ')
            print("".join(f"{v:>6.3f}" for v in np.asarray(operations["TauC"][i][0].reshape(-1), float)))
            
            
            print("".join(f"{v:>6.3f}" for v in np.asarray(spin_rot_list[i][1, :], float)),end=' ')
            print("".join(f"{v:>2d}" for v in np.asarray(operations["RotC"][i][1, :].reshape(-1), int)),end=' ')
            print("".join(f"{v:>6.3f}" for v in np.asarray(operations["TauC"][i][1].reshape(-1), float)))
            
            print(" "*13,end='')
            print("".join(f"{v:>2d}" for v in np.asarray(operations["RotC"][i][2 :].reshape(-1), int)),end=' ')
            print("".join(f"{v:>6.3f}" for v in np.asarray(operations["TauC"][i][2].reshape(-1), float)))


        else:
            print("".join(f"{v:>6.3f}" for v in np.asarray(spin_rot_list[i][0, :], float)),end=' ')
            print("".join(f"{v:>2d}" for v in np.asarray(operations["RotC"][i][0, :].reshape(-1), int)),end=' ')
            print("".join(f"{v:>6.3f}" for v in np.asarray(operations["TauC"][i][0].reshape(-1), float)))
            
            
            print("".join(f"{v:>6.3f}" for v in np.asarray(spin_rot_list[i][1, :], float)),end=' ')
            print("".join(f"{v:>2d}" for v in np.asarray(operations["RotC"][i][1, :].reshape(-1), int)),end=' ')
            print("".join(f"{v:>6.3f}" for v in np.asarray(operations["TauC"][i][1].reshape(-1), float)))
            
            print("".join(f"{v:>6.3f}" for v in np.asarray(spin_rot_list[i][2, :], float)),end=' ')
            print("".join(f"{v:>2d}" for v in np.asarray(operations["RotC"][i][2 :].reshape(-1), int)),end=' ')
            print("".join(f"{v:>6.3f}" for v in np.asarray(operations["TauC"][i][2].reshape(-1), float)))

def get_ssg(file_name,ssg_list=None,tol=1e-3,tolm=1e-4):
    
    if ssg_list is None:
        ssg_list = load_ssg_list()
        
    if os.path.exists(file_name):

        cell,elements = read_poscar_no_elements(file_name)

    else:
        # print('Have no input files !!!')
        raise ValueError('cant find input file!!!')
    
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
        
        operations = findAllOp(cell, tol=tolerance,tolm=magtolerance)
        np.save('ssgop.npy',operations,allow_pickle=True)
        dim_mag = findDimension(cell,tolm=magtolerance)
        
        if dim_mag == 1:
            axis_vector = [line_normal_vector(mag,tolm=magtolerance)]
            if abs(axis_vector[0][-1]) > 0.7:
                axis_vector[0] = axis_vector[0] / axis_vector[0][-1]* abs(axis_vector[0][-1])
            axis_vector.append(orthonormal_basis_from_vector(axis_vector[0])[0])
            spin_rot_list = [1 if det(i) > 0 else -1 for i in operations['spin']]
            
        elif dim_mag == 2:
            axis_vector = [generate_normal_vector(mag,tolm=magtolerance)]
            if abs(axis_vector[0][-1]) > 0.7:
                axis_vector[0] = axis_vector[0] / axis_vector[0][-1]* abs(axis_vector[0][-1])
                
            if operations['QLabel'] == '1' or operations['QLabel'] == 'm':
                axis_vector.append(orthonormal_basis_from_vector(axis_vector[0])[0])
            else:
                axis_vector1 = spin_axis(operations['spin'],operations['QLabel'])
                axis_vector.append(axis_vector1[1])
                
            
            spin_rot_list = [change_basis_O3(i,axis_vector[1],axis_vector[0],tol=magtolerance)[0:2,0:2] for i in operations['spin']]
            
            
            
        else:
            spin_rot_list = copy.deepcopy(operations['spin'])
            axis_vector = spin_axis(spin_rot_list,operations['QLabel'])
            
        
        
        cell_addelements = (lattice,position,numbers,elements,np.array(mag))
        cell_sc_output_acc,wyckoff = get_swyckoff(cell_addelements, operations, tol=tolerance)
        output_wyckoff(cell_sc_output_acc,wyckoff)
        
        
        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        pg_path = os.path.join(pkg_dir, 'ssg_data', 'PG_dic_list.npy')
        if not os.path.exists(pg_path):
            raise FileNotFoundError(f"PG_dic_list.npy not found: {pg_path}")
        loaded_list = np.load(pg_path, allow_pickle=True)
        loaded_list = loaded_list.tolist()
        for pg in loaded_list:
            if pg['HM_label'] == operations['QLabel']:
                pg_op_num = len(pg['ops'])
            

        ssgnum = search4ssg(cell, ssg_list, tol = tolerance,tolm=magtolerance)

        if ssgnum[-1] == 'L':
            lps = 'Collinear magnetic configuration'
        elif ssgnum[-1] == 'P':
            lps = 'Coplanar magnetic configuration'
        elif ssgnum == 'need more loop':
            lps = ''
        else:
            lps = 'Spatial magnetic configuration'
        
        if ssgnum != 'need more loop':
            format_ssg = get_SSG_label(ssgnum,ssg_list_)
        else:
            format_ssg = 'Cannot find SSG number and international symbol!!!'
        
        generate_irssg_in(operations['Gnum'], format_ssg, cell, mag, operations, tolm=magtolerance)
        
        format_output(dim_mag,axis_vector,spin_rot_list,operations,lps,pg_op_num,nonmag_sym,ssgnum,format_ssg)
        
        if ssgnum != 'need more loop':
            print('More information about this SSG can be found at')
            print('https://cmpdc.iphy.ac.cn/ssg/ssgs/'+ssgnum)
        
            if standardize_flag and vasp_input:
                standardize_ssg_cell(input_file,ssgnum,cell,operations,ssg_list_,tol=tolerance,tolm=magtolerance,symm_flag=True,check_flag=True)

if __name__ == '__main__':
    main()
