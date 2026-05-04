#!/usr/bin/env python3
"""
POS2SSG - Convert POSCAR to Spin Space Group (SSG) with magnetic moments

This script analyzes a crystal structure and generates a magnetic configuration
that corresponds to a specific spin space group (SSG).

Author: irssg Team
"""

import numpy as np
from spglib import *
import argparse
import sys
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from irssg import load_ssg_list, generate_irssg_in, get_SSG_label
import warnings
from .mag_atom import construct_magmom_dict, ELEMENT_TABLE
from .poscar_io import read_poscar, write_poscar
from .deal_cell import std_cell_pos2ssg, prim_cell, extend_cell
from .generate_mag import ssg2magmom
from .generate_structure import find_possible_ssg_numbers, generate_ssg_structure

warnings.filterwarnings('ignore')


def build_argparser():
    """Build command line argument parser."""
    parser = argparse.ArgumentParser(
        description='Convert POSCAR to Spin Space Group (SSG) with magnetic moments',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Find possible SSG numbers for a structure
  POS2SSG -c POSCAR

  # Find SSG numbers with specific Ik value
  POS2SSG -c POSCAR -Ik 2

  # Generate SPOSCAR for specific SSG
  POS2SSG -c POSCAR -ssg 43.2.4.1.P

  # Use custom magnetic atoms and moments
  POS2SSG -c POSCAR -ssg 43.2.4.1.P --magatom Fe Cu --magmom 3.0 2.0

  # Adjust tolerance for symmetry analysis
  POS2SSG -c POSCAR -ssg 43.2.4.1.P --tolerance 1e-3
        """
    )
    
    parser.add_argument('-c', '--poscar', type=str, default='POSCAR',
                       help='Input POSCAR file (default: POSCAR)')
    parser.add_argument('-Ik', type=int, required=False,
                       help='Filter by specific Ik value')
    parser.add_argument('-ssg', type=str, required=False,
                       help='Target spin space group number (e.g., 43.2.4.1.P)')
    parser.add_argument('--tolerance', type=float, default=1e-2,
                       help='Tolerance for symmetry analysis (default: 1e-2)')
    parser.add_argument('--magatom', type=str, nargs='*', default=[],
                       help='List of magnetic atom types (e.g., Fe Cu)')
    parser.add_argument('--magmom', type=float, nargs='*', default=[],
                       help='List of magnetic moments (e.g., 3.0 2.0)')
    parser.add_argument('--output', '-o', type=str, default='SPOSCAR',
                       help='Output file name (default: SPOSCAR)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    return parser


def validate_arguments(args):
    """Validate command line arguments."""
    if args.magatom and args.magmom and len(args.magatom) != len(args.magmom):
        print("Error: --magatom and --magmom lists must have the same length!")
        print(f"  magatom: {len(args.magatom)} items")
        print(f"  magmom: {len(args.magmom)} items")
        sys.exit(1)
    
    if args.tolerance <= 0:
        print("Error: tolerance must be positive!")
        sys.exit(1)





def analyze_magnetic_atoms(cell, magmom_dict, verbose=False):
    """Analyze magnetic atoms in the structure."""
    unique_numbers = sorted(set(cell[2]))  # cell[2] is numbers
    elements = [ELEMENT_TABLE[atomic_num] for atomic_num in unique_numbers]
    
    magnetic_elements = []
    for element in elements:
        if element in magmom_dict.keys():
            magnetic_elements.append(element)
    
    if verbose:
        print(f"Structure contains elements: {elements}")
        print(f"Magnetic elements: {magnetic_elements}")
        print(f"Magnetic moments: {magmom_dict}")
    
    return len(magnetic_elements) > 0, magnetic_elements


def main():
    """Main function."""
    # Parse arguments
    parser = build_argparser()
    args = parser.parse_args()
    
    # Validate arguments
    validate_arguments(args)
    
    # Load SSG database
    try:
        ssg_list = load_ssg_list()
    except Exception as e:
        print(f"Error loading SSG database: {e}")
        sys.exit(1)
    
    # Read structure
    try:
        structure = Structure.from_file(args.poscar)
        cell = read_poscar(file_name=args.poscar)
    except Exception as e:
        print(f"Error reading structure file {args.poscar}: {e}")
        sys.exit(1)
    
    # Construct magnetic moment dictionary
    magmom_dict = construct_magmom_dict(args.magatom, args.magmom)
    
    if args.verbose:
        print(f"Input file: {args.poscar}")
        print(f"Tolerance: {args.tolerance}")
        print(f"Magnetic atoms: {args.magatom}")
        print(f"Magnetic moments: {args.magmom}")
        print()
    
    if args.ssg is None:
        # Find possible SSG numbers
        find_possible_ssg_numbers(structure, ssg_list, args.Ik, args.tolerance)
    else:
        # Generate specific SSG structure
        # Check space group compatibility
        sga = SpacegroupAnalyzer(structure, symprec=args.tolerance)
        sgnum = sga.get_space_group_number()
        sg_in = int(args.ssg.split('.')[0])
        
        if sg_in != sgnum:
            print(f'Warning: The SSG number ({sg_in}) is not the same as the SG of the structure ({sgnum})!')
            print('Continuing with the analysis...')
        
        # Check for magnetic atoms
        has_mag_atoms, mag_elements = analyze_magnetic_atoms(cell, magmom_dict, args.verbose)
        
        if not has_mag_atoms:
            print("Error: No magnetic atoms found in the structure!")
            print("Please specify magnetic atoms using --magatom and --magmom options.")
            sys.exit(1)
        
        # Generate SSG structure
        success, output_cell, ssg_ops = generate_ssg_structure(
            args.ssg, cell, magmom_dict, ssg_list, 
            args.tolerance, args.verbose
        )
        
        
        if success:
            write_poscar(output_cell, args.output)
            print(f'Structure has been generated successfully in {args.output}.')
            format_ssg = get_SSG_label(args.ssg, ssg_list)
            # generate_irssg_in(ssg_ops['Gnum'], format_ssg, output_cell, output_cell[-1], ssg_ops, tolm=1e-4)
        else:
            print('Failed to generate SSG structure.')
            sys.exit(1)
            

if __name__ == "__main__":
    main()
