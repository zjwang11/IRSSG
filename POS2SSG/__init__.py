"""
POS2SSG - Convert POSCAR to Spin Space Group (SSG) with magnetic moments

This package provides tools for analyzing crystal structures and generating
magnetic configurations corresponding to specific spin space groups (SSG).
"""

__version__ = "1.0.1"
__author__ = "irssg Team"
__email__ = "songziyin@iphy.ac.cn, zhangsheng221@mails.ucas.ac.cn, wzj@iphy.ac.cn"

from .POS2SSG import main
from .mag_atom import construct_magmom_dict, ELEMENT_TABLE, analyze_magnetic_atoms
from .poscar_io import read_poscar, write_poscar
from .deal_cell import std_cell_pos2ssg, prim_cell, extend_cell
from .generate_mag import ssg2magmom, identify_wyckoff_positions

__all__ = [
    "main",
    "construct_magmom_dict",
    "ELEMENT_TABLE", 
    "analyze_magnetic_atoms",
    "read_poscar",
    "write_poscar",
    "std_cell_pos2ssg",
    "prim_cell",
    "extend_cell",
    "ssg2magmom",
    "identify_wyckoff_positions",
]
