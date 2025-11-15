"""Convenience imports for IRSSG.

Expose commonly used functions directly under the ``irssg`` namespace
so external programs can do ``from irssg import ...``:

- ``findAllOp``: compute all symmetry operations for a magnetic cell
- ``load_ssg_list``: load the bundled Spin Space Group (SSG) dataset
- ``identify_SG_lattice(gid)``: return lattice type and primitive vectors (use [1] prim_vec)
- ``load_one_ssg(ssgnum)``: load one SSG entry
- ``pos2abr(cell)``: convert positions to A/B/R centered settings, return (cell_pos2abr, shift)
- ``generate_irssg_in(...)``: generate IRSSG input files from cell/args
- ``get_SSG_label(ssgnum, ssg_list)``: get formatted SSG label

This package provides the Python interface and data for IRSSG.
"""

from .ssg.find_ssg_operation import findAllOp  # re-export
from .ssg.load_ssgdata import load_ssg_list, load_one_ssg  # re-export
from .ssg.SG_utils import identify_SG_lattice  # re-export
from .ssg.pos2abr import pos2abr  # re-export
from .ssg.lirssg import generate_irssg_in  # re-export
from .ssg.get_ssg_symbol import get_SSG_label  # re-export

__all__ = [
    "findAllOp",
    "load_ssg_list",
    "identify_SG_lattice",
    "load_one_ssg",
    "pos2abr",
    "generate_irssg_in",
    "get_SSG_label",
]
