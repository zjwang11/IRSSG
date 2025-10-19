"""Convenience imports for IRSSG.

Expose commonly used functions directly under the ``irssg`` namespace:

- ``findAllOp``: compute all symmetry operations for a magnetic cell
- ``load_ssg_list``: load the bundled Spin Space Group (SSG) dataset

This package provides the Python interface and data for IRSSG.
"""

from .ssg.find_ssg_operation import findAllOp  # re-export
from .ssg.load_ssgdata import load_ssg_list  # re-export

__all__ = ["findAllOp", "load_ssg_list"]
