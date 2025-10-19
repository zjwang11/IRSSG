"""SSG analysis module vendored under irssg_pkg.ssg.

Convenience re-exports for common entry points:
    - findAllOp: compute all symmetry operations for a magnetic cell
    - load_ssg_list: load the bundled Spin Space Group (SSG) dataset
"""

from .find_ssg_operation import findAllOp, findAllOp_v2  # noqa: F401
from .load_ssgdata import (
    load_ssg_list,
    load_one_ssg,
    construct_std_operations,
)  # noqa: F401

__all__ = [
    "findAllOp",
    "findAllOp_v2",
    "load_ssg_list",
    "load_one_ssg",
    "construct_std_operations",
]
