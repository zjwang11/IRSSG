"""Compatibility shim for the legacy package layout.

When the build installs modules under ``irssg_pkg`` (legacy), this module
provides the modern ``irssg`` import path by forwarding symbols and submodules.
"""

from importlib import import_module
import sys

# Re-export commonly used functions from legacy package
_find_ssg_op = import_module('irssg_pkg.ssg.find_ssg_operation')
_load_ssgdata = import_module('irssg_pkg.ssg.load_ssgdata')

findAllOp = getattr(_find_ssg_op, 'findAllOp')
load_ssg_list = getattr(_load_ssgdata, 'load_ssg_list')

# Expose submodule alias so that `import irssg.ssg.*` works
_ssg = import_module('irssg_pkg.ssg')
sys.modules[__name__ + '.ssg'] = _ssg

__all__ = ["findAllOp", "load_ssg_list"]

