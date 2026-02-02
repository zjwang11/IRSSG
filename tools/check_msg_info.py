#!/usr/bin/env python3
"""Check MSG identifiers (BNS/OG) using spglib for a POSCAR file.

Usage:
  python tools/check_msg_info.py POSCAR
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _add_repo_root_to_path() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root))


def main() -> int:
    _add_repo_root_to_path()

    from src_mom2ssg.msg import format_msg_label, get_msg_info_from_cell
    from src_mom2ssg.poscar_io import read_poscar_no_elements
    import spglib

    parser = argparse.ArgumentParser(description="Check MSG identifiers from a POSCAR.")
    parser.add_argument("poscar", nargs="?", default="POSCAR", help="POSCAR file path")
    parser.add_argument("--symprec", type=float, default=1e-3)
    parser.add_argument("--angle-tolerance", type=float, default=5.0)
    parser.add_argument("--mag-symprec", type=float, default=1e-2)
    args = parser.parse_args()

    cell, _ = read_poscar_no_elements(args.poscar)
    info = get_msg_info_from_cell(
        cell,
        symprec=args.symprec,
        angle_tolerance=args.angle_tolerance,
        mag_symprec=args.mag_symprec,
    )

    print(f"spglib_version: {spglib.__version__}")
    print(f"msg_info: {info}")

    if info:
        bns_symbol, bns_number, og_number = format_msg_label(info, return_pair=True)
        print(f"bns_symbol: {bns_symbol}")
        print(f"bns_number: {bns_number}")
        print(f"og_number: {og_number}")
        print(f"label: {format_msg_label(info)}")
    else:
        print("msg_info is None (no MSG detected)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
