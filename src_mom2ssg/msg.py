"""
Utilities to detect magnetic space group (MSG) identifiers using spglib.
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np
import spglib


def _get_field(obj, name: str):
    """
    spglib is migrating from dict-like dataset to attribute-like dataset.
    This getter supports both styles.
    """
    if hasattr(obj, name):
        return getattr(obj, name)
    try:
        return obj[name]
    except Exception:
        return None


def get_msg_info_from_cell(
    cell,
    symprec: float = 1e-3,
    angle_tolerance: float = 5.0,
    mag_symprec: float = 1e-2,
) -> Optional[Dict[str, object]]:
    """
    Return magnetic space group identifiers from a cell tuple.

    Args:
        cell: (lattice, positions, numbers, magmoms)
        symprec: spglib tolerance for atomic positions
        angle_tolerance: spglib tolerance for angles (degrees)
        mag_symprec: spglib tolerance for magnetic moments
    """
    if cell is None or len(cell) < 4:
        raise ValueError("cell must be (lattice, positions, numbers, magmoms)")

    lattice, positions, numbers, magmoms = cell[:4]
    magmoms = np.asarray(magmoms, dtype=float)

    spg_cell = (
        np.asarray(lattice, dtype=float),
        np.asarray(positions, dtype=float),
        np.asarray(numbers, dtype=int),
        magmoms,
    )

    try:
        dataset = spglib.get_magnetic_symmetry_dataset(
            spg_cell,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            mag_symprec=mag_symprec,
        )
    except Exception:
        return None

    if dataset is None:
        return None

    uni_number = _get_field(dataset, "uni_number")
    if uni_number is None:
        return None

    try:
        uni_number = int(uni_number)
    except Exception:
        pass

    try:
        msg = spglib.get_magnetic_spacegroup_type(int(uni_number))
    except Exception:
        msg = None

    info = {"uni_number": uni_number}
    if msg is not None:
        for key in (
            "bns_number",
            "og_number",
            "litvin_number",
            "type",
            "bns_symbol",
            "og_symbol",
        ):
            info[key] = _get_field(msg, key)

    return info


def format_msg_label(info: Optional[Dict[str, object]]) -> str:
    """
    Format MSG identifiers into a compact string for logging/output.
    """
    if not info:
        return "MSG: unknown"

    def clean(val):
        if val is None:
            return None
        text = str(val).strip()
        return text if text else None

    og_symbol = clean(info.get("og_symbol"))
    if og_symbol:
        return og_symbol

    og_number = clean(info.get("og_number"))
    if og_number:
        return og_number

    return "MSG: unknown"


def _main():
    import argparse
    try:
        from .poscar_io import read_poscar_no_elements
    except Exception:
        from poscar_io import read_poscar_no_elements

    parser = argparse.ArgumentParser(
        description="Detect magnetic space group identifiers from a POSCAR."
    )
    parser.add_argument("poscar", help="POSCAR-like file (may include magmom).")
    parser.add_argument("--symprec", type=float, default=1e-3)
    parser.add_argument("--angle_tolerance", type=float, default=5.0)
    parser.add_argument("--mag_symprec", type=float, default=1e-2)
    args = parser.parse_args()

    cell, _ = read_poscar_no_elements(args.poscar)
    info = get_msg_info_from_cell(
        cell,
        symprec=args.symprec,
        angle_tolerance=args.angle_tolerance,
        mag_symprec=args.mag_symprec,
    )
    print(format_msg_label(info))
    if info:
        for key in (
            "uni_number",
            "bns_number",
            "bns_symbol",
            "og_number",
            "og_symbol",
            "litvin_number",
            "type",
        ):
            if key in info:
                print(f"{key:>15s} : {info[key]}")


if __name__ == "__main__":
    _main()
