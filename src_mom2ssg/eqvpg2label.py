import re

def hm_to_schoenflies(hm: str) -> str:
    """
    Convert Hermann–Mauguin (international) point-group symbol to Schoenflies.
    Prefer using pymatgen for parsing; fall back to a lightweight mapping if unavailable.
    """
    if hm is None:
        raise ValueError("hm must be a non-empty string")
    if not isinstance(hm, str):
        # For non-string labels, behave as unsupported (KeyError).
        raise KeyError(f"Unknown or unsupported HM point-group label: '{hm}'")

    # ---------- 1) Preprocessing and synonyms ----------
    raw = hm.strip()
    if raw == "":
        raise KeyError("Unknown or unsupported HM point-group label: ''")
    s = raw.replace(" ", "").replace(".", "").lower()

    # Support "bar" notation (e.g., 4bar3m -> -43m) and simple reordering
    if "bar" in s:
        s = s.replace("bar", "-")
        # Reorder tokens like "3-"/"4-"/"6-"/"1-" into "-3"/"-4"/"-6"/"-1"
        for d in ("1", "3", "4", "6"):
            s = s.replace(d + "-", "-" + d)

    # Normalize common synonyms (formatting only; same meaning) to help library matching
    if s == "m3m":
        s = "m-3m"
    # Orthorhombic mm2 family synonyms
    if s in {"2mm", "m2m"}:
        s = "mm2"
    # Normalize to commonly preferred forms
    if s in {"-6m2", "-62m"}:
        s = "-6m2"
    if s in {"-4m2", "-42m"}:
        s = "-42m"
    # Common trigonal variants (same 3m point group)
    if s in {"31m", "3m1"}:
        s = "3m"

    # ---------- 2) Resolve via pymatgen first ----------
    try:
        from pymatgen.symmetry.groups import PointGroup

        # 2a) Direct construction; handle attribute-name differences across versions
        for candidate in (s, raw):
            try:
                pg = PointGroup(candidate)
            except Exception:
                continue
            for attr in ("sch_symbol", "schoenflies_symbol", "schoenflies"):
                if hasattr(pg, attr):
                    val = getattr(pg, attr)
                    if val:
                        return str(val)

        # 2b) Iterate all point groups and match by international symbol
        try:
            from pymatgen.symmetry.groups import all_pointgroups

            def norm(x: str) -> str:
                return x.replace(" ", "").replace(".", "").lower()

            target = norm(s)
            for pg in all_pointgroups:
                intl = getattr(pg, "intl_symbol", None)
                if intl and norm(str(intl)) == target:
                    for attr in ("sch_symbol", "schoenflies_symbol", "schoenflies"):
                        if hasattr(pg, attr):
                            val = getattr(pg, attr)
                            if val:
                                return str(val)
        except Exception:
            pass
    except Exception:
        pass

    # ---------- 3) General HM pattern parsing (molecular + crystallographic) ----------
    # Handle families like Cn/Cnv/Cnh, Dn/Dnh/Dnd, Sn, and icosahedral variants.
    import re as _re

    # 3a) If not one of the explicit cubic/icosahedral symbols, try family patterns
    if s not in {"23", "m-3", "-43m", "-4m3", "432", "m-3m", "532", "m-3-5", "-5-3m"}:
        # Cnv: nmm (even n) or nm (odd n) -> Cnv
        m = _re.fullmatch(r"(\d+)mm", s)
        if m:
            n = int(m.group(1))
            return f"C{n}v"
        m = _re.fullmatch(r"(\d+)m", s)
        if m:
            n = int(m.group(1))
            return f"C{n}v"

        # Cnh: n/m (even n) or -(2n) (odd n) -> Cnh
        m = _re.fullmatch(r"(\d+)/m", s)
        if m:
            n = int(m.group(1))
            return f"C{n}h"

        # Dn: n22 (even n) or n2 (odd n) -> Dn
        m = _re.fullmatch(r"(\d+)22", s)
        if m:
            n = int(m.group(1))
            return f"D{n}"
        m = _re.fullmatch(r"(\d+)2", s)
        if m:
            n = int(m.group(1))
            return f"D{n}"

        # Dnh: n/mmm (even n) or -(2n)m2 (odd n)
        m = _re.fullmatch(r"(\d+)/mmm", s)
        if m:
            n = int(m.group(1))
            return f"D{n}h"
        m = _re.fullmatch(r"-(\d+)m2", s)
        if m:
            k = int(m.group(1))
            if k % 2 == 0 and k >= 2:
                return f"D{k//2}h"

        # Dnd: -n2m (even n) or -nm (odd n)
        m = _re.fullmatch(r"-(\d+)2m", s)
        if m:
            n = int(m.group(1))
            return f"D{n}d"
        m = _re.fullmatch(r"-(\d+)m", s)
        if m:
            n = int(m.group(1))
            return f"D{n}d"

        # Pure cyclic: n -> Cn
        m = _re.fullmatch(r"(\d+)", s)
        if m:
            return f"C{int(m.group(1))}"

        # Improper rotation or Cnh (odd n): -k
        m = _re.fullmatch(r"-(\d+)", s)
        if m:
            k = int(m.group(1))
            if k == 1:
                return "Ci"
            if k == 3:
                return "C3i"  # prefer C3i over S3
            if k % 2 == 0:
                # If k/2 is odd, prefer C(k/2)h (e.g., -6 -> C3h, -10 -> C5h); else Sn
                return f"C{(k // 2)}h" if (k // 2) % 2 == 1 else f"S{k}"
            else:
                return f"S{k}"

    # ---------- 4) Lightweight fallback mapping (common cases, plus icosahedral) ----------
    # Explicit map for the 32 crystallographic groups and key icosahedral symbols
    hm2sch = {
        # triclinic
        "1": "C1",
        "-1": "Ci",

        # monoclinic
        "2": "C2",
        "m": "Cs",
        "2/m": "C2h",

        # orthorhombic
        "222": "D2",
        "mm2": "C2v",
        "m2m": "C2v",
        "mmm": "D2h",

        # tetragonal
        "4": "C4",
        "-4": "S4",
        "4/m": "C4h",
        "422": "D4",
        "4mm": "C4v",
        "-42m": "D2d",
        "-4m2": "D2d",
        "4/mmm": "D4h",

        # trigonal
        "3": "C3",
        "-3": "C3i",     # = S6
        "32": "D3",
        "3m": "C3v",
        "31m": "C3v",    # synonym
        "3m1": "C3v",    # synonym
        "-3m": "D3d",

        # hexagonal
        "6": "C6",
        "-6": "C3h",
        "6/m": "C6h",
        "622": "D6",
        "6mm": "C6v",
        "-62m": "D3h",
        "-6m2": "D3h",
        "6/mmm": "D6h",

        # cubic
        "23": "T",
        "m-3": "Th",
        "-43m": "Td",
        "-4m3": "Td",
        "432": "O",
        "m-3m": "Oh",

        # icosahedral (non-crystallographic, molecular convention)
        "532": "I",
        "m-3-5": "Ih",
        "-5-3m": "Ih",
    }

    key = s
    if key in hm2sch:
        return hm2sch[key]

    # If neither path works, raise a unified error
    raise KeyError(
        f"Unsupported or ambiguous HM point-group symbol: {raw!r}. "
        "Try installing pymatgen (pip install pymatgen) for broader parsing."
    )


def schoenflies_to_hm(sch: str) -> str:
    """
    Convert Schoenflies notation to Hermann–Mauguin (international) symbol.
    - Prefer crystallographic canonical HM where applicable.
    - Handle common non-crystallographic groups (I, Ih) approximately.

    Exceptions align with hm_to_schoenflies:
      - None -> ValueError("...non-empty string")
      - Non-string/unsupported -> KeyError
    """
    if sch is None:
        raise ValueError("sch must be a non-empty string")
    if not isinstance(sch, str):
        raise KeyError(f"Unknown or unsupported Schoenflies label: '{sch}'")

    raw = sch.strip()
    if raw == "":
        raise KeyError("Unknown or unsupported Schoenflies label: ''")

    # Canonicalize basic casing like C2v, D3h, Td, Oh, etc.
    s = raw.replace(" ", "").replace(".", "")

    # Exact named groups first (icosahedral, tetrahedral, octahedral, reflections, inversion)
    exact = {
        "I": "532",
        "Ih": "-5-3m",
        "T": "23",
        "Td": "-43m",
        "Th": "m-3",
        "O": "432",
        "Oh": "m-3m",
        "Cs": "m",
        "Ci": "-1",
        # Useful synonym: C3i is equivalent to S6 (HM -3)
        "C3i": "-3",
    }
    if s in exact:
        return exact[s]

    # Try to use pymatgen if available by inverting the point-group list
    try:
        from pymatgen.symmetry.groups import all_pointgroups

        def norm(x: str) -> str:
            return x.replace(" ", "").replace(".", "")

        target = norm(s)
        for pg in all_pointgroups:
            for attr in ("sch_symbol", "schoenflies_symbol", "schoenflies"):
                if hasattr(pg, attr):
                    sch_sym = getattr(pg, attr)
                    if sch_sym and norm(str(sch_sym)) == target:
                        intl = getattr(pg, "intl_symbol", None)
                        if intl:
                            return str(intl)
        # If not matched, continue to manual rules
    except Exception:
        pass

    # Manual family rules (generic to molecular point groups as well).
    # Support: Cn, Cnv, Cnh, Sn (any n), Dn, Dnh, Dnd.
    import re as _re
    m = _re.fullmatch(r"([CDS])(\d+)([hvd]?|i)", s)
    if m:
        prefix, n_str, suf = m.groups()
        n = int(n_str)
        suf = (suf or "").lower()
        P = prefix.upper()

        if P == "C":
            if suf == "":
                return str(n)
            if suf == "v":
                if n == 2:
                    return "mm2"  # canonical for C2v
                if n == 3:
                    return "3m"   # canonical for C3v
                return f"{n}mm" if n % 2 == 0 else f"{n}m"
            if suf == "h":
                return f"{n}/m" if n % 2 == 0 else f"-{2*n}"
            if suf == "i":
                # Cni only relevant crystallographically for n=3 (C3i == S6)
                if n == 3:
                    return "-3"
                raise KeyError(f"Unknown or unsupported Schoenflies label: '{sch}'")
            raise KeyError(f"Unknown or unsupported Schoenflies label: '{sch}'")

        if P == "S":
            # Generalized improper-rotation mapping for any n
            if n == 2:
                return "-1"  # S2 == Ci
            if n == 4:
                return "-4"
            if n == 6:
                return "-3"  # S6 == C3i
            return f"-{2*n}" if n % 2 == 0 else f"-{n}"

        if P == "D":
            if suf == "":
                # Dn
                return f"{n}22" if n % 2 == 0 else f"{n}2"
            if suf == "h":
                # Dnh
                return f"{n}/mmm" if n % 2 == 0 else f"-{2*n}m2"
            if suf == "d":
                # Dnd
                if n == 2:
                    return "-42m"  # special canonical for D2d
                return f"-{n}2m" if n % 2 == 0 else f"-{n}m"
            raise KeyError(f"Unknown or unsupported Schoenflies label: '{sch}'")

    # No rule matched
    raise KeyError(f"Unknown or unsupported Schoenflies label: '{sch}'")




def get_std_pg(pg_label):
    cry_TF = True
    check_set = {"C", "S", "D", "T", "O", "I"}
    if any(ch in pg_label for ch in check_set):
        cry_TF = False
    # print(pg_label, cry_TF)
    if cry_TF:
        sch_label = hm_to_schoenflies(pg_label)
    else:
        sch_label = pg_label
    # print(sch_label)
    # get the H-M label
    if cry_TF:
        HM_label =pg_label
    else:
        HM_label = schoenflies_to_hm(pg_label)
    return HM_label, sch_label



if __name__ == '__main__':
    print(get_std_pg('32'))