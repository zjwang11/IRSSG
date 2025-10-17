import re


def hm_to_schoenflies(hm):
    if hm is None:
        raise ValueError("hm must be a non-empty string")


    mapping = {
        # Triclinic
        "1": "C1",
        "-1": "Ci",

        # Monoclinic
        "2": "C2",
        "m": "Cs",
        "2/m": "C2h",

        # Orthorhombic
        "222": "D2",
        "mm2": "C2v",
        "mmm": "D2h",

        # Tetragonal
        "4": "C4",
        "-4": "S4",
        "4/m": "C4h",
        "422": "D4",
        "4mm": "C4v",
        "-42m": "D2d",
        "4/mmm": "D4h",

        # Trigonal
        "3": "C3",
        "-3": "C3i",   # a.k.a. S6
        "32": "D3",
        "3m": "C3v",
        "-3m": "D3d",

        # Hexagonal
        "6": "C6",
        "-6": "C3h",   # crystallographic '-6' corresponds to Schoenflies C3h
        "6/m": "C6h",
        "622": "D6",
        "6mm": "C6v",
        "-6m2": "D3h",
        "6/mmm": "D6h",

        # Cubic
        "23": "T",
        "m-3": "Th",
        "432": "O",
        "-43m": "Td",
        "m-3m": "Oh",
    }

    if hm in mapping:
        return mapping[hm]
    else:
        raise KeyError(f"Unknown or unsupported HM point-group label: '{hm}'")




def schoenflies_to_hm(sch: str) -> str:
    """
    Convert Schoenflies notation to Hermann-Mauguin (approximate rules for non-crystallographic groups).
    Rules:
      - Cn   -> n
      - Cnv  -> nmm (if n even), nm (if n odd)
      - Cnh  -> n/m (if n even), -2n (if n odd)
      - Sn   -> -2n (if n even), -n (if n odd)
      - Dn   -> n22 (if n even), n2 (if n odd)
      - Dnh  -> n/mmm (if n even), -(2n)m2 (if n odd)
      - Dnd  -> -(2n)2m (if n even), -nm (if n odd)
      - I    -> 532
      - Ih   -> -5-3m
      - T    -> 23
      - Td   -> -43m
      - Th   -> m-3
      - O    -> 432
      - Oh   -> m-3m
      - Cs   -> m
      - Ci   -> -1
    """

    s = sch.strip()

    # Special cases
    if s == "I":
        return "532"
    if s == "Ih":
        return "-5-3m"
    
    # T groups (tetrahedral)
    if s == "T":
        return "23"
    if s == "Td":
        return "-43m"
    if s == "Th":
        return "m-3"
    
    # O groups (octahedral)
    if s == "O":
        return "432"
    if s == "Oh":
        return "m-3m"
    
    # C groups (cyclic)
    if s == "Cs":
        return "m"
    if s == "Ci":
        return "-1"

    # General matching
    m = re.match(r"([CDS])(\d+)([hvd]*)", s, re.IGNORECASE)
    if not m:
        raise ValueError(f"Unsupported Schoenflies label: {sch}")

    prefix, n_str, suffix = m.groups()
    n = int(n_str)
    suffix = suffix.lower()

    # Handle cases
    if prefix.upper() == "C":
        if suffix == "":
            return str(n)
        elif suffix == "v":
            return f"{n}mm" if n % 2 == 0 else f"{n}m"
        elif suffix == "h":
            return f"{n}/m" if n % 2 == 0 else f"-{2*n}"
        else:
            raise ValueError(f"Unsupported Cn suffix: {suffix}")

    elif prefix.upper() == "S":
        return f"-{2*n}" if n % 2 == 0 else f"-{n}"

    elif prefix.upper() == "D":
        if suffix == "":
            return f"{n}22" if n % 2 == 0 else f"{n}2"
        elif suffix == "h":
            return f"{n}/mmm" if n % 2 == 0 else f"-{2*n}m2"
        elif suffix == "d":
            return f"-{2*n}2m" if n % 2 == 0 else f"-{n}m"
        else:
            raise ValueError(f"Unsupported Dn suffix: {suffix}")

    else:
        raise ValueError(f"Unsupported Schoenflies prefix: {prefix}")



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