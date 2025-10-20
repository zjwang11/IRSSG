import re
import sys

# Separator: a line consisting of only asterisks (and optional whitespace)
SEP_RE = re.compile(r"^\s*\*+\s*$")


def _is_sep(line: str) -> bool:
    return bool(SEP_RE.match(line))


def _extract_name_from_block(block_lines):
    # Find first non-empty line in the block
    first = ""
    for ln in block_lines:
        if ln.strip():
            first = ln
            break
    # Extract name from a line like: "kname = ... Fractional" (case-insensitive)
    m = re.search(r"kname\s*=\s*(.*?)\s*Fractional", first, re.IGNORECASE)
    return (m.group(1).strip() if m else first.strip())


def dedup_fort154(in_path="fort.154", out_path="chart.dat", comprel_path="comprel.log"):
    with open(in_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # Parse into blocks: content between two separator lines is a block
    blocks = []  # list of tuples: (start_sep, block_lines, end_sep)
    start_sep = None
    cur_lines = None

    for line in lines:
        if _is_sep(line):
            if start_sep is None:
                start_sep = line
                cur_lines = []
            else:
                # End current block
                blocks.append((start_sep, cur_lines or [], line))
                # This separator becomes the start of the next block
                start_sep = line
                cur_lines = []
        else:
            if cur_lines is not None:
                cur_lines.append(line)

    # If file doesn't end with a separator but we have content, finalize last block
    if start_sep is not None and cur_lines is not None and cur_lines:
        blocks.append((start_sep, cur_lines, ""))

    # Deduplicate consecutive blocks with the same extracted name and keep names
    prev_name = None
    kept = []            # tuples: (start, blines, end)
    kept_names = []      # corresponding names
    for start, blines, end in blocks:
        name = _extract_name_from_block(blines)
        if prev_name is not None and name == prev_name:
            pass
        else:
            kept.append((start, blines, end))
            kept_names.append(name)
        prev_name = name

    # Parse comprel.log and map unordered kpoint pairs to relation blocks
    pair_to_blocks = {}
    try:
        with open(comprel_path, "r", encoding="utf-8", errors="ignore") as f:
            clines = f.readlines()
    except FileNotFoundError:
        clines = []

    if clines:
        eq_re = re.compile(r"^\s*=+\s*$")

        def norm_kname(s: str) -> str:
            s = s.strip()
            s = re.sub(r"\d+", "", s)
            s = re.sub(r"\s+", "", s)
            s = re.sub(r"[^A-Za-z]", "", s)
            return s.upper()

        cblocks = []  # list of list[str]
        cur = None
        for ln in clines:
            if eq_re.match(ln):
                if cur is None:
                    cur = []
                else:
                    cblocks.append(cur)
                    cur = []
            else:
                if cur is not None:
                    cur.append(ln)
        if cur:
            cblocks.append(cur)

        for b in cblocks:
            first_nonempty = ""
            for ln in b:
                if ln.strip():
                    first_nonempty = ln.strip()
                    break
            if not first_nonempty or "->" not in first_nonempty:
                continue
            left, right = first_nonempty.split("->", 1)
            k3 = norm_kname(left)
            k4 = norm_kname(right)
            if not k3 or not k4:
                continue
            key = frozenset((k3, k4))
            pair_to_blocks.setdefault(key, []).append(b)

    def norm_name_for_chart(name: str) -> str:
        s = re.sub(r"\d+", "", name)
        s = re.sub(r"\s+", "", s)
        s = re.sub(r"[^A-Za-z]", "", s)
        return s.upper()

    for i in range(len(kept)):
        if i == 0:
            continue
        k1 = norm_name_for_chart(kept_names[i])
        k2 = norm_name_for_chart(kept_names[i - 1])
        key = frozenset((k1, k2))
        if key in pair_to_blocks:
            start, blines, end = kept[i]
            if blines and blines[-1] and not blines[-1].endswith("\n"):
                blines[-1] = blines[-1] + "\n"
            blines.append("Compatibility relations:\n")
            for b in pair_to_blocks[key]:
                blines.extend(b)
                if not blines or blines[-1].strip():
                    blines.append("\n")
            kept[i] = (start, blines, end)

    # Write output
    with open(out_path, "w", encoding="utf-8") as out:
        for i, (start, blines, end) in enumerate(kept):
            if i == 0:
                out.write(start)
            out.writelines(blines)
            if end:
                out.write(end)


def main(argv=None) -> int:
    args = list(argv) if argv is not None else sys.argv[1:]
    in_path = args[0] if len(args) > 0 else "fort.154"
    out_path = args[1] if len(args) > 1 else "chart.dat"
    comprel_path = args[2] if len(args) > 2 else "comprel.log"
    dedup_fort154(in_path, out_path, comprel_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
