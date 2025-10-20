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


def dedup_fort154(in_path="fort.154", out_path="chart.dat"):
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

    # Deduplicate consecutive blocks with the same extracted name
    prev_name = None
    kept = []
    for start, blines, end in blocks:
        name = _extract_name_from_block(blines)
        if prev_name is not None and name == prev_name:
            # Drop this block
            pass
        else:
            kept.append((start, blines, end))
        prev_name = name

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
    dedup_fort154(in_path, out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
