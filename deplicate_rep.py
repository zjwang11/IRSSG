#!/usr/bin/env python3
import re
import sys

# Separator: a line that contains only asterisks
SEP_RE = re.compile(r'^\s*\*+\s*$')


def is_sep(line: str) -> bool:
    return bool(SEP_RE.match(line))


def extract_name_from_block(block_lines):
    # Find the first non-empty line within the block
    first = ''
    for ln in block_lines:
        if ln.strip():
            first = ln
            break
    # Extract name from "kname = ... Fractional" (case-insensitive)
    m = re.search(r'kname\s*=\s*(.*?)\s*Fractional', first, re.IGNORECASE)
    return (m.group(1).strip() if m else first.strip())


def extract_frac_tokens_from_block(block_lines):
    # From the block header (or the first line containing "Fractional coordinate"),
    # extract the three coordinate tokens
    line = ''
    for ln in block_lines:
        if 'Fractional coordinate' in ln:
            line = ln
            break
    if not line:
        # Fallback: try the first non-empty line
        for ln in block_lines:
            if ln.strip():
                line = ln
                break
    if not line:
        return None
    m = re.search(r'Fractional\s+coordinate:\s*(.*?)\s*\(given\s+in\s+the\s+conventional\s+basis\)', line, re.IGNORECASE)
    if not m:
        return None
    inner = m.group(1).strip()
    toks = inner.split()
    if len(toks) != 3:
        return None
    return toks


def _parse_number(tok):
    tok = tok.strip()
    # Pure floating point number
    if re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', tok):
        try:
            return float(tok)
        except Exception:
            return None
    # Fraction a/b
    if re.fullmatch(r'[+-]?\d+/\d+', tok):
        sgn = -1.0 if tok.startswith('-') else 1.0
        if tok[0] in '+-':
            tok2 = tok[1:]
        else:
            tok2 = tok
        a, b = tok2.split('/')
        try:
            return sgn * (float(a) / float(b))
        except Exception:
            return None
    return None


def _pure_var(tok):
    # Letters-only token (alphanumeric allowed but must start with a letter),
    # e.g., x, u, k1 -> normalize to uppercase
    if re.fullmatch(r'[A-Za-z][A-Za-z0-9]*', tok.strip()):
        return tok.strip().upper()
    return None


def next_point_in_manifold(prev_tokens, curr_tokens, tol=0.005):
    # Linear check: round terms to 2 decimals and eliminate to test for inconsistency.
    # Supported terms: constants/fractions, num*var, var*num, pure variables,
    # and sums/differences of these (e.g., 1-u, 1-2u, 1-2*u, 0.5-2*u, u+v).
    if not prev_tokens or not curr_tokens or len(prev_tokens) != 3 or len(curr_tokens) != 3:
        return None

    def q(x: float) -> float:
        return round(float(x), 2)

    num_re = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|[+-]?\d+/\d+'

    def parse_term(term: str):
        # Return (const, {var: coef}) with values quantized to 2 decimals
        term = term.strip()
        if term == '':
            return 0.0, {}
        if '*' in term:
            parts = term.split('*')
            if len(parts) != 2:
                return None
            a, b = parts[0].strip(), parts[1].strip()
            aval = _parse_number(a)
            bval = _parse_number(b)
            avar = _pure_var(a)
            bvar = _pure_var(b)
            if aval is not None and bvar is not None:
                return 0.0, {bvar: q(aval)}
            if bval is not None and avar is not None:
                return 0.0, {avar: q(bval)}
            if aval is not None and bval is not None:
                return q(aval * bval), {}
            return None  # var*var not supported
        # Pure numeric literal
        val = _parse_number(term)
        if val is not None:
            return q(val), {}
        # Forms like 2u or 1/2u
        m = re.fullmatch(rf'({num_re})?([A-Za-z][A-Za-z0-9]*)', term)
        if m:
            coef_s, var = m.group(1), m.group(2)
            coef = 1.0 if coef_s is None or coef_s == '' else _parse_number(coef_s)
            if coef is None:
                return None
            return 0.0, {var.upper(): q(coef)}
        # Pure variable
        v = _pure_var(term)
        if v is not None:
            return 0.0, {v: 1.0}
        return None

    def parse_linear(expr: str):
        # expr -> (const, dict[var]=coef). Quantized to 2 decimals; no parentheses
        # or var*var.
        if expr is None:
            return 0.0, {}, False
        s = expr.strip()
        if s == '':
            return 0.0, {}, True
        s = s.replace('−', '-')
        s = s.replace(' ', '')
        terms = re.findall(r'[+-]?[^+-]+', s)
        const = 0.0
        coeffs = {}
        ok = True
        for t in terms:
            if t in ('+', '-'):
                continue
            parsed = parse_term(t)
            if parsed is None:
                ok = False
                tag = f'EXPR_{hash(t) & 0xffff:04x}'
                coeffs[tag] = coeffs.get(tag, 0.0) + 1.0
                continue
            c, m = parsed
            const += c
            for k, v in m.items():
                coeffs[k] = coeffs.get(k, 0.0) + v
        const = q(const)
        coeffs = {k: q(v) for k, v in coeffs.items() if q(v) != 0.0}
        return const, coeffs, ok

    # Build 3 equations: prev[i] - curr[i] = 0
    rows = []
    for i in range(3):
        c1, m1, _ = parse_linear(prev_tokens[i])
        c2, m2, _ = parse_linear(curr_tokens[i])
        row = {'c': q(c1 - c2)}
        for k, v in m1.items():
            row[k] = row.get(k, 0.0) + v
        for k, v in m2.items():
            row[k] = row.get(k, 0.0) - v
        # Drop near-zero terms
        for k in [kk for kk in list(row.keys()) if kk != 'c']:
            if q(row[k]) == 0.0:
                del row[k]
        row['c'] = q(row['c'])
        rows.append(row)

    # Approximate Gaussian elimination
    def eliminate(rows):
        rows = [dict(r) for r in rows]
        used = set()
        while True:
            piv_idx = -1
            piv_var = None
            # Choose pivot variable
            for i, r in enumerate(rows):
                vars_in_r = [k for k in r.keys() if k != 'c']
                if not vars_in_r:
                    continue
                for v in vars_in_r:
                    if v not in used and abs(r[v]) >= tol:
                        piv_idx = i
                        piv_var = v
                        break
                if piv_idx != -1:
                    break
            if piv_idx == -1:
                break
            used.add(piv_var)
            # Normalize pivot row
            piv_coef = rows[piv_idx][piv_var]
            for k in list(rows[piv_idx].keys()):
                rows[piv_idx][k] = q(rows[piv_idx][k] / piv_coef)
            # Eliminate other rows
            for j, rj in enumerate(rows):
                if j == piv_idx:
                    continue
                if piv_var in rj and abs(rj[piv_var]) >= tol:
                    f = rj[piv_var]
                    for k, v in rows[piv_idx].items():
                        rj[k] = q(rj.get(k, 0.0) - f * v)
                    if piv_var in rj:
                        del rj[piv_var]
                # Drop near-zero terms
                for k in [kk for kk in list(rj.keys()) if kk != 'c']:
                    if q(rj[k]) == 0.0:
                        del rj[k]
                rj['c'] = q(rj.get('c', 0.0))
        return rows

    red = eliminate(rows)
    # If we get 0 = nonzero, it's inconsistent
    for r in red:
        if all(k == 'c' for k in r.keys()) or len([k for k in r.keys() if k != 'c']) == 0:
            if abs(r.get('c', 0.0)) >= tol:
                return False
    return True


def dedup_fort154(in_path='fort.154', out_path='chart.dat', comprel_path='comprel.log'):
    with open(in_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    # Parse into blocks: each pair of asterisk lines forms one block
    blocks = []  # [(start_sep, block_lines, end_sep)]
    start_sep = None
    cur_lines = None

    for line in lines:
        if is_sep(line):
            if start_sep is None:
                start_sep = line
                cur_lines = []
            else:
                # End current block
                blocks.append((start_sep, cur_lines or [], line))
                # This asterisk line starts the next block
                start_sep = line
                cur_lines = []
        else:
            if cur_lines is not None:
                cur_lines.append(line)
            else:
                # Ignore content before the first separator
                pass

    # If the file doesn't end with a separator and has content, add the final block
    if start_sep is not None and cur_lines is not None and cur_lines:
        blocks.append((start_sep, cur_lines, ''))

    # Drop duplicate kname blocks globally (case/whitespace-insensitive),
    # and record order/names/coordinates
    def norm_block_name(s: str) -> str:
        return re.sub(r"\s+", "", s).upper()

    seen_names = set()
    kept = []            # list of tuples (start, blines, end)
    kept_names = []      # corresponding kpoint names
    kept_frac_tokens = []  # corresponding fractional coordinate tokens (list[str] or None)
    for start, blines, end in blocks:
        name = extract_name_from_block(blines)
        nname = norm_block_name(name)
        if nname in seen_names:
            continue
        seen_names.add(nname)
        kept.append((start, blines, end))
        kept_names.append(name)
        kept_frac_tokens.append(extract_frac_tokens_from_block(blines))

    # Parse comprel.log and collect directed block metadata (left_base -> right_base)
    cblocks_meta = []
    try:
        with open(comprel_path, 'r', encoding='utf-8', errors='ignore') as f:
            clines = f.readlines()
    except FileNotFoundError:
        clines = []

    if clines:
        eq_re = re.compile(r'^\s*=+\s*$')

        def base_from_side(s: str) -> str:
            # Extract first letter sequence as base name (ignore digits,+,*),
            # e.g., "GM1"->"GM", "SM1+SM2"->"SM"
            m = re.search(r'[A-Za-z]+', s)
            return m.group(0).upper() if m else ''

        cblocks = []  # list of list[str] (content between '=' lines)
        cur = None
        for ln in clines:
            if eq_re.match(ln):
                if cur is None:
                    cur = []
                else:
                    # End current block
                    cblocks.append(cur)
                    cur = []
            else:
                if cur is not None:
                    cur.append(ln)
        # If the file doesn't end with '=' and has content, add the last block
        if cur:
            cblocks.append(cur)

        # Extract two k-point base names (directed) from the first line of each block
        for b in cblocks:
            first_nonempty = ''
            for ln in b:
                if ln.strip():
                    first_nonempty = ln.strip()
                    break
            if not first_nonempty or '->' not in first_nonempty:
                continue
            left, right = first_nonempty.split('->', 1)
            lbase = base_from_side(left)
            rbase = base_from_side(right)
            if lbase and rbase:
                cblocks_meta.append({'left': lbase, 'right': rbase, 'lines': b})

    # Decide whether the next k-point lies on the current block's manifold,
    # and append matching comprel.log blocks to the corresponding chart block (k1)
    def norm_name_for_chart(name: str) -> str:
        # Normalize aligned with comprel: drop digits/whitespace, letters only, uppercase
        s = re.sub(r'\d+', '', name)
        s = re.sub(r'\s+', '', s)
        s = re.sub(r'[^A-Za-z]', '', s)
        return s.upper()

    # For each block starting from the second, decide if it belongs to the
    # previous block's manifold
    belongs_prev = [None] * len(kept)
    for i in range(1, len(kept)):
        prev_tokens = kept_frac_tokens[i - 1]
        curr_tokens = kept_frac_tokens[i]
        verdict = next_point_in_manifold(prev_tokens, curr_tokens)
        belongs_prev[i] = verdict

    # Collect compatibility relations for output, grouped by unordered pairs
    # in appearance order. Note: {k1,k2} and {k2,k1} are the same group; key
    # uses sorted (a,b).
    compat_groups = {}  # (a,b) -> list[list[str]]
    groups_order = []   # list of (a,b)
    seen_blocks = set()  # Deduplicate

    for i in range(len(kept)):
        if i == 0:
            continue  # First block has no previous k-point
        k1 = norm_name_for_chart(kept_names[i])
        k2 = norm_name_for_chart(kept_names[i - 1])
        if belongs_prev[i] is True:
            # Collect comprel blocks matching unordered pair {k2,k1}
            # (grouped by unordered pair)
            for m in cblocks_meta:
                if {m['left'], m['right']} == {k1, k2}:
                    key_pair = tuple(sorted((m['left'], m['right'])))
                    if key_pair not in compat_groups:
                        compat_groups[key_pair] = []
                        groups_order.append(key_pair)
                    b = m['lines']
                    sig = (key_pair[0], key_pair[1], re.sub(r'\s+', ' ', ''.join(b)).strip())
                    if sig in seen_blocks:
                        continue
                    seen_blocks.add(sig)
                    compat_groups[key_pair].append(b)

    # Write results
    with open(out_path, 'w', encoding='utf-8') as out:
        for i, (start, blines, end) in enumerate(kept):
            # Do not output separator lines; only block content
            out.writelines(blines)
            # Two blank lines between blocks (except after the last block)
            if i != len(kept) - 1:
                out.write('\n\n')
        # Append compatibility relations at the end, with group headers
        if groups_order:
            if kept:
                out.write('\n\n')
            for src, dst in groups_order:
                # Deduplicate lines within the group (ignore blank lines,
                # extra whitespace; keep first occurrence)
                seen = set()
                unique_lines = []
                for b in compat_groups[(src, dst)]:
                    for ln in b:
                        key = re.sub(r'\s+', ' ', ln).strip()
                        if not key:
                            continue
                        if key in seen:
                            continue
                        seen.add(key)
                        unique_lines.append(ln if ln.endswith('\n') else ln + '\n')
                if unique_lines:
                    out.write(f"{src}->{dst} Compatibility relations:\n")
                    out.writelines(unique_lines)
                    out.write('\n')


if __name__ == '__main__':
    in_path = sys.argv[1] if len(sys.argv) > 1 else 'fort.154'
    out_path = sys.argv[2] if len(sys.argv) > 2 else 'chart.dat'
    comprel_path = sys.argv[3] if len(sys.argv) > 3 else 'comprel.log'
    dedup_fort154(in_path, out_path, comprel_path)
