#!/usr/bin/env python3
import re
import sys

# 分隔符：只包含星号的行
SEP_RE = re.compile(r'^\s*\*+\s*$')


def is_sep(line: str) -> bool:
    return bool(SEP_RE.match(line))


def extract_name_from_block(block_lines):
    # 找到块内第一行非空行
    first = ''
    for ln in block_lines:
        if ln.strip():
            first = ln
            break
    # 从 "kname = ... Fractional" 提取名字（不区分大小写）
    m = re.search(r'kname\s*=\s*(.*?)\s*Fractional', first, re.IGNORECASE)
    return (m.group(1).strip() if m else first.strip())


def extract_frac_tokens_from_block(block_lines):
    # 从块的首行（或首个包含 Fractional coordinate 的行）提取三分量字符串
    line = ''
    for ln in block_lines:
        if 'Fractional coordinate' in ln:
            line = ln
            break
    if not line:
        # 回退：尝试首个非空行
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
    # 纯浮点
    if re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', tok):
        try:
            return float(tok)
        except Exception:
            return None
    # 分数 a/b
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
    # 仅字母（可含数字但首字符为字母），例如 x, u, k1 -> 归一化为大写
    if re.fullmatch(r'[A-Za-z][A-Za-z0-9]*', tok.strip()):
        return tok.strip().upper()
    return None


def next_point_in_manifold(prev_tokens, curr_tokens, tol=0.005):
    # 线性判定：把表达式近似到两位小数并消元检查是否存在矛盾。
    # 支持项：常数/分数、num*var、var*num、纯变量、以及这些项的加减组合：
    # 例如 1-u, 1-2u, 1-2*u, 0.5-2*u, u+v 等。
    if not prev_tokens or not curr_tokens or len(prev_tokens) != 3 or len(curr_tokens) != 3:
        return None

    def q(x: float) -> float:
        return round(float(x), 2)

    num_re = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|[+-]?\d+/\d+'

    def parse_term(term: str):
        # 返回 (const, {var: coef})，两位小数量化
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
            return None  # var*var 等不支持
        # 纯数字
        val = _parse_number(term)
        if val is not None:
            return q(val), {}
        # 形如 2u 或 1/2u
        m = re.fullmatch(rf'({num_re})?([A-Za-z][A-Za-z0-9]*)', term)
        if m:
            coef_s, var = m.group(1), m.group(2)
            coef = 1.0 if coef_s is None or coef_s == '' else _parse_number(coef_s)
            if coef is None:
                return None
            return 0.0, {var.upper(): q(coef)}
        # 纯变量
        v = _pure_var(term)
        if v is not None:
            return 0.0, {v: 1.0}
        return None

    def parse_linear(expr: str):
        # expr -> (const, dict[var]=coef). 两位小数量化，不支持括号/乘变量。
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

    # 构造 3 条等式：prev[i] - curr[i] = 0
    rows = []
    for i in range(3):
        c1, m1, _ = parse_linear(prev_tokens[i])
        c2, m2, _ = parse_linear(curr_tokens[i])
        row = {'c': q(c1 - c2)}
        for k, v in m1.items():
            row[k] = row.get(k, 0.0) + v
        for k, v in m2.items():
            row[k] = row.get(k, 0.0) - v
        # 清理接近 0
        for k in [kk for kk in list(row.keys()) if kk != 'c']:
            if q(row[k]) == 0.0:
                del row[k]
        row['c'] = q(row['c'])
        rows.append(row)

    # 近似高斯消元
    def eliminate(rows):
        rows = [dict(r) for r in rows]
        used = set()
        while True:
            piv_idx = -1
            piv_var = None
            # 选主元变量
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
            # 规范化主元行
            piv_coef = rows[piv_idx][piv_var]
            for k in list(rows[piv_idx].keys()):
                rows[piv_idx][k] = q(rows[piv_idx][k] / piv_coef)
            # 消去其他行
            for j, rj in enumerate(rows):
                if j == piv_idx:
                    continue
                if piv_var in rj and abs(rj[piv_var]) >= tol:
                    f = rj[piv_var]
                    for k, v in rows[piv_idx].items():
                        rj[k] = q(rj.get(k, 0.0) - f * v)
                    if piv_var in rj:
                        del rj[piv_var]
                # 清理接近 0
                for k in [kk for kk in list(rj.keys()) if kk != 'c']:
                    if q(rj[k]) == 0.0:
                        del rj[k]
                rj['c'] = q(rj.get('c', 0.0))
        return rows

    red = eliminate(rows)
    # 若出现 0 = 非零，则矛盾
    for r in red:
        if all(k == 'c' for k in r.keys()) or len([k for k in r.keys() if k != 'c']) == 0:
            if abs(r.get('c', 0.0)) >= tol:
                return False
    return True


def dedup_fort154(in_path='fort.154', out_path='chart.dat', comprel_path='comprel.log'):
    with open(in_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    # 解析为块：每两行星号之间为一块
    blocks = []  # [(start_sep, block_lines, end_sep)]
    start_sep = None
    cur_lines = None

    for line in lines:
        if is_sep(line):
            if start_sep is None:
                start_sep = line
                cur_lines = []
            else:
                # 结束当前块
                blocks.append((start_sep, cur_lines or [], line))
                # 此星号作为下一块的起始分隔符
                start_sep = line
                cur_lines = []
        else:
            if cur_lines is not None:
                cur_lines.append(line)
            else:
                # 在首个分隔符前的内容忽略
                pass

    # 若最后未以分隔符结束且有内容，则补上一块
    if start_sep is not None and cur_lines is not None and cur_lines:
        blocks.append((start_sep, cur_lines, ''))

    # 去掉重复 kname 的块（全局去重，忽略大小写与多余空白），并记录顺序名称与坐标表达
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

    # 解析 comprel.log，收集有向块元数据（left_base -> right_base）
    cblocks_meta = []
    try:
        with open(comprel_path, 'r', encoding='utf-8', errors='ignore') as f:
            clines = f.readlines()
    except FileNotFoundError:
        clines = []

    if clines:
        eq_re = re.compile(r'^\s*=+\s*$')

        def base_from_side(s: str) -> str:
            # 提取侧字符串的首个字母串作为基名（忽略数字、+、* 等），如 "GM1"->"GM"，"SM1+SM2"->"SM"
            m = re.search(r'[A-Za-z]+', s)
            return m.group(0).upper() if m else ''

        cblocks = []  # list of list[str] (content between '=' lines)
        cur = None
        for ln in clines:
            if eq_re.match(ln):
                if cur is None:
                    cur = []
                else:
                    # 结束当前块
                    cblocks.append(cur)
                    cur = []
            else:
                if cur is not None:
                    cur.append(ln)
        # 若文件未以'='结尾且有内容，补最后一块
        if cur:
            cblocks.append(cur)

        # 提取每个块首行的两个 k 点基名（有向）
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

    # 判断“下一个 k 点是否在本块的集合上”，并把 comprel.log 中匹配的块追加到对应的 chart 块（k1 的块）
    def norm_name_for_chart(name: str) -> str:
        # 与 comprel 的归一化对齐：忽略数字与空白，仅字母，大写
        s = re.sub(r'\d+', '', name)
        s = re.sub(r'\s+', '', s)
        s = re.sub(r'[^A-Za-z]', '', s)
        return s.upper()

    # 计算每个块是否属于上一个块的集合（从第二块起）
    belongs_prev = [None] * len(kept)
    for i in range(1, len(kept)):
        prev_tokens = kept_frac_tokens[i - 1]
        curr_tokens = kept_frac_tokens[i]
        verdict = next_point_in_manifold(prev_tokens, curr_tokens)
        belongs_prev[i] = verdict

    # 收集末尾输出的相容性关系，按有向对分组并保持出现顺序
    compat_groups = {}  # (src_base, dst_base) -> list[list[str]]
    groups_order = []   # list of (src_base, dst_base)
    seen_blocks = set() # 去重

    for i in range(len(kept)):
        if i == 0:
            continue  # 第一块无上一个 k 点
        k1 = norm_name_for_chart(kept_names[i])
        k2 = norm_name_for_chart(kept_names[i - 1])
        if belongs_prev[i] is True:
            # 收集与无序对 {k2,k1} 匹配的 comprel 块（按 comprel 的方向分组）
            for m in cblocks_meta:
                if {m['left'], m['right']} == {k1, k2}:
                    key_dir = (m['left'], m['right'])
                    if key_dir not in compat_groups:
                        compat_groups[key_dir] = []
                        groups_order.append(key_dir)
                    b = m['lines']
                    sig = (key_dir[0], key_dir[1], ''.join(b))
                    if sig in seen_blocks:
                        continue
                    seen_blocks.add(sig)
                    compat_groups[key_dir].append(b)

    # 写出结果
    with open(out_path, 'w', encoding='utf-8') as out:
        for i, (start, blines, end) in enumerate(kept):
            # 不输出星号分隔行，只输出块内容
            out.writelines(blines)
            # 块之间空两行（最后一块除外）
            if i != len(kept) - 1:
                out.write('\n\n')
        # 文件末尾追加相容性关系，按组输出标题
        if groups_order:
            if kept:
                out.write('\n\n')
            for src, dst in groups_order:
                # 对该组内容逐行去重（忽略空行、忽略多余空白，保留首次出现）
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
