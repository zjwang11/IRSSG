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
        # 形如 u/2
        m = re.fullmatch(rf'([A-Za-z][A-Za-z0-9]*)/({num_re})', term)
        if m:
            var, denom = m.group(1), m.group(2)
            dval = _parse_number(denom)
            if dval in (None, 0.0):
                return None
            return 0.0, {var.upper(): q(1.0 / dval)}
        return None

    def parse_expr(expr: str):
        # 拆成带符号的项后累加，得到 (const, {var: coef})
        s = expr.replace(' ', '')
        # 用 '+' 分割，先把负号变成 '+-'
        s = s.replace('-', '+-')
        parts = s.split('+')
        const = 0.0
        coefs = {}
        for p in parts:
            p = p.strip()
            if not p:
                continue
            t = parse_term(p)
            if t is None:
                return None
            c0, vc = t
            const = q(const + c0)
            for v, c in vc.items():
                coefs[v] = q(coefs.get(v, 0.0) + c)
        return const, coefs

    # 解析三分量表达式
    exprs = []  # [(const, {var:coef}), ...]
    vars_all = set()
    for tok in prev_tokens:
        pc = parse_expr(tok)
        if pc is None:
            return None
        const, coefs = pc
        exprs.append((const, coefs))
        vars_all.update(coefs.keys())

    # 下一个点必须是数值坐标
    tgt = []
    for t in curr_tokens:
        v = _parse_number(t)
        if v is None:
            return None
        tgt.append(q(v))

    # 若没有变量，则直接数值比较
    if not vars_all:
        for (c0, coefs), val in zip(exprs, tgt):
            if coefs:
                return None
            if q(c0) != q(val):
                return False
        return True

    # 变量不超过 3 个
    if len(vars_all) > 3:
        return None
    vars_list = sorted(vars_all)

    # 构造增广矩阵 [A|b]，用浮点近似高斯消元判断是否有解
    A = []
    b = []
    for (c0, coefs), val in zip(exprs, tgt):
        row = [q(coefs.get(v, 0.0)) for v in vars_list]
        A.append(row)
        b.append(q(val - c0))

    # 高斯消元（容忍 1e-6 级别误差）
    m = len(A)
    n = len(vars_list)
    # 构造增广矩阵
    M = [row + [bb] for row, bb in zip(A, b)]
    r = 0
    eps = 1e-6
    for c in range(n):
        # 选主元
        pivot = None
        for i in range(r, m):
            if abs(M[i][c]) > eps:
                pivot = i
                break
        if pivot is None:
            continue
        if pivot != r:
            M[r], M[pivot] = M[pivot], M[r]
        pc = M[r][c]
        M[r] = [x / pc for x in M[r]]
        for i in range(m):
            if i == r:
                continue
            factor = M[i][c]
            if abs(factor) > eps:
                M[i] = [xi - factor * xr for xi, xr in zip(M[i], M[r])]
        r += 1
        if r == m:
            break
    # 判矛盾：系数全 0 但 RHS 非 0
    for i in range(m):
        if all(abs(M[i][j]) <= eps for j in range(n)) and abs(M[i][n]) > eps:
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

    # 去掉与上一块同名的块，并记录顺序名称、frac tokens 和规范名
    prev_name = None
    kept = []            # list of tuples (start, blines, end)
    kept_names = []      # matching kpoint names (raw)
    kept_frac_tokens = []
    seen_names = set()
    for start, blines, end in blocks:
        name = extract_name_from_block(blines)
        # 连续相同名字的块去重
        if prev_name is not None and name == prev_name:
            prev_name = name
            continue
        prev_name = name
        # 全局同名的首个保留，其余去重
        nname = re.sub(r'\s+', '', name).upper()
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
        # 收集与无序对 {k2,k1} 匹配的 comprel 块（按 comprel 的方向分组）
        for m in cblocks_meta:
            if {m['left'], m['right']} == {k1, k2}:
                key_dir = (m['left'], m['right'])
                if key_dir not in compat_groups:
                    compat_groups[key_dir] = []
                    groups_order.append(key_dir)
                b = m['lines']
                # 去重仅按整个块内容，避免误删重复行
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
                blocks = compat_groups[(src, dst)]
                if not blocks:
                    continue
                out.write(f"{src}->{dst} Compatibility relations:\n")
                for b in blocks:
                    for ln in b:
                        out.write(ln if ln.endswith('\n') else ln + '\n')
                    # separate blocks by a blank line
                    out.write('\n')


if __name__ == '__main__':
    in_path = sys.argv[1] if len(sys.argv) > 1 else 'fort.154'
    out_path = sys.argv[2] if len(sys.argv) > 2 else 'chart.dat'
    comprel_path = sys.argv[3] if len(sys.argv) > 3 else 'comprel.log'
    dedup_fort154(in_path, out_path, comprel_path)
