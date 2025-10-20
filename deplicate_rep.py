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

    # 去掉与上一块同名的块，并记录顺序名称
    prev_name = None
    kept = []            # list of tuples (start, blines, end)
    kept_names = []      # corresponding kpoint names
    for start, blines, end in blocks:
        name = extract_name_from_block(blines)
        if prev_name is not None and name == prev_name:
            # 丢弃该块
            pass
        else:
            kept.append((start, blines, end))
            kept_names.append(name)
        prev_name = name

    # 解析 comprel.log，构建无序对到块内容的映射
    pair_to_blocks = {}
    try:
        with open(comprel_path, 'r', encoding='utf-8', errors='ignore') as f:
            clines = f.readlines()
    except FileNotFoundError:
        clines = []

    if clines:
        eq_re = re.compile(r'^\s*=+\s*$')

        def norm_kname(s: str) -> str:
            # 忽略数字与空白，仅保留字母，转为大写
            s = s.strip()
            s = re.sub(r'\d+', '', s)
            s = re.sub(r'\s+', '', s)
            s = re.sub(r'[^A-Za-z]', '', s)
            return s.upper()

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

        # 提取每个块首行的两个 k 点（忽略数字）
        for b in cblocks:
            first_nonempty = ''
            for ln in b:
                if ln.strip():
                    first_nonempty = ln.strip()
                    break
            if not first_nonempty or '->' not in first_nonempty:
                continue
            left, right = first_nonempty.split('->', 1)
            k3 = norm_kname(left)
            k4 = norm_kname(right)
            if not k3 or not k4:
                continue
            key = frozenset((k3, k4))
            pair_to_blocks.setdefault(key, []).append(b)

    # 将 comprel.log 中匹配的块追加到对应的 chart 块（k1 的块）
    def norm_name_for_chart(name: str) -> str:
        # 与 comprel 的归一化对齐：忽略数字与空白，仅字母，大写
        s = re.sub(r'\d+', '', name)
        s = re.sub(r'\s+', '', s)
        s = re.sub(r'[^A-Za-z]', '', s)
        return s.upper()

    for i in range(len(kept)):
        if i == 0:
            continue  # 第一块无上一个 k 点
        k1 = norm_name_for_chart(kept_names[i])
        k2 = norm_name_for_chart(kept_names[i - 1])
        key = frozenset((k1, k2))
        if key in pair_to_blocks:
            # 在该块末尾追加匹配的 comprel 块（内容本身，不含 '=' 分隔线）
            start, blines, end = kept[i]
            # 确保上一行以换行结束，并插入标题
            if blines and blines[-1] and not blines[-1].endswith('\n'):
                blines[-1] = blines[-1] + '\n'
            blines.append('Compatibility relations:\n')
            for b in pair_to_blocks[key]:
                blines.extend(b)
                # 块间再加一个空行以便区分
                if not blines or blines[-1].strip():
                    blines.append('\n')
            kept[i] = (start, blines, end)

    # 写出结果
    with open(out_path, 'w', encoding='utf-8') as out:
        for i, (start, blines, end) in enumerate(kept):
            if i == 0:
                out.write(start)
            out.writelines(blines)
            if end:
                out.write(end)


if __name__ == '__main__':
    in_path = sys.argv[1] if len(sys.argv) > 1 else 'fort.154'
    out_path = sys.argv[2] if len(sys.argv) > 2 else 'chart.dat'
    comprel_path = sys.argv[3] if len(sys.argv) > 3 else 'comprel.log'
    dedup_fort154(in_path, out_path, comprel_path)
