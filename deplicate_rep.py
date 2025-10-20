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


def dedup_fort154(in_path='fort.154', out_path='chart.dat'):
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

    # 去掉与上一块同名的块
    prev_name = None
    kept = []
    for start, blines, end in blocks:
        name = extract_name_from_block(blines)
        if prev_name is not None and name == prev_name:
            # 丢弃该块
            pass
        else:
            kept.append((start, blines, end))
        prev_name = name

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
    dedup_fort154(in_path, out_path)

