"""mjrm_config.py — BEAST2 Markov Jumps & Rewards Matrix 配置生成（自研）

等效 VirPhyKit MJRM Generator (src/MakovMJump/function_mmj.py)，不依赖 VirPhyKit。
用途：为 BEAST2 离散性状分析（Markov Jumps）生成迁移矩阵 + rewards 配置，
     或把矩阵/rewards 插入到已有 BEAST2 XML 模板中。

用法:
  python mjrm_config.py China,USA,Japan out_dir config.txt
  python mjrm_config.py China,USA,Japan out_dir beast_with_mjrm.xml --xml template.xml
"""
from __future__ import annotations

import os


def wrtcfg(*migr: str) -> str:
    """生成迁移矩阵文本（对角线 0，单向 1）。"""
    tlen = len(migr)
    str0 = ['0'] * tlen
    output = []
    for ii in range(tlen):
        for jj in range(tlen):
            if ii == jj:
                continue
            str1 = str0[:]
            str1[jj] = '1'
            output.append('<parameter id="%2s-to-%2s" value="\n' % (migr[ii], migr[jj]))
            for kk in range(ii):
                output.append(' '.join(str0) + '\n')
            output.append(' '.join(str1) + '\n')
            for kk in range(ii + 1, tlen):
                output.append(' '.join(str0) + '\n')
            output.append('"/>\n')
    return ''.join(output)


def wrt_rewards(*migr: str) -> str:
    """生成 rewards 配置（每性状 1.0/0.0 向量）。"""
    tlen = len(migr)
    output = ['<rewards>\n']
    for ii in range(tlen):
        values = ['0.0'] * tlen
        values[ii] = '1.0'
        output.append(f'    <parameter id="{migr[ii]}_reward" value="{" ".join(values)}" />\n')
    output.append('</rewards>\n')
    return ''.join(output)


_TARGET_MARKERS = (
    "<!--  END Ancestral state reconstruction",
    "<!-- END Ancestral state reconstruction",
    "</markovJumpsTreeLikelihood>",
)


def process_xml(input_xml_path: str, migr: list, output_xml_path: str):
    """把矩阵 + rewards 插入 XML 模板（在祖先状态重建标记前）。"""
    try:
        with open(input_xml_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except UnicodeDecodeError:
        return False, 'Error: XML file encoding is not UTF-8'
    except OSError as e:
        return False, f'Error: {e}'

    comment_index = -1
    for i, line in enumerate(lines):
        if any(marker in line for marker in _TARGET_MARKERS):
            comment_index = i
            break
    if comment_index == -1:
        return False, ('Error: Target markers "END Ancestral state reconstruction" '
                       'or "</markovJumpsTreeLikelihood>" not found in XML file')

    matrix_content = wrtcfg(*migr)
    rewards_content = wrt_rewards(*migr)
    output_lines = lines[:comment_index] + [matrix_content, rewards_content] + lines[comment_index:]
    with open(output_xml_path, 'w', encoding='utf-8') as f:
        f.writelines(output_lines)
    return True, f'Success: Modified XML with matrix and rewards saved to {output_xml_path}'


def generate_config(migr_text: str, file_path: str, filename: str,
                    input_xml_path: str = None):
    """生成 config 文件（.txt）或处理 XML（给 input_xml_path）。

    返回 (success: bool, message: str)
    """
    if not migr_text:
        return False, 'Error: No migr parameters provided'
    if not file_path:
        return False, 'Error: No save directory provided'
    if not filename:
        return False, 'Error: No filename provided'
    try:
        if ',' not in migr_text and len(migr_text.strip().split()) > 1:
            return False, ('Error: The Migr Parameters input is incorrect. '
                           'Please separate different parameters with English commas.')
        migr = [x.strip() for x in migr_text.split(',') if x.strip()]
        if not migr or len(migr) < 2:
            return False, 'Error: No valid migr parameters provided or fewer than 2 parameters'
        output_file = os.path.join(file_path.rstrip('/\\'), filename)
        if input_xml_path:
            if not output_file.endswith('.xml'):
                output_file += '.xml'
            return process_xml(input_xml_path, migr, output_file)
        if not output_file.endswith('.txt'):
            output_file += '.txt'
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(wrtcfg(*migr))
            f.write(wrt_rewards(*migr))
        return True, f'Success: File with matrix and rewards saved to {output_file}'
    except Exception as e:
        return False, f'Error: {str(e)}'


def main():
    import argparse
    ap = argparse.ArgumentParser(description='BEAST2 Markov Jumps & Rewards 配置生成')
    ap.add_argument('traits', help='逗号分隔离散性状, 如 China,USA,Japan')
    ap.add_argument('output_dir', help='输出目录')
    ap.add_argument('filename', help='输出文件名')
    ap.add_argument('--xml', default=None, help='BEAST2 XML 模板路径 (给则插入矩阵到 XML)')
    args = ap.parse_args()
    ok, msg = generate_config(args.traits, args.output_dir, args.filename, args.xml)
    print(msg)
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())