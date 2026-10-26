# -*- coding: utf-8 -*-
"""修复 virus_vcf_pipeline.py 的 Mantel 样本 ID 匹配缺陷。

问题: compute_geo_distance() 要求样本名与元数据索引完全相等，
但管线样本名形如 <Run>_<Acc>，元数据索引只有 <Run> -> 恒不命中，
Mantel 检验（Isolation by Distance）被静默跳过。

修复: 精确匹配优先，未命中则按 SRA run 号匹配。向后兼容。
"""
import io
import os

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_REPO   = f"{_SERVER}/MMPV-RNA"


P = f"{_REPO}/virome_analysis_pipeline/virus_vcf_pipeline.py"
src = io.open(P, encoding='utf-8').read()

# ---------- 1) 确保 import re ----------
if '\nimport re\n' not in src:
    anchor = '\nimport shutil\n'
    assert src.count(anchor) == 1, 'anchor import shutil not unique'
    src = src.replace(anchor, anchor + 'import re\n', 1)
    print('OK: added import re')
else:
    print('SKIP: import re already present')

# ---------- 2) 替换 compute_geo_distance 的匹配逻辑 ----------
OLD = '''def compute_geo_distance(meta_df, samples):
    """
    \u4ece\u5143\u6570\u636e\u8ba1\u7b97\u6837\u672c\u95f4\u5730\u7406\u8ddd\u79bb\u3002
    meta_df: \u542b Sample, Latitude, Longitude (\u6216 Province/Country)
    samples: \u6837\u672c\u540d\u5217\u8868
    \u8fd4\u56de: (n\u00d7n \u5730\u7406\u8ddd\u79bb\u77e9\u9635, labels)
    """
    common = [s for s in samples if s in meta_df.index]
    if len(common) < 5:
        return None, None
    lat = meta_df.loc[common, "Latitude"].values.astype(float)
    lon = meta_df.loc[common, "Longitude"].values.astype(float)
'''

NEW = '''_RUN_ID_RE = re.compile(r"^(?:sample_)?((?:SRR|CRR|ERR|DRR)\\d+)(?=[._-]|$)")


def compute_geo_distance(meta_df, samples):
    """
    \u4ece\u5143\u6570\u636e\u8ba1\u7b97\u6837\u672c\u95f4\u5730\u7406\u8ddd\u79bb\u3002
    meta_df: \u542b Sample, Latitude, Longitude (\u6216 Province/Country)
    samples: \u6837\u672c\u540d\u5217\u8868
    \u8fd4\u56de: (n\u00d7n \u5730\u7406\u8ddd\u79bb\u77e9\u9635, labels)

    \u6837\u672c\u540d\u5339\u914d: \u5148\u6309\u7d22\u5f15\u7cbe\u786e\u5339\u914d; \u672a\u547d\u4e2d\u5219\u6309 SRA run \u53f7\u5339\u914d
    (\u7ba1\u7ebf\u6837\u672c\u540d\u5f62\u5982 <Run>_<Acc>, \u800c\u5143\u6570\u636e\u7d22\u5f15\u901a\u5e38\u53ea\u6709 <Run>)\u3002
    2026-10-09 barbarum: \u4fee\u590d\u524d\u6052\u4e0d\u547d\u4e2d, Mantel \u68c0\u9a8c\u88ab\u9759\u9ed8\u8df3\u8fc7\u3002
    """
    coord, run_to_key = {}, {}
    for _key, _row in meta_df.iterrows():
        _k = str(_key)
        if _k in coord:
            continue
        try:
            _la = float(_row["Latitude"])
            _lo = float(_row["Longitude"])
        except Exception:
            continue
        if not (np.isfinite(_la) and np.isfinite(_lo)):
            continue
        coord[_k] = (_la, _lo)
        _m = _RUN_ID_RE.match(_k)
        if _m:
            run_to_key.setdefault(_m.group(1), _k)

    pairs = []
    for s in samples:
        _k = str(s)
        if _k not in coord:
            _m = _RUN_ID_RE.match(_k)
            if _m:
                _k = run_to_key.get(_m.group(1))
        if _k is not None and _k in coord:
            pairs.append((s, coord[_k]))
    if len(pairs) < 5:
        return None, None
    common = [s for s, _ in pairs]
    lat = np.array([p[0] for _, p in pairs], dtype=float)
    lon = np.array([p[1] for _, p in pairs], dtype=float)
'''

cnt = src.count(OLD)
assert cnt == 1, 'OLD anchor count = %d (expect 1)' % cnt
src = src.replace(OLD, NEW, 1)
io.open(P, 'w', encoding='utf-8').write(src)
print('OK: compute_geo_distance patched')
