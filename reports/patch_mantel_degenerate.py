# -*- coding: utf-8 -*-
"""Mantel 退化保护补丁（补丁 2）。

场景: 当所有样本来自同一坐标（地理距离矩阵全零）时，Mantel r 数学上无定义
(nan)，但 mantel_test() 仍返回一个看似"极显著"的 p（如 0.0001）。
若原样写入 mantel_test.tsv，会被误读为显著 IBD；同时 plot_mantel() 的
np.polyfit 会抛 "SVD did not converge"。

修复: 判定退化（r 非有限 / 地理距离零方差 / 无有效配对）时，
写入 r=nan, p=nan 并跳过作图与趋势线。
"""
import io
import os

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_REPO   = f"{_SERVER}/MMPV-RNA"


P = f"{_REPO}/virome_analysis_pipeline/virus_vcf_pipeline.py"
src = io.open(P, encoding='utf-8').read()

OLD = '''                                if r_val is not None:
                                    print(f"[Mantel] r={r_val:.4f}, p={p_val:.6f} (permutations={n_perm}), N={len(common_samples)}")
                                    mantel_out = os.path.join(stats_dir, "mantel_test.tsv")
                                    with open(mantel_out, "w") as mf:
                                        mf.write(f"r\\tp_value\\tpermutations\\tn_samples\\tn_pairs\\n")
                                        iu = np.triu_indices(len(common_samples), k=1)
                                        g_flat = D_gen_sub[iu]; e_flat = D_geo_sub[iu]
                                        mask = np.isfinite(g_flat) & np.isfinite(e_flat)
                                        mf.write(f"{r_val:.6f}\\t{p_val:.6f}\\t{n_perm}\\t{len(common_samples)}\\t{mask.sum()}\\n")
                                    # scatter plot
                                    if HAS_VIZ:
                                        plot_mantel(g_flat[mask], e_flat[mask], r_val, p_val,
                                                    os.path.join(fig_dir, "mantel_ibd.png"))
                                else:
                                    print(f"[Mantel] skipped: {p_val}")
'''

NEW = '''                                if r_val is not None:
                                    print(f"[Mantel] r={r_val:.4f}, p={p_val:.6f} (permutations={n_perm}), N={len(common_samples)}")
                                    mantel_out = os.path.join(stats_dir, "mantel_test.tsv")
                                    iu = np.triu_indices(len(common_samples), k=1)
                                    g_flat = D_gen_sub[iu]; e_flat = D_geo_sub[iu]
                                    mask = np.isfinite(g_flat) & np.isfinite(e_flat)
                                    # \u9000\u5316\u4fdd\u62a4: r \u65e0\u5b9a\u4e49 / \u5730\u7406\u8ddd\u79bb\u96f6\u65b9\u5dee / \u65e0\u6709\u6548\u914d\u5bf9 -> p \u503c\u65e0\u610f\u4e49, \u4e0d\u5f97\u5448\u73b0\u4e3a\u663e\u8457
                                    _degenerate = False
                                    if (not np.isfinite(r_val)) or mask.sum() == 0:
                                        _degenerate = True
                                    else:
                                        _e = e_flat[mask]
                                        if float(np.nanmax(_e) - np.nanmin(_e)) == 0.0:
                                            _degenerate = True
                                    with open(mantel_out, "w") as mf:
                                        mf.write(f"r\\tp_value\\tpermutations\\tn_samples\\tn_pairs\\n")
                                        if _degenerate:
                                            mf.write(f"nan\\tnan\\t{n_perm}\\t{len(common_samples)}\\t{mask.sum()}\\n")
                                        else:
                                            mf.write(f"{r_val:.6f}\\t{p_val:.6f}\\t{n_perm}\\t{len(common_samples)}\\t{mask.sum()}\\n")
                                    if _degenerate:
                                        print("[Mantel] skipped: degenerate distance matrix (zero variance) -> r undefined")
                                    elif HAS_VIZ:
                                        plot_mantel(g_flat[mask], e_flat[mask], r_val, p_val,
                                                    os.path.join(fig_dir, "mantel_ibd.png"))
                                else:
                                    print(f"[Mantel] skipped: {p_val}")
'''

cnt = src.count(OLD)
assert cnt == 1, 'OLD anchor count = %d (expect 1)' % cnt
src = src.replace(OLD, NEW, 1)
io.open(P, 'w', encoding='utf-8').write(src)
print('OK: mantel degenerate guard patched')
