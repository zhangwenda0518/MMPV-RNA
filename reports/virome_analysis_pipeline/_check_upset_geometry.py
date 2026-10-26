#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""_check_upset_geometry.py — 分类 UpSet 图「几何体检」

本图（`Fig_GojiVirome_UpsetClassification[_Family]`）已迭代 5 次，每次都要回答同一类问题，
故固化成判据（**全部量像素，不看"出图成功"**）：

  P1 柱↔点 **逐列对齐残差 = 0 px**
     （patchwork 保证；但改主题 / 改宽高比会破）
  P2 点阵**绝对 x 位置 + 列间距**（判「点阵是否被整体平移 / 拉伸」）
  P3 **矩阵行标签槽（x≈820..1040）内深色文字段数 = 0**
     （2026-10-14 需求：集合名只保留最左 set-size 面板那一处）
  P4 画布四边各 4 px **零墨迹**（无 TEXT-OFF-CANVAS）

用法:
  python _check_upset_geometry.py <fig.png> [<参考图.png>]
     带参考图时额外打印**逐列差异剖面**（直接指出「哪些 x 区间变了」）。

⚠⚠ 必须喂 `--no-timing` 的**干净底图**：左上耗时条会污染「饱和色柱」检测
   （实测把耗时条蓝柱当成 x541 的"交集柱"）。脚本会自检并报错。
⚠ 文字带**次像素抗锯齿**（暖色边缘）→ 用「棕色」掩膜找点阵时会被标签文字误命中，
   故点阵检测限定 x≥960 且段宽 20..45 px。
"""
import sys
import numpy as np
from PIL import Image

# 面板 x 分区（画布 2796 x 1800 @300dpi，9.32 x 6.0 in）
# ⚠⚠ 2026-10-14 第 7 次修订后**布局变了**：set-size 面板加宽（`width_ratio` 0.26 → 0.312）
#    → 面板 x 254..904（原 254..757），矩阵面板**不动**（1054..2488），中间留 150 px 空白
#    （= 顶部柱面板 y 轴刻度文字 + 两侧 margin，结构性，压不掉）。
#    故 `X_SETSIZE_END` / `X_TAGSLOT` 必须同步改，否则 P3 会把 set-size 的**灰柱**当文字误报。
X_LABEL_END = 260      # 最左集合名标签区: 0..260
X_SETSIZE_END = 906    # set-size 灰条面板: 260..906（第 7 次修订后）
X_MATRIX_END = 2480    # 矩阵点阵面板: 1054..2480
# 矩阵「行标签槽」原位置：第 6 次修订时是 x≈820..1040 的一串透明文字；
# 第 7 次修订改 `element_blank()` 后**槽已归零** → 该处变成纯空白（904..1054 的 gutter）。
X_TAGSLOT = (912, 1048)
Y_TOP = (350, 1050)    # 交集柱行
Y_BOT = (1150, 1780)   # 点阵 / set-size 行


def _segs(idx, gap):
    if len(idx) == 0:
        return []
    out, s = [], idx[0]
    for i in range(1, len(idx)):
        if idx[i] - idx[i - 1] > gap:
            out.append((s, idx[i - 1])); s = idx[i]
    out.append((s, idx[-1]))
    return out


def profile(p, tag):
    a = np.asarray(Image.open(p).convert("RGB")).astype(np.int16)
    H, W, _ = a.shape
    R, G, B = a[..., 0], a[..., 1], a[..., 2]
    print(f"\n===== {tag}  {W}x{H} =====")

    # 自检：左上是否残留耗时条（x<950 且 y350..1050 有饱和色）
    sat_all = (a.max(axis=2) - a.min(axis=2)) > 60
    z = np.zeros_like(sat_all); z[Y_TOP[0]:Y_TOP[1], :950] = True
    n_overlay = int((sat_all & z).sum())
    if n_overlay > 200:
        print(f"  ⚠⚠ 检测到左上残留 {n_overlay} 个饱和像素 —— 这不是 --no-timing 干净底图，"
              f"P1/P2 不可信！")

    # 交集柱（饱和色，x 限定在柱面板内）
    # ⚠⚠ 下界必须 ≥1050：顶部柱图的 **y 轴刻度文字**（"5,000"…"0"）在 x≈998..1028，
    #    黑字经 **ClearType 子像素抗锯齿**会产出 (204,139,95)/(95,139,204) 这类**彩边**
    #    → `max-min>60` 判据会把它当"饱和柱"，凭空多出一根假柱（实测 x≈1010，零报错）。
    z = np.zeros_like(sat_all); z[Y_TOP[0]:Y_TOP[1], 1050:X_MATRIX_END] = True
    cs = (sat_all & z).sum(axis=0)
    bars = _segs(np.where(cs >= 3)[0], 8)
    bars = [s for s in bars if s[1] - s[0] >= 20]
    bcen = [(s0 + s1) // 2 for s0, s1 in bars]
    print(f"  P1 交集柱 {len(bcen)} 列；前 5 中心 = {bcen[:5]}")

    # 点阵棕点（限定矩阵面板 + 段宽 20..45）
    dot = (R > 150) & (G > 40) & (G < 150) & (B < 95) & ((R - B) > 80)
    z = np.zeros_like(dot); z[Y_BOT[0]:Y_BOT[1], 960:X_MATRIX_END] = True
    cd = (dot & z).sum(axis=0)
    dots = [s for s in _segs(np.where(cd >= 5)[0], 8) if 20 <= s[1] - s[0] <= 45]
    dcen = [(s0 + s1) // 2 for s0, s1 in dots]
    print(f"  P2 点阵 {len(dcen)} 列；前 5 中心 = {dcen[:5]}")
    if len(dcen) >= 2:
        sp = [dcen[i + 1] - dcen[i] for i in range(len(dcen) - 1)]
        print(f"     首列 x = {dcen[0]}  末列 x = {dcen[-1]}  列间距 {min(sp)}..{max(sp)} px")
    else:
        sp = []

    # 柱↔点 逐列对齐
    if len(bcen) == len(dcen) and len(bcen):
        res = [abs(b - d) for b, d in zip(bcen, dcen)]
        ok = "PASS" if max(res) == 0 else f"FAIL(最大 {max(res)} px)"
        print(f"  P1 柱↔点 对齐残差 max = {max(res)} px  → {ok}  ({len(res)} 列)")
    else:
        print(f"  P1 柱 {len(bcen)} 列 vs 点 {len(dcen)} 列 —— 数量不等，无法逐列对拍")

    # 矩阵「行标签槽」内是否还有文字
    dark = (R < 130) & (G < 130) & (B < 130)
    z = np.zeros_like(dark); z[Y_BOT[0]:Y_BOT[1], X_TAGSLOT[0]:X_TAGSLOT[1]] = True
    c3 = np.where((dark & z).sum(axis=0) > 0)[0]
    seg3 = [s for s in _segs(c3, 25) if s[1] - s[0] >= 8]
    print(f"  P3 矩阵行标签槽 x{X_TAGSLOT[0]}..{X_TAGSLOT[1]} 内深色文字段 = {len(seg3)} 处 "
          f"→ {'PASS（行标签已只在最左）' if not seg3 else 'FAIL ' + str(seg3)}")

    # 最左标签区
    z = np.zeros_like(dark); z[Y_BOT[0]:Y_BOT[1], :X_LABEL_END] = True
    c4 = np.where((dark & z).sum(axis=0) > 0)[0]
    seg4 = [s for s in _segs(c4, 25) if s[1] - s[0] >= 8]
    print(f"      最左标签区 x0..{X_LABEL_END}: {len(seg4)} 段 "
          + (f"x {seg4[0][0]}..{seg4[-1][1]}" if seg4 else "（空！）"))

    # 四边 4 px 零墨迹
    ink = (a.min(axis=2) < 245)
    bad = {k: int(v.sum()) for k, v in
           {"上": ink[:4, :], "下": ink[-4:, :], "左": ink[:, :4], "右": ink[:, -4:]}.items()
           if v.sum() > 0}
    print(f"  P4 四边 4 px 墨迹: " + ("PASS（全 0）" if not bad else f"FAIL {bad}"))

    # P5 左上「耗时图框」净空（框 = x94..946 / y265..1076，见 fig_upset_classification.R）
    # ⚠ 判据 = **该框内零墨迹**（底图干净才能叠预渲耗时图）；⚠⚠ 改 set-size 面板宽
    #   （`width_ratio`）会让面板右缘右移，**必须重验**框右缘与面板/刻度文字的间距。
    z = np.zeros_like(ink); z[265:1077, 94:947] = True
    n_box = int((ink & z).sum())
    print(f"  P5 耗时图框 x94..946 / y265..1076 内墨迹 = {n_box} → "
          + ("PASS" if n_box == 0 else "FAIL"))
    # 框右缘与柱面板 y 轴刻度文字的左缘间距
    # ⚠ `.sum(axis=0)` 的索引是**画布绝对 x**（不是窗口内偏移）→ 直接用 `c5.min()`。
    z2 = np.zeros_like(ink); z2[265:1077, 947:1050] = True
    c5 = np.where((ink & z2).sum(axis=0) > 0)[0]
    print(f"      框右缘 946 → 刻度文字左缘 {int(c5.min()) if len(c5) else '（无）'}"
          f"  间距 {int(c5.min()) - 946 if len(c5) else '-'} px")
    return a


a_new = profile(sys.argv[1], "TARGET")
if len(sys.argv) > 2:
    a_ref = profile(sys.argv[2], "REFERENCE")
    if a_ref.shape == a_new.shape:
        d = (np.abs(a_ref - a_new).max(axis=2) > 25)
        print(f"\n===== 逐列差异剖面（REF vs TARGET）总差异像素 {int(d.sum())} =====")
        for name, band in [("上半 y<1150（交集柱行）", d[:1150]),
                           ("下半 y>=1150（点阵/set-size 行）", d[1150:])]:
            print(f"  {name}:")
            col = band.sum(axis=0)
            for s0, s1 in _segs(np.where(col > 0)[0], 12):
                if s1 - s0 >= 2:
                    print(f"      x {s0:>5}..{s1:<5} (宽 {s1-s0+1:>4})")
    else:
        print(f"\n⚠ 尺寸不同（{a_ref.shape} vs {a_new.shape}），跳过差异剖面")
