# 辅助脚本 / Utility Scripts

> preprocess.py, viroid_circular_detect.py — independent helpers for data preprocessing and viroid detection.

## 图件四格式导出 — `mmpv_common/figure_export.py` / `.R`

所有出图脚本默认同时产出 **SVG + PNG(300 dpi) + PDF + PNG(1200 dpi)**，
高清副本落在同目录 `_hires_1200dpi/` 下。既有 `fig.savefig()` / `ggsave()` 写法
无需改动（脚本头已注入引导块），新代码用 `export_figure()` / `save_figure()`。

完整口径、环境变量与未纳入范围见 **[FIGURE_EXPORT.md](FIGURE_EXPORT.md)**。

```bash
# 关闭自动补齐（只写脚本原本要的那一种格式）
MMPV_FIG_EXPORT=0 python <script>.py

# 只产出 PNG + PDF
MMPV_FIG_FORMATS=png,pdf python <script>.py
```

---

## preprocess.py — 预处理独立入口 / Combined Preprocessing

```bash
python data_preprocessing_pipeline/preprocess.py --input <dir> --output <dir> --threads 40 --jobs 10
```

合并 clean-data.py + host_depletion.py 的独立入口（三脚本同位于 `data_preprocessing_pipeline/`）。

---

## viroid_circular_detect.py — 类病毒环状检测

独立的类病毒 (circular RNA) 检测脚本。

```bash
python viroid_circular_detect.py -i <fasta> -o out/ [-t 8]
```

---

## classify_contigs.py — ICTV 参考宿主查找

由 run_host_prediction.py 内部调用, 查找 ICTV 官方分类库中的宿主信息。

```bash
python classify_contigs.py \
    -i taxonomy.tsv -o out/ --prob_dir cross_analysis/
```

参数:

| 参数 | 默认 | 说明 |
|------|------|------|
| `-i` | 必需 | 分类 TSV |
| `-o` | 必需 | 输出目录 |
| `--prob_dir` | cross_analysis/ | 宿主概率表目录 |
