# 辅助脚本 / Utility Scripts

> preprocess.py, viroid_circular_detect.py — independent helpers for data preprocessing and viroid detection.

## preprocess.py — 预处理独立入口 / Combined Preprocessing

```bash
python preprocess.py --input <dir> --output <dir> --threads 40 --jobs 10
```

合并 clean-data.py + host_depletion.py 的独立入口。

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
