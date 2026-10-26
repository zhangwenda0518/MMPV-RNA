# 部署特征规格（修复版，dim=1030）

与旧版 `deploy_final/` 的差别：删除 `coding_density`（与 `stop_density`
完全共线）、删除 `orf_max_frac`（在 ORF 模式阴性侧按定义≈1）与 `len_nt`
（配对设计下逐对等长、train 内为常数）。

## 特征向量顺序

下标 0–5：

0. `gc`
1. `gc_skew`
2. `stop_density`
3. `atg_frac`
4. `at_frac`
5. `dinuc_entropy`

下标 6–1029：`5-mer` 频率，按 A=0,C=1,G=2,T=3、权重 4^(k-1-i) 编码，
含 N 等非 ACGT 的 k-mer 丢弃后归一化（与 `eve_genome_primary.kmer_vec` 一致）。

## 评分公式

```python
z = (x - mu) / sd - mean_A
score = z @ w + ym
```

`mu`/`sd`/`mean_A`/`w`/`ym` 全部存在 npz 里，`dim` 校验长度。

## 阈值

`index.json` 的 `threshold` 是 OOF 上 Youden 最优点，**不是**真实先验下的
最佳切点；部署时若先验稀疏，应下调或改用排序用法。
