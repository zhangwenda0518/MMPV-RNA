# 模型位置

16 科 deploy_genome 模型已移至仓库级数据目录：

```
database/eve_ml/
├── models/        16 科 .npz 模型
├── index.json     OOF AUC + Youden J 阈值
└── FEATURE_SPEC.md  1,030 维特征规格
```

打分命令中的 `--models` 和 `--index` 参数指向该目录：

```bash
python3 eve_score_genome.py \
    --candidates cand.fa \
    --meta labels.tsv \
    --models ../../../database/eve_ml/models \
    --index  ../../../database/eve_ml/index.json \
    --out scores.tsv
```
