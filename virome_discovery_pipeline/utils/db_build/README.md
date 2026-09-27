# db_build/ — 自写建库辅助脚本收编目录

> 这些脚本是 MMPV 自写的建库工具，此前只存在于服务器（`~/bin`、`~/database/cdd/`），仓库里没有任何版本化，丢失后无法再生。
> 本目录用于收编它们。`build_virus_db.py` 的 `check_dependencies()` 会**自动把本目录插到 PATH 最前**——脚本拷进来即可用，服务器 PATH 作为兜底。

## 待收编清单

| 脚本 | 服务器源路径 | 调用方 | 功能 |
|---|---|---|---|
| `SearchAccessionIdToTaxId.py` | `~/bin/` | build_virus_db.py（check_dependencies 强制要求） | accession → taxid 映射 |
| `db_seqid2taxid_add_legth.py` | `~/bin/` | build_virus_db.py | seqid2taxid 增补序列长度 |
| `make_ktaxonomy.py` | `~/bin/` | build_virus_db.py | 生成 ktaxonomy 分类映射文件 |
| `make_genus-length.py` | `~/bin/` | DATABASE_SETUP.md §5（genus_lens 生成） | 从 Assembly Reports 提取病毒属平均基因组长度 → `~/database/virus-db/genus_lens` |
| `build_cdd_databases.sh` | `~/database/cdd/` | 生成 CDD 全套：mmseqs `cdd_db` + 白名单 `cdd_virus_final_v4.txt` + 分类表 `cdd_classified_taxid.tsv`（filter_virus.py / gen_final_judgement.py 硬依赖后两者） | CDD 主库与病毒域白名单/分类表构建 |

## 收编方法（在服务器上执行）

```bash
cd /path/to/MMPV-RNA/virome_discovery_pipeline/utils/db_build/
scp 17711@<server>:~/bin/SearchAccessionIdToTaxId.py ./
scp 17711@<server>:~/bin/db_seqid2taxid_add_legth.py ./
scp 17711@<server>:~/bin/make_ktaxonomy.py ./
scp 17711@<server>:~/bin/make_genus-length.py ./
scp 17711@<server>:~/database/cdd/build_cdd_databases.sh ./
chmod +x *.py *.sh
```

## 同批建议拷回仓库的 CDD 产物（防单点）

```bash
mkdir -p database/cdd
scp 17711@<server>:~/database/cdd/cdd_virus_final_v4.txt  database/cdd/   # 278KB, 可直接入库
scp 17711@<server>:~/database/cdd/cdd_classified_taxid.tsv database/cdd/  # 47MB, 入库或注明服务器为唯一源
```

> 注意 `.gitignore` 当前为 `database/*` 仅白名单 `final.cluster.ref.fasta`，入库需加 `!database/cdd/` 白名单或 `git add -f`。
