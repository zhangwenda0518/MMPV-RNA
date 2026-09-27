# data_preprocessing_pipeline — 数据清洗管线（质控 + 去宿主）

**定位**: 链路第 2 步。质控清理 → 去宿主去 rRNA，输出可直接进病毒发现/已知病毒分析的 clean reads。
脚本独立、检查点续传；支持一条龙或分步运行。

## 脚本

| 脚本 | 角色 |
|---|---|
| `clean-data.py` | 数据清理: Fastp (QC+去接头) → Seqkit (FASTQ→FASTA+统计) → Clumpify (k-mer 聚类重排，提升压缩率、加速组装)，自动识别 PE/SE |
| `host_depletion.py` | 去宿主 + 去rRNA: Kraken2 → 比对 (Bowtie2/HISAT2/Minimap2) → Ribodetector/SILVA |
| `data_preprocessing.py` | 一条龙便捷入口（可配置版）: 清洗 → 去宿主 → 报告, `--stage clean/deplete/bbnorm/report/all` |
| `preprocess.py` | 一条龙便捷入口（轻量版）: 清洗 → 去宿主 |
| `run_bbnorm.py` | 可选覆盖度归一化: BBNorm (target=70, mindepth=2)，自动 PE/SE；由 `--stage bbnorm` 显式调用（`--stage all` 不含），输出 `00c_BBnorm/` |
| `preprocess_report.py` | report 阶段: 汇总各步统计 → 交接清单（纯解析上游已落盘的统计文件） |

## 用法

```bash
# 质控
python data_preprocessing_pipeline/clean-data.py -i raw_fastqs/ -o out/00a_CleanData/ -j 20 -t 16

# 去宿主 + 去rRNA (Ribodetector, 仅 rna-short)
python data_preprocessing_pipeline/host_depletion.py \
    -k out/host_reference/hostdb/kraken2 \
    -x out/host_reference/hostdb/bowtie2/host \
    -I out/00a_CleanData/ -O out/00b_HostDepletion/ \
    --tool bowtie2 --seq-type rna-short --rrna -t 40 -j 10

# 一条龙: 清洗 → 去宿主 → 报告 (--host_db 自动推导 kraken2/比对索引子目录)
python data_preprocessing_pipeline/data_preprocessing.py --stage all --input_reads raw/ --output_dir out/ --host_db out/host_reference/hostdb/ -t 40 -j 10

# 单独补跑 report (汇总 + 交接清单)
python data_preprocessing_pipeline/data_preprocessing.py --stage report --output_dir out/

# 可选: 覆盖度归一化 (深度差异大的多样本共组装前才需要; --stage all 不含此步)
# 输入为 deplete 产物 00b_HostDepletion/, 输出 00c_BBnorm/ 可直接作 --input_reads
python data_preprocessing_pipeline/data_preprocessing.py --stage bbnorm --output_dir out/ -t 16 -j 4
```

输出流转: clean 产出 `1.fastp/ → 2.fasta/ → 3.clumpify/`，deplete 输入优先取 `3.clumpify/`，回退 `2.fasta/`；
bbnorm（可选）输入取 `00b_HostDepletion/`，输出 `00c_BBnorm/`。

## 宿主数据库: 三种指定方式（deplete 阶段）

| 方式 | 参数 | 说明 |
|---|---|---|
| ① 现成库目录 | `--host_db <dir>` | 自动查找子索引；`<dir>` 支持三种指向：hostdb/ 本身、host_reference 根（自动下钻）、直接的 Kraken2 库目录（含 hash.k2d） |
| ② 分开精确指定 | `--kraken2_db <dir> --host_align_db <prefix>` | Kraken2 库与比对索引前缀分别给定 |
| ③ 现场构建 | `--host_fasta <FASTA>` 或 `--host_species X --host_taxid N` | 无现成库时直接调 `download_host_genome.py`（四通道）+ `build_hostbase.py` 构建到 `{output_dir}/host_reference/`（只建 kraken2 + 当前比对工具；复用已有结果，`--force` 重建） |

`--host_db` 兼容上游产物: 公共数据管道 `--stage hostref` 建好的 `public_data_pipeline_output/host_reference/` 可直接作为 `--host_db` 传入。

## report 阶段交付物（上下游交接）

`--stage all` / `--stage report` 在输出根目录产出三件套：

| 文件 | 用途 |
|---|---|
| `preprocessing_summary.tsv` | 每样本 Raw → Fastp → Kraken2 → 去宿主 → 去rRNA 各步 reads 数、留存率、比对率、PASS/LOW_RETAINED 状态（阈值 `--warn_retained`，默认 20%） |
| `assembly_ready.list` | 最终 reads 文件绝对路径清单（逐行一条）→ 交棒 `virome_pipeline.py --input_reads` |
| `preprocessing_report.html` | 交互式 HTML（留存率图 + 明细表，Chart.js 内嵌离线可用） |

统计源（纯解析、缺失留空不报错）：clean 的 fastp JSON；deplete 的 `host_depletion_seqkit_summary.tsv`、`ribodetector.report.txt`、bowtie2/hisat2 stderr 日志（`overall alignment rate`）。

## 在整体链路中的位置

```
public_metadata_pipeline/  →  data_preprocessing_pipeline/  →  virome_discovery_pipeline/
   (测序数据检索/下载/转换       (质控 → 去宿主+去rRNA)          (病毒发现, de novo)
    + 宿主参考基因组获取/建库)              ↘ (已知病毒轨: 仅质控、免去宿主)
```

宿主参考库（Kraken2/Bowtie2/HISAT2/Minimap2 四索引）由上游 `public_metadata_pipeline/build_host_pipeline.py`
构建（NCBI 四通道下载 datasets/ngd/FTP/gget，或 `--genome-fasta` 自有基因组），详见其 `doc.md` 第 6-8 节。

**上下游衔接**:
- **上游**: `public_metadata_pipeline/sra2fastx.py` 产出的 FASTQ.GZ（其 report 阶段的 `sample_handoff.csv` 即本管道文件级输入清单）；宿主索引来自 `build_host_pipeline.py` 的 `hostdb/`
- **下游**: depleted reads 作为 `virome_discovery_pipeline/virome_pipeline.py --input_reads` 进入组装；
  `virome_pipeline.py` 的 `clean`/`deplete` 阶段自动调用本目录两个主脚本；
  已知病毒分析轨（auto_known_virus）直接以 clean 输出为输入；
  `assembly_ready.list` / `preprocessing_summary.tsv` 即人工核对与下游取数的交接点
