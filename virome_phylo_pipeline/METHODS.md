# 系统发育与分子进化分析（virome_phylo_pipeline）—— 论文级方法章

> **📌 文档定位：本文件是「论文里怎么写」的权威**（科学表述、阈值、参考文献）。
> 工具调用与参数名以本文件为准；stage 内部细节见 `STAGE_REFERENCE.md`；
> 「怎么跑」见 `RUN_GUIDE.md`；文档总索引见 `README.md`。
>
> 本方法章可直接粘贴进投稿 Methods；工具名、参数、阈值均取自 `virome_phylo_pipeline/` 脚本实际调用。
> 软件版本号标 ⊙ 者待从服务器运行环境确认。文末附参考文献与占位符清单。
>
> **同步状态（2026-09-16）**：§1 建树支持值、§4 地理解析器已按代码更新；
> §2 的 DRT 随机化次数与本文件末「投稿前需确认」条目一并列出。改动前备份
> `_docconv_20260916/backup/METHODS.md`。

---

## 1. 序列比对与系统发育树构建

对每个高置信病毒，将全长基因组序列用 MAFFT v7.525（Katoh & Standley, 2013）做多序列比对（`--auto --reorder`）。比对结果用 IQ-TREE v3.1.2（Minh et al., 2020）构建最大似然树，模型经 ModelFinder 自动选择（`-m MFP`；Kalyaanamoorthy et al., 2017），分支支持度以 ultrafast bootstrap 1000 次重复（`-B 1000`）联合 SH-aLRT 1000 次重复（`-alrt 1000`）评估；两类支持值均写入 IQ-TREE 树标签，并另存 `iqtree.support.json` 记录实际生效的支持值类型（若 UFBoot 失败会降级为仅 SH-aLRT 并在该文件标注 `fallback`，报告不将降级结果表述为双支持值）。比对前对低覆盖序列质控（coverage ≥ 85%），剔除低深度组装产生的 40–85% 覆盖嵌合序列，避免虚高差异。本地与公共数据库合并后，按（采样日期、采样地点、序列完全一致）三元组去除克隆重复序列，组内保留一条代表性序列；不同时间或地点采得的相同序列予以保留，以维持跨时空传播信号。

## 2. 时间信号验证（定年前置检验）

定年之前先检验序列是否携带充分的时间信号：

- **Root-to-tip 回归**：用 TreeTime（Sagulenko et al., 2018）将各序列到根节点的遗传距离对采样日期做线性回归，计算决定系数 R² 与进化速率 β。分级阈值：R² > 0.85 可用严格分子钟，> 0.70 推荐放松分子钟，> 0.10 为边际信号（须先做 DRT 确认），> 0.05 为弱信号，< 0.01 视为无时间信号。负 β 提示生根方向错误，用 TempEst（Rambaut et al., 2016）复核。
- **日期随机化检验（DRT）**（Duchêne et al., 2015）：将采样日期随机打乱并重跑 RTT 回归（`--drt_randomizations`，管线默认 20；本研究取 100），真实数据的 R² 需超过随机分布 95% 百分位方判定通过（`passed = percentile ≥ 95`）。
- **定年决策**：以 DRT 经验 p 值为准（百分位 ≥ 95 即通过，见 `utils/temporal_signal.py`），通过者进行 BEAST 定年并报告 TMRCA；R² 只作参考（经验参考线 R² > 0.3，非硬门槛）。DRT 未通过则不报 TMRCA（短序列类病毒因位点独立突变与弱时间信号，分子钟不可靠）。

> **方法学说明（实际执行与阴性结论）**：本节的时间信号分级阈值先行用于分基因筛选，随后以 DRT 作正式判定；DRT 未通过者不进入 BEAST 定年、不报告 TMRCA，相关基因在论文中仅以拓扑与分化特征陈述。对枸杞相关病毒（GCVA）六基因分段的实测中，6 个基因段的真实 R² 为 0.010–0.050，DRT 百分位 70–95，仅 N 基因段恰好达到 95 阈值，其余 5 段判定为未通过（`DRT_TT3_ALL_GENES.tsv`）。据此本研究不对该数据集作分子钟定年解读，而以遗传分化与空间扩散描述其进化格局。
>
> **关于 BETS**：贝叶斯时间信号评估（Duchêne et al., 2020）曾作为候选检验纳入流程，因其在 1.x 平台需注入 marginal likelihood estimator 的路径采样框架，本研究的实现对复杂树先验不收敛（多轮参数与结构验证后边际似然发散），已整体停用且不纳入本文结果与方法陈述。判定时间信号的依据为 RTT 与 DRT 两项。

## 3. 分子钟定年（BEAST 1.x）

对通过时间信号检验的病毒，用 BEAST v1.10.4（Suchard et al., 2018）配合 BEAGLE 库（Ayres et al., 2019）做贝叶斯定年。分子钟用 uncorrelated lognormal relaxed clock（UCLN）；树先验用 coalescent skyline（默认）、constant 或 birth-death skyline（birth-death skyline 需 classpath 覆盖注册 `BirthDeathSerialSkylineModelParser`）；替代模型用 HKY 或 GTR 加 gamma 位点速率异质性（4 类）。MCMC 以多链独立运行（默认 5 链 × 8 线程；每链 states 数以 config.yaml 为准——2026-09-16 修订：`beast.chain` 默认 5×10⁶、`phylogeo.chain` 默认 1.5×10⁷，此前文档误写 5×10⁷ 比实际高 5–10 倍；**论文引用时以该轮运行的 run_status.json 记录为准，不写统一值**），不同随机种子初始化，burnin 去除前 10% states 后，LogCombiner 合并各链。链混合度用 Geyer (1992) initial positive sequence 估计量计算 ESS（与 Tracer/beast-mcmc TraceCorrelation 严格一致：有偏自协方差 + 成对和>0 截断，maxLag 2000；ess_results.csv 由 utils/ess.py 生成，2026-08-27 经 trace.jar 官方库逐参数对账验证，残余偏差<15% 来自 burnin 取整），收敛判据 ESS > 200。TreeAnnotator 生成最大分支可信度（MCC）树与节点 95% HPD 区间。采样日期无法解析的样本（如标注 "not applicable"）在 XML 构建阶段剔除并记录，不参与定年；非编码 RNA 病毒（类病毒）因参考基因组无 CDS，不进行分基因定年与密码子水平选择压力分析。

## 4. 离散系统地理学分析

在 BEAST 1.x 中同时重建祖先地理状态：以采样地点为离散性状，用连续时间马尔可夫链（CTMC）模型加贝叶斯随机搜索变量选择（BSSVS；Lemey et al., 2009）推断地点间迁移率与显著迁移路径。采样地点按实际保留、简化至省/州级，超过 12 个地点时其余归并。迁移路径以 Bayes factor（BF）≥ 5 判定显著。辅以 RRT（随机化组间距离检验地理聚类显著性）与 TempMig（时序迁移追踪），并用 Mantel 检验评估遗传距离与地理距离相关（9999 次置换）。采样地点坐标由分层地理解析器获得：入口治理阶段使用内置全球地名表（Natural Earth 1:10m，共 11960 条：255 个国家、4410 个一级行政区、7295 个城市），支持 `国家:行政区:城市` 多级拆分与中英文及拼音的规范化匹配，中国港台澳条目归一至中国；地理分析阶段的解析器在离线省/市表之外，对未命中地名回退查询 Photon/Nominatim 公共地理编码服务并将结果持久缓存，以消除重复网络请求。两条路径**均只采用实际命中的地名条目，不做坐标插补**；删去 `Unknown` 层级后仅剩国家名的记录判为无地理分辨率，既不生成国家质心等推测坐标，亦不参与地理分析（因坐标不可解析而被排除的样本数在日志中如实记录）。地理距离以 Haversine 球面距离计算。传播路径用 SpreaD3（Bielejec et al., 2016）交互式地图可视化。

## 5. 群体遗传学分析

用 pypopart（PopART 的纯 Python 实现；Leigh & Bryant, 2015）做单倍型与多态性分析：coverage ≥ 85% 质控、单倍型鉴定、多样性指标（π、θ、分离位点数、Tajima's D、Fu's Fs）、单倍型网络（MJN 无向与 VirNA 最小生成网络有向）与群体分化 Fst（位点法加 permutation）。全唯一单倍型（Hd = 1.0）视为保守病毒的特征而非缺陷，星状网络如实呈现；单倍型网络视觉分区与 Fst 显著性不等价。

## 6. 选择压力分析

对蛋白编码序列做密码子水平选择压力检测：DnaSP 6（Rozas et al., 2017）计算 Ka/Ks、Fu & Li's D*/F*、最小重组数 Rm、连锁不平衡与位点频谱；codeml（PAML 4；Yang, 2007）与 HyPhy 2.5（Kosakovsky Pond et al., 2020）对基因树做位点级与基因级正选择检验（FEL、MEME、BUSTED）。

## 7. 重组检测

用 RDP5（Martin et al., 2021，wine 封装）对多序列比对做重组断点检测（RDP、GENECONV、BootScan、MaxChi、Chimaera、SiScan、3Seq），对检出事件做解析验证、IQ-TREE 拓扑切换验证与 SimPlot 相似度图。比对长度 < 500 bp 的超短基因组（如类病毒）因滑动窗口无足够位点不适用此分析，如实跳过。下游建树采用两种策略处理重组区：将重组区间 mask 为 N，或删除高置信（≥ 3 方法）重组序列。断点分布与基因映射经 RBDP/Rgrapher 可视化，DnaSP Rm 交叉验证。

## 8. 软件与版本

| 工具 | 版本 | 引用 |
|---|---|---|
| MAFFT | 7.525 | Katoh & Standley, 2013 |
| IQ-TREE | 3.1.2 | Minh et al., 2020 |
| ModelFinder | ⊙ | Kalyaanamoorthy et al., 2017 |
| BEAST | 1.10.4 | Suchard et al., 2018 |
| BEAGLE | ⊙ | Ayres et al., 2019 |
| TreeTime | ⊙ | Sagulenko et al., 2018 |
| TempEst | ⊙ | Rambaut et al., 2016 |
| TreeDater | ⊙ | Volz & Frost, 2017 |
| pypopart | ⊙（GitHub） | Leigh & Bryant, 2015 |
| DnaSP | 6 | Rozas et al., 2017 |
| RDP5 | ⊙（wine 封装） | Martin et al., 2021 |
| codeml（PAML） | ⊙ | Yang, 2007 |
| HyPhy | 2.5 | Kosakovsky Pond et al., 2020 |
| VirNA | ⊙（GitHub） | 无正式文献 |
| SpreaD3 | ⊙ | Bielejec et al., 2016 |

分析在 Linux 计算服务器（246）执行，Python 3.10（numpy、pandas、scipy、Biopython）。BEAST 经 `java -Xmx2048m` 调用，并行链数与线程数由 `--chains`/`--threads` 控制。

---

## 参考文献

- Ayres DL, Cummings MP, Baele G, et al. 2019. BEAGLE 3: Improved performance, scaling, and usability for a high-performance computing library for statistical phylogenetics. Systematic Biology 68:1052–1061.
- Bielejec F, Baele G, Vrancken B, Suchard MA, Rambaut A, Lemey P. 2016. SpreaD3: Interactive visualization of spatiotemporal history and trait evolutionary processes. Molecular Biology and Evolution 33:2167–2169.
- Duchêne S, Duchêne D, Holmes EC, Ho SYW. 2015. The performance of the date-randomization test in phylogenetic analyses of time-structured virus data. Molecular Biology and Evolution 32:1895–1906.
- Duchêne S, Lemey P, Stadler T, et al. 2020. Bayesian evaluation of temporal signal in measurably evolving populations. Molecular Biology and Evolution 37:3363–3379.
- Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. Nature Methods 14:587–589.
- Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution 30:772–780.
- Kosakovsky Pond SL, Poon AFY, Velazquez R, et al. 2020. HyPhy 2.5—A customizable platform for evolutionary hypothesis testing using phylogenies. Molecular Biology and Evolution 37:295–299.
- Leigh JW, Bryant D. 2015. PopART: Full-feature software for haplotype network construction. Methods in Ecology and Evolution 6:1110–1116.
- Lemey P, Rambaut A, Drummond AJ, Suchard MA. 2009. Bayesian phylogeography finds its roots. PLoS Computational Biology 5:e1000520.
- Martin DP, Varsani A, Roumagnac P, et al. 2021. RDP5: a computer program for analyzing recombination in, and removing signals of recombination from, nucleotide sequence datasets. Virus Evolution 7:veaa087.
- Minh BQ, Schmidt HA, Chernomor O, et al. 2020. IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution 37:1530–1534.
- Rambaut A, Lam TT, Max Carvalho L, Pybus OG. 2016. Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen). Virus Evolution 2:vew007.
- Rozas J, Ferrer-Mata A, Sánchez-DelBarrio JC, et al. 2017. DnaSP 6: DNA sequence polymorphism analysis of large data sets. Molecular Biology and Evolution 34:3299–3302.
- Sagulenko P, Puller V, Neher RA. 2018. TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evolution 4:vex042.
- Suchard MA, Lemey P, Baele G, Ayres DL, Drummond AJ, Rambaut A. 2018. Bayesian phylogenetic and phylodynamic data integration using BEAST 1.10. Virus Evolution 4:vey016.
- Volz EM, Frost SDW. 2017. Scalable relaxed clock phylogenetic dating. Virus Evolution 3:vex025.
- Yang Z. 2007. PAML 4: phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution 24:1586–1591.

---

## 占位符清单（投稿前需替换 / 核对）

- [ ] 病毒名与序列数（如 "Cytorhabdovirus sp. 'lycii'，n = 186"）
- [ ] 标 ⊙ 的软件版本号
- [ ] 定年数值（R²、进化速率、TMRCA 及 95% HPD）——注：GCVA 数据集 DRT 未通过，此项不报；若其他病毒分析通过则填此表
- [x] MCMC 链长与 burnin 比例（以 config.yaml 默认为基准、run_status.json 实际值为准；默认 beast 5×10⁶ / phylogeo 1.5×10⁷，burnin 10%）
- [ ] 单倍型网络着色维度（location 或 host）
- [ ] 重组事件数量与涉及病毒
- [x] 时间信号阴性结论回填（第 2 节方法学说明已写，2026-09-10）
- [ ] **DRT 随机化次数**：正文写「本研究取 100」，而管线默认 20（`--drt_randomizations`）——
      需与服务器实际运行值核对后统一，勿凭记忆填写
- [ ] **建树支持值口径**：正文写 `-B 1000 -alrt 1000`（双支持值）。
      需按实际运行产出的 `phylogeny/iqtree.support.json` 核对 `support_type`；
      若发生降级（`fallback` / `degraded`），正文须如实改为对应的单一支持值或无支持值
- [ ] **coverage ≥ 85% 序列质控的来源**：该质控原由 `pstvd_qc.py` 承担，**该脚本已移除**；
      现 `clean` stage 只做长度比例/N 占比/大小写/非法字符过滤，**不含 coverage 过滤**。
      需确认该质控实际由上游提取流程（`06_extraction`）完成；若确无来源，本节应改写或删除该句
