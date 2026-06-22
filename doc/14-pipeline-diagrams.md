# MMPV-RNA v2.3 流程框架图

> 基于昨日实际运行命令行历史绘制 | 2026-06-13

---

## 一、virome_pipeline.py 全架构总览

```mermaid
flowchart TB
    accTitle: MMPV-RNA virome_pipeline.py 10阶段全架构
    accDescr: 从 Raw FASTQ 到 HQ vOTU 的 10 阶段端到端流水线架构, 包含数据清洗/去宿主/组装/鉴定/COBRA延伸/聚类/分类/宿主预测/CheckV/拯救

    subgraph Input["输入层"]
        RAW["Raw FASTQ<br/>--input_reads $RAW/"]
    end

    subgraph Preprocess["预处理层"]
        CL["1. Clean<br/>clean-data.py<br/>Fastp → Seqkit → Clumpify<br/>-t 120 -j 20"]
        DP["2. Deplete<br/>host_depletion.py<br/>Kraken2 → Bowtie2 → rRNA<br/>--host_db ~/database/host_db/<br/>--rrna -t 120 -j 20"]
    end

    subgraph Core["核心分析层"]
        AS["3. Assembly<br/>assembly_pipeline.py<br/>MEGAHIT<br/>-t 20 -j 20 -m 256"]
        ID["4. Identification<br/>virus_identification.py<br/>6工具并行鉴定<br/>--virus_db + 全DB路径<br/>-t 30 -j 10"]
        CO["5. COBRA<br/>cobra_pipeline.py<br/>BWA-MEM2→CoverM→COBRA<br/>auto detect 组装工具<br/>-t 30 -j 10"]
        CLU["6. Cluster<br/>cluster_pipeline.py<br/>CD-HIT ref-guide + vclust Leiden<br/>--ref-genomes ×2<br/>-t 60"]
    end

    subgraph Postprocess["后处理层"]
        TX["7. Taxonomy<br/>virus_classifier.py + R<br/>9工具分类 → R共识<br/>--virus_db + 全DB路径<br/>-t 60"]
        HO["8. Host<br/>run_host_prediction.py<br/>ICTV > RNAVirHost > PhaBOX2<br/>-t 120"]
        CV["9. CheckV<br/>checkv completeness<br/>按宿主预评估<br/>--checkv_db checkv-db-v1.7<br/>-t 120"]
        RE["10. Rescue<br/>rescue_pipeline.py<br/>三支路A→C→D级联拯救<br/>--host-filter Plant<br/>-t 20 -j 10"]
    end

    subgraph Output["输出层"]
        VLD["validate_novel_viruses.py<br/>★known / ★★novel / ★★★truly<br/>→ final_virus_catalog.fasta"]
    end

    RAW --> CL --> DP --> AS --> ID --> CO --> CLU
    CLU --> TX --> HO --> CV --> RE --> VLD

    classDef input fill:#e0f2fe,stroke:#0284c7,color:#0c4a6e
    classDef pre fill:#dcfce7,stroke:#16a34a,color:#14532d
    classDef core fill:#fef9c3,stroke:#ca8a04,color:#713f12
    classDef post fill:#fce7f3,stroke:#db2777,color:#831843
    classDef out fill:#ede9fe,stroke:#7c3aed,color:#4c1d95

    class RAW input
    class CL,DP pre
    class AS,ID,CO,CLU core
    class TX,HO,CV,RE post
    class VLD out
```

---

## 二、virome_pipeline.py 实际运行流程 (昨日命令行)

```mermaid
flowchart LR
    accTitle: 昨日实际运行流程图
    accDescr: 展示 2026-06-13 从 clean 到 rescue 的 10 阶段实际执行顺序及关键参数

    subgraph S1["14:36 — Stage: clean"]
        direction TB
        c1["python virome_pipeline.py --stage clean<br/>--input_reads $RAW/ --output_dir $OUT/<br/>-t 120 -j 20"]
        c1_sub["→ clean-data.py<br/>400样本 PE Fastq → Clean FASTA"]
        c1 --> c1_sub
    end

    subgraph S2["15:27 — Stage: deplete"]
        direction TB
        c2["python virome_pipeline.py --stage deplete<br/>--output_dir $OUT/<br/>--host_db ~/database/host_db/<br/>--aligner bowtie2 --seq_type rna-short --rrna<br/>-t 120 -j 20"]
        c2_sub["→ host_depletion.py<br/>Kraken2 + Bowtie2 + Ribodetector"]
        c2 --> c2_sub
    end

    subgraph S3["16:06 — Stage: assembly"]
        direction TB
        c3["python virome_pipeline.py --stage assembly<br/>--output_dir $OUT/<br/>--input_reads $OUT/00b_HostDepletion/<br/>--assembler megahit<br/>-t 20 -j 20 -m 256 -l 0 --force"]
        c3_sub["→ assembly_pipeline.py<br/>MEGAHIT k-list 21..99"]
        c3 --> c3_sub
    end

    subgraph S4["16:48 — Stage: identification"]
        direction TB
        c4["python virome_pipeline.py --stage identification<br/>--output_dir $OUT/ --virus_db ~/database/virus-db<br/>--virus_protein_db ... --uniprot_db ... --viroids_db ...<br/>--virsorter_db ... --viralverify_hmm ... --virhunter_path ...<br/>--virhunter_weights ... --metabuli_db ... --virbot_path ...<br/>--viralm_path ... --virus_taxid ...<br/>--blast_mode both --blast_evalue 1e-5 --blast_top_n 5<br/>-t 30 -j 10"]
        c4_sub["→ virus_identification.py<br/>6工具并行 + Venn图"]
        c4 --> c4_sub
    end

    subgraph S5["18:03 — Stage: cobra"]
        direction TB
        c5["python virome_pipeline.py --stage cobra<br/>--output_dir $OUT/<br/>--input_reads $OUT/00b_HostDepletion/<br/>-t 30 -j 10"]
        c5_sub["→ cobra_pipeline.py<br/>auto detect 组装工具 → BWA-MEM2 + COBRA"]
        c5 --> c5_sub
    end

    subgraph S6["20:34 — Stage: cluster"]
        direction TB
        c6["python virome_pipeline.py --stage cluster<br/>--output_dir $OUT/ -t 30<br/>--ref-genomes final.complete_ref.fasta viral.1.1.genomic.fna"]
        c6_sub["→ cluster_pipeline.py<br/>CD-HIT ref-guide → vclust Leiden → centroids"]
        c6 --> c6_sub
    end

    subgraph S7["20:31 — Stage: taxonomy"]
        direction TB
        c7["python virome_pipeline.py --stage taxonomy<br/>--output_dir $OUT/ --virus_db ~/database/virus-db/<br/>--uniprot_db ... --genomad_db ... --metabuli_db ...<br/>--cat_db ... --cat_tax ... --mmseqs_db ...<br/>--vitap_db ... --acvirus_db ... --vcontact3_db ...<br/>-t 60"]
        c7_sub["→ virus_classifier.py + R consensus<br/>9工具并行 → 8级taxonomy"]
        c7 --> c7_sub
    end

    subgraph S8["22:18 — Stage: host"]
        direction TB
        c8["python virome_pipeline.py --stage host<br/>--output_dir $OUT/<br/>--phabox-db ~/database/virus-db/phabox_db_v2_2/<br/>--prob-dir .../cross_analysis/<br/>-t 120"]
        c8_sub["→ run_host_prediction.py<br/>ICTV > RNAVirHost > PhaBOX2 决策树"]
        c8 --> c8_sub
    end

    subgraph S9["22:27 — Stage: checkv"]
        direction TB
        c9["python virome_pipeline.py --stage checkv<br/>--output_dir $OUT/<br/>--checkv_db ~/database/virus-db/checkv-db-v1.7<br/>-t 120"]
        c9_sub["→ checkv completeness<br/>按宿主分组评估"]
        c9 --> c9_sub
    end

    subgraph S10["23:15 — Stage: rescue"]
        direction TB
        c10["python virome_pipeline.py --stage rescue<br/>--output_dir $OUT/<br/>--input_reads $OUT/00b_HostDepletion/<br/>--checkv_db ~/database/virus-db/checkv-db-v1.7<br/>--blast-db $DB/ref.fasta<br/>--host-filter Plant -t 20 -j 10"]
        c10_sub["→ rescue_pipeline.py<br/>三支路 A→C→D → HQ vOTU"]
        c10 --> c10_sub
    end

    S1 --> S2 --> S3 --> S4 --> S5 --> S6 --> S7 --> S8 --> S9 --> S10

    classDef stage fill:#dbeafe,stroke:#2563eb,color:#1e3a5f
    class S1,S2,S3,S4,S5,S6,S7,S8,S9,S10 stage
```

---

## 三、virome_pipeline.py 参数传递框架

```mermaid
flowchart TB
    accTitle: virome_pipeline.py 参数传递关系
    accDescr: 展示编排器 virome_pipeline.py 如何通过 _build_parser 收集参数并分发到各阶段的子脚本

    subgraph Orchestrator["virome_pipeline.py 编排器"]
        PARSER["_build_parser()<br/>60+ 参数 10 个组"]
        MAIN["main()<br/>--dry-run → scan<br/>--stage → dispatch"]
        PARSER --> MAIN
    end

    subgraph GlobalCtrl["编排器自身控制"]
        G1["--stage all/clean/..."]
        G2["--dry-run"]
        G3["--log-level DEBUG/INFO"]
        G4["--stop-on-error"]
        G5["--force"]
    end

    subgraph StageParams["阶段专用参数组"]
        direction TB
        SP1["Clean<br/>--dedup --no_compress<br/>--clumpify_memory --clean_debug"]
        SP2["Deplete<br/>--kraken2_confidence --keep_rrna<br/>--rrna_chunk_size --deplete_steps<br/>--align_config --deplete_debug"]
        SP3["Assembly<br/>--refinec_threads --refinec_min_id<br/>--refinec_min_cov --asm_tmp_dir<br/>--asm_keep_temp"]
        SP4["Identification<br/>--nr_db --skip_uniprot_filter<br/>--skip_nr_filter --skip_id_plots<br/>--clean_failed --ident_ext<br/>--virus_protein_db --uniprot_db<br/>--viroids_db --virsorter_db<br/>--viralverify_hmm --metabuli_db<br/>--virus_taxid --virhunter_path<br/>--virhunter_weights --virbot_path<br/>--viralm_path --blast_mode"]
        SP5["COBRA<br/>--cobra_mink --cobra_maxk<br/>--cobra_linkage_mismatch<br/>--cobra_verbose --virus_mode"]
        SP6["Cluster<br/>--skip_vclust --vclust_cluster_file<br/>--min-length --ani --qcov<br/>--cdhit_ani --cdhit_qcov<br/>--ref-genomes"]
        SP7["Taxonomy<br/>--tax_jobs --tax_ext<br/>--tax_remove_suffix<br/>--genomad_db --cat_db --cat_tax<br/>--mmseqs_db --vitap_db<br/>--acvirus_db --vcontact3_db"]
        SP8["Host<br/>--skip_rnavirhost --skip_phabox<br/>--skip_ictv --phabox-db<br/>--prob-dir"]
        SP9["Rescue<br/>--host-filter --blast-db<br/>--virseqimprover-path<br/>--salmon-bin"]
    end

    subgraph CommonParams["公共参数"]
        CP["-t/--threads<br/>-m/--memory<br/>-j/--jobs<br/>--force<br/>--input_reads<br/>--output_dir"]
    end

    subgraph Dispatch["阶段调度 → 子脚本"]
        direction LR
        R1["run_clean() → clean-data.py"]
        R2["run_depletion() → host_depletion.py"]
        R3["run_assembly() → assembly_pipeline.py"]
        R4["run_identification() → virus_identification.py"]
        R5["run_cobra() → cobra_pipeline.py"]
        R6["run_cluster() → cluster_pipeline.py"]
        R7["run_taxonomy() → virus_classifier.py + R"]
        R8["run_host() → run_host_prediction.py"]
        R9["run_checkv_stage() → checkv completeness"]
        R10["run_rescue() → rescue_pipeline.py"]
    end

    GlobalCtrl --> MAIN
    SP1 --> R1
    SP2 --> R2
    SP3 --> R3
    SP4 --> R4
    SP5 --> R5
    SP6 --> R6
    SP7 --> R7
    SP8 --> R8
    SP9 --> R10
    CP --> Dispatch

    classDef orch fill:#dbeafe,stroke:#2563eb,color:#1e3a5f
    classDef ctrl fill:#dcfce7,stroke:#16a34a,color:#14532d
    classDef param fill:#fef3c7,stroke:#d97706,color:#78350f
    classDef common fill:#ede9fe,stroke:#7c3aed,color:#4c1d95
    classDef dispatch fill:#fce7f3,stroke:#db2777,color:#831843

    class Orchestrator orch
    class GlobalCtrl ctrl
    class SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9 param
    class CP common
    class R1,R2,R3,R4,R5,R6,R7,R8,R9,R10 dispatch
```

---

## 四、virome_pipeline.py 每个阶段内部调用流程

### 4.1 Clean 阶段

```mermaid
flowchart LR
    accTitle: Clean 阶段内部流程
    accDescr: clean-data.py 内部三步流水线 Fastp 质控 → Seqkit 转FASTA → Clumpify 去重

    subgraph Clean["run_clean()"]
        direction LR
        A1["1. Fastp<br/>--qualified_quality_phred 20<br/>-g --poly_g_min_len 10<br/>--length_required 50"]
        A2["2. Seqkit fq2fa<br/>FASTQ → FASTA<br/>-w 0 (不换行)"]
        A3["3. Clumpify<br/>BBMap 光学去重<br/>--skip-clumpify 跳过"]
        A1 --> A2 --> A3
    end

    IN[/"00a_CleanData/<br/>1.fastp_tmp/ → 2.fasta/ → 3.clumpify/"/]
    OUT[/"reads_dir → 3.clumpify/"/]

    Clean --> IN --> OUT

    classDef clean fill:#dcfce7,stroke:#16a34a,color:#14532d
    class Clean clean
```

### 4.2 Deplete 阶段

```mermaid
flowchart LR
    accTitle: Deplete 阶段内部流程
    accDescr: host_depletion.py 三阶段 Kraken2 初筛 → Bowtie2 精筛 → Ribodetector/SILVA 去rRNA

    subgraph Deplete["run_depletion()"]
        direction LR
        B1["1. Kraken2<br/>--confidence 0.4<br/>物种分类标记宿主 reads"]
        B2["2. Bowtie2<br/>--end-to-end<br/>精确比对去宿主<br/>samtools view -f 12 -F 256"]
        B3["3. rRNA 去除<br/>Ribodetector (默认)<br/>SILVA Bowtie2 (可选)"]
        B1 --> B2 --> B3
    end

    IN2[/"reads_dir<br/>→ 00b_HostDepletion/"/]
    OUT2[/"清洁 reads<br/>*_clean_1.fa.gz<br/>*_clean_2.fa.gz"/]

    Deplete --> IN2 --> OUT2

    classDef dp fill:#dbeafe,stroke:#2563eb,color:#1e3a5f
    class Deplete dp
```

### 4.3 Assembly 阶段

```mermaid
flowchart TB
    accTitle: Assembly 阶段内部流程
    accDescr: assembly_pipeline.py 支持 MEGAHIT/rnaviralSPAdes/Penguin/All 四种组装模式

    subgraph Assembly["run_assembly()"]
        direction TB
        TOOLS["--assembler<br/>megahit | rnaviralspades | penguin | all"]
        AUTO["auto detect 已有输出<br/>→ 跳过已完成样本"]
        ASM["并行组装<br/>ProcessPoolExecutor<br/>-m 256 GB -t 20"]

        TOOLS --> AUTO --> ASM

        subgraph RefineCMode["≥2工具模式"]
            RC1["refineC split<br/>--frag-min-len 1000"]
            RC2["merge_all split"]
            RC3["refineC merge<br/>--min-id 0.95 --min-cov 0.50"]
            RC1 --> RC2 --> RC3
        end
        ASM --> RefineCMode
    end

    OUT3[/"01_Assembly/{sample}/<br/>{sample}_{tool}.contig.fasta"/]
    Assembly --> OUT3

    classDef asm fill:#fef3c7,stroke:#d97706,color:#78350f
    class Assembly asm
```

### 4.4 Identification 阶段

```mermaid
flowchart LR
    accTitle: Identification 阶段内部流程
    accDescr: virus_identification.py 6工具并行鉴定 + UniProt/NR 后置过滤

    subgraph Identify["run_identification()"]
        direction TB
        INPUT["input: 01_Assembly/{sample}/*.contig.fasta"]
        TOOLS6["6工具并行鉴定<br/>--------------------<br/>Genomad · Diamond BLASTX<br/>VirSorter2 · ViralVerify<br/>VirHunter · Metabuli"]
        FILTER["后置过滤<br/>UniProt strict/filter<br/>NR 对抗验证"]
        VENN["Venn/Upset 可视化<br/>--skip_id_plots 跳过"]

        INPUT --> TOOLS6 --> FILTER --> VENN
    end

    OUT4[/"02_Identification/{sample}/<br/>*_virus.all.candidate.fasta"/]
    Identify --> OUT4

    classDef ident fill:#fce7f3,stroke:#db2777,color:#831843
    class Identify ident
```

### 4.5 COBRA 阶段

```mermaid
flowchart LR
    accTitle: COBRA 阶段内部流程
    accDescr: cobra_pipeline.py 自动匹配三元组 → BWA-MEM2 比对 → CoverM 覆盖度 → COBRA 延伸

    subgraph COBRA["run_cobra()"]
        direction LR
        MATCH["自动匹配三元组<br/>reads + contig + virus<br/>--virus-mode strict/filter/raw"]
        BWAMEM["BWA-MEM2 比对<br/>samtools sort"]
        COVERM["CoverM 覆盖度<br/>covered_fraction/mean/rpkm"]
        COBRA_M["COBRA 延伸<br/>--mink 21 --maxk 141<br/>--linkage-mismatch 2"]
        MATCH --> BWAMEM --> COVERM --> COBRA_M
    end

    OUT5[/"03_COBRA/{sample}/<br/>cobra_{tool}_result/*.cobra.fa"/]
    COBRA --> OUT5

    classDef cobra fill:#ede9fe,stroke:#7c3aed,color:#4c1d95
    class COBRA cobra
```

### 4.6 Cluster 阶段

```mermaid
flowchart TB
    accTitle: Cluster 阶段内部流程
    accDescr: cluster_pipeline.py CD-HIT 参考引导预聚类 + vclust Leiden 聚类

    subgraph Cluster["run_cluster()"]
        direction TB
        COLLECT["自动收集 COBRA *.cobra.fa<br/>或 --cluster_input 直接指定"]
        SEQKIT["seqkit 长度过滤<br/>--min-length 500"]
        CDHIT["CD-HIT 参考引导预聚类<br/>vclust dedup → cd-hit<br/>--ref-genomes ICTV+NCBI<br/>拆分 known/novel"]
        VCLUST["vclust Leiden<br/>prefilter → align → cluster<br/>--ani 0.95 --qcov 0.85"]
        CENTROID["centroids 产出<br/>known_linked + vclust novel<br/>→ final_centroids.fasta"]

        COLLECT --> SEQKIT --> CDHIT --> VCLUST --> CENTROID
    end

    OUT6[/"04_CLUSTER/<br/>centroids/final_centroids.fasta<br/>centroids/known_association.tsv<br/>3_vclust/vclust_clusters.tsv"/]
    Cluster --> OUT6

    classDef clu fill:#fef9c3,stroke:#ca8a04,color:#713f12
    class Cluster clu
```

### 4.7 Rescue 阶段 (三支路级联)

```mermaid
flowchart TB
    accTitle: Rescue 三支路级联拯救
    accDescr: A路 CheckV评估 → C路 VSI延伸 → D路 BLASTN+VSI 三级联拯救

    subgraph Rescue["run_rescue() — 三支路级联拯救"]
        direction TB
        HOSTFILTER["宿主过滤<br/>--host-filter Plant<br/>CD-HIT known + CheckV pass → 免拯救"]
        
        subgraph BranchA["分支 A: CheckV"]
            A1["checkv completeness<br/>分块并行评估 centroids"]
            A2["completeness ≥ 90%<br/>→ pass"]
            A3["< 90% → 进入分支 C"]
            A1 --> A2
            A1 --> A3
        end

        subgraph BranchC["分支 C: Virseqimprover"]
            C1["cluster 多样本 reads 聚合<br/>Salmon 定量 → BBMap 提取"]
            C2["SPAdes 组装 → CheckV"]
            C3["pass → 输出"]
            C4["fail → 进入分支 D"]
            C1 --> C2 --> C3
            C2 --> C4
        end

        subgraph BranchD["分支 D: BLASTN + VSI"]
            D1["BLASTN megablast<br/>--blast-db ref.fasta"]
            D2["CheckV + VSI"]
            D3["pass → 输出"]
            D1 --> D2 --> D3
        end

        MERGE["合并 A+C+D pass<br/>vclust 最终去重<br/>→ HQ vOTU"]
        
        HOSTFILTER --> BranchA
        A3 --> BranchC
        C4 --> BranchD
        A2 --> MERGE
        C3 --> MERGE
        D3 --> MERGE
    end

    OUT7[/"08_Rescue/{host}/<br/>centroids/final_centroids.fasta ★"/]
    Rescue --> OUT7

    classDef rescue fill:#fee2e2,stroke:#dc2626,color:#7f1d1d
    class Rescue rescue
```

---

## 五、auto_known_virus.py 全架构总览

```mermaid
flowchart TB
    accTitle: auto_known_virus.py 已知病毒三阶段分析架构
    accDescr: 从清洁 reads 出发, 依次执行快速检测/变异分析/全长组装三阶段

    subgraph Input2["输入"]
        RD["清洁 reads<br/>--reads_dir $OUT/00b_HostDepletion/"]
        REF["参考数据库<br/>--ref_info ref_info.tsv<br/>--reference ref.fasta"]
    end

    subgraph Orchestrator2["auto_known_virus.py 编排器"]
        CTRL2["自身控制<br/>--stage all/detect/variants/full<br/>--dry-run --force<br/>--log-level INFO"]
    end

    subgraph Stage1["1. detect — 快速检测"]
        direction TB
        D1["batch_virus_depth.py<br/>双引擎: Salmon/Kallisto vs Bowtie2/BWA<br/>Poisson_Ratio 去假阳性<br/>双轨过滤: A轨(全基因组) + B轨(基因区)"]
        D1_OUT["输出: 1_FastViromeExplorer/<br/>summary/summary.tsv<br/>summary/best.summary.tsv"]
        D1 --> D1_OUT
    end

    subgraph Stage2["2. variants — 变异分析"]
        direction TB
        D2["batch_virus_variants.py<br/>三池解耦并行:<br/>  [1] 提取reads → BAM<br/>  [2] 共识序列 → FreeBayes/iVar<br/>  [3] 变异检出 → VCF<br/>  [4] SnpEff 注释<br/>  [5] SnpGenie dN/dS"]
        D2_OUT["输出: 2_Virus_variants_Results/<br/>summary/all_summary.tsv<br/>virus-variants/ virus-SnpEff/<br/>virus-SNPGenie/ virus-consensus/"]
        D2 --> D2_OUT
    end

    subgraph Stage3["3. full — 全长组装"]
        direction TB
        D3["batch_virus_full.py<br/>多工具组装: SPAdes/IVA/...<br/>按 {Taxonomy}_{Accession}/<br/>   {Sample}_{Accession}/ 归档<br/>--extra_args '--iter 3 --vc-min-depth 1'"]
        D3_OUT["输出: 3_Virus_assemblies_final/<br/>{tax}_{acc}/{sample}_{acc}/<br/>final.fasta / scaffolds.fasta"]
        D3 --> D3_OUT
    end

    RD --> Orchestrator2
    REF --> Orchestrator2
    Orchestrator2 --> Stage1 --> Stage2 --> Stage3

    classDef known_input fill:#e0f2fe,stroke:#0284c7,color:#0c4a6e
    classDef known_ctrl fill:#dcfce7,stroke:#16a34a,color:#14532d
    classDef known_s1 fill:#fef3c7,stroke:#d97706,color:#78350f
    classDef known_s2 fill:#fce7f3,stroke:#db2777,color:#831843
    classDef known_s3 fill:#ede9fe,stroke:#7c3aed,color:#4c1d95

    class RD,REF known_input
    class Orchestrator2 known_ctrl
    class Stage1,D1_OUT known_s1
    class Stage2,D2_OUT known_s2
    class Stage3,D3_OUT known_s3
```

---

## 六、auto_known_virus.py 参数传递关系

```mermaid
flowchart TB
    accTitle: auto_known_virus.py 参数传递
    accDescr: 编排器收集三阶段专用参数并分发到 batch_virus_depth / batch_virus_variants / batch_virus_full

    subgraph AutoOrch["auto_known_virus.py"]
        AP["argparse<br/>★★★ 必需: reads_dir output_dir ref_info reference<br/>★★☆ 可选: 40+ 参数覆盖三阶段"]
        AM["main()<br/>--dry-run → 仅显示概览<br/>--force → 覆盖 --resume<br/>--stage → 选择阶段"]
    end

    subgraph DetectParams["Step1: batch_virus_depth.py"]
        direction LR
        DP1["--tool bowtie2/salmon/...<br/>--coverage 10.0 --ratio 0.3<br/>--sp_thresh 95.0<br/>--genes_cov (双轨B轨)<br/>--use_coverm (CoverM清洗)<br/>--single_end --keep_tmp<br/>--verbose"]
    end

    subgraph VariantsParams["Step2: batch_virus_variants.py"]
        direction LR
        VP1["--variant_caller freebayes<br/>--snpeff --snpgenie<br/>--snpeff_jar/--snpeff_config<br/>--vc_qual 20 --vc_depth 5<br/>--vc_freq 0.5 --vc_ambig N<br/>--bam (已有BAM)<br/>--no_extract_reads<br/>--no_consensus<br/>--no_call_variants<br/>--disable_dynamic_vcf"]
    end

    subgraph FullParams["Step3: batch_virus_full.py"]
        direction LR
        FP1["--assembly_tools all<br/>--min_covered 10.0<br/>--extra_args '--iter 3'<br/>--virus_full_script<br/>--gb (GenBank)"]
    end

    subgraph Common["公共参数"]
        CP2["--threads -t (默认 40)<br/>--jobs -j (默认 4)<br/>--align_threads (默认 8)<br/>--resume --force"]
    end

    AP --> AM
    AM --> DetectParams
    AM --> VariantsParams
    AM --> FullParams
    Common --> AP

    classDef orch fill:#dbeafe,stroke:#2563eb,color:#1e3a5f
    classDef params fill:#fef3c7,stroke:#d97706,color:#78350f
    classDef common fill:#dcfce7,stroke:#16a34a,color:#14532d

    class AutoOrch orch
    class DetectParams,VariantsParams,FullParams params
    class Common common
```

---

## 七、两脚本对比总览

```mermaid
mindmap
    accTitle: MMPV-RNA 双流水线对比
    accDescr: virome_pipeline.py 10阶段新病毒发现 vs auto_known_virus.py 3阶段已知病毒分析

    MMPV-RNA_v2.3
        virome_pipeline.py[新病毒发现]
            预处理
                Clean::icon(fa fa-broom)
                Deplete::icon(fa fa-filter)
            核心
                Assembly::icon(fa fa-cogs)
                Identification::icon(fa fa-search)
                COBRA::icon(fa fa-arrow-right)
                Cluster::icon(fa fa-project-diagram)
            后处理
                Taxonomy::icon(fa fa-tags)
                Host::icon(fa fa-leaf)
                CheckV::icon(fa fa-check-circle)
                Rescue::icon(fa fa-life-ring)
            产出
                新种/新属/新科
                HQ_vOTU
        auto_known_virus.py[已知病毒分析]
            检测
                detect::icon(fa fa-magnifying-glass)
                双引擎+Poisson打假
            变异
                variants::icon(fa fa-dna)
                SnpEff+SnpGenie
            组装
                full::icon(fa fa-puzzle-piece)
                多工具全长组装
            产出
                定量丰度表
                变异谱/dNdS
                全长基因组
```

---

## 八、Rescue 三支路详细序列

```mermaid
sequenceDiagram
    accTitle: Rescue 三支路执行序列
    accDescr: 宿主过滤后依次执行分支A CheckV、分支C VSI、分支D BLASTN+VSI, 最后合并去重

    participant Orch as virome_pipeline<br/>run_rescue()
    participant A as 分支A<br/>CheckV
    participant C as 分支C<br/>Virseqimprover
    participant D as 分支D<br/>BLASTN+VSI
    participant Merge as vclust合并

    Note over Orch: 加载宿主预测 + centroids

    Orch->>Orch: 按 Final_Host 分离<br/>目标宿主 / Unknown / 其他
    Orch->>Orch: 标记免拯救<br/>CD-HIT known + CheckV pass(≥90%)
    Orch->>Orch: 写入 known/centroids

    Note over Orch,A: 分支 A: CheckV 并行评估
    Orch->>A: target_novel centroids
    A->>A: 分块并行 checkv completeness
    A->>A: completeness ≥ 90% → pass
    A->>A: < 90% → fail → 分支 C

    Note over Orch,C: 分支 C: Virseqimprover reads延伸
    Orch->>C: branchA_fail centroids
    C->>C: cluster多样本reads聚合
    C->>C: Salmon定量 → BBMap提取
    C->>C: SPAdes组装 → CheckV
    C->>C: pass → 输出
    C->>C: fail → 分支 D

    Note over Orch,D: 分支 D: BLASTN 最后拯救
    Orch->>D: branchB_fail centroids
    D->>D: BLASTN megablast ref.fasta
    D->>D: CheckV → VSI
    D->>D: pass → 输出

    Note over Orch,Merge: 合并 + 最终去重
    A-->>Merge: branchA_pass
    C-->>Merge: branchB_pass
    D-->>Merge: branchC_pass
    Merge->>Merge: vclust prefilter→align→cluster
    Merge-->>Orch: final_centroids.fasta ★
    Orch->>Orch: CheckV 质量报告 汇总
```

---

*文档生成时间: 2026-06-14 | MMPV-RNA v2.3*
