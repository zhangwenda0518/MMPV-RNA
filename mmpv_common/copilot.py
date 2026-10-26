#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
copilot.py — MMPV 命令规划副驾驶 (四层递进架构 · L4 编排层)
=============================================================

自然语言 → 一条「经过白名单校验」的 MMPV 编排器命令草稿。

三层防线 (缺一不可):
  1. 脚本白名单: 只认 COMMAND_SPECS 里登记的编排器, 其余一律拒绝;
  2. 参数白名单: flag 必须在该编排器允许集合内; --stage 词元必须在
     该编排器真实 choices 里; 带子命令的编排器 (⑦ known/discovery)
     还必须给出登记过的子命令 (全部取自源码 argparse, 见各 spec 注释);
  3. shell 元字符零容忍: ; & | ` $() > < 换行 → 直接拒绝, 不解析。

红线: 本模块【不执行】任何命令 —— 只产出校验过的计划 + dry-run 预览,
执行权永远在人的终端里 (doc/AI_DRIVEN_DESIGN.md §红线)。危险 flag
(--force 等) 一律置 requires_confirmation=True 并显著提示。

两条规划路径, 结果都过同一套 validate():
  * suggest_intents()   纯规则关键词匹配, 离线可用;
  * plan_with_ai()      LLM 意图映射 → 候选 JSON → 白名单校验。
"""

from __future__ import annotations

import json
import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# 仓库根 (copilot.py 在 mmpv_common/ 下)。仅用于「脚本是否还在仓库里」的存在性判定,
# 不参与任何路径拼接执行 —— 本模块【不执行】命令。
_REPO_ROOT = Path(__file__).resolve().parent.parent

try:
    from mmpv_common import stage_cards
    from mmpv_common.ai import chat, log_ai_event, parse_json, render_prompt
except ImportError:
    from mmpv_common import stage_cards  # type: ignore
    from mmpv_common.ai import chat, log_ai_event, parse_json, render_prompt  # type: ignore

# ═══════════════════════════════════════════════════════════════════
# 命令规格 (2026-10 按各编排器真实 --help / argparse 逐条重建:
# 幻影 flag 清零、必填项与 choices 对齐; 见 tests 里的 spec×真实 CLI 回归)
# ═══════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class SubCommand:
    """顶层子命令 (仅 ⑦ submission 用): 各自的必填与 flag 面。

    `submission_pipeline.py {known,discovery}` 是 argparse subparsers (都不带
    --stage): validate 只认登记过的子命令, plan_with_ai 也只按登记生成草稿。
    """
    name: str
    desc: str = ""
    required: Tuple[str, ...] = ()
    any_required: Tuple[Tuple[str, ...], ...] = ()   # 每组内全要, 组间满足其一
    optional: Tuple[str, ...] = ()
    aliases: Tuple[Tuple[str, str], ...] = ()        # 短选项 → 长选项 (canonical 名)


@dataclass(frozen=True)
class CommandSpec:
    key: str                    # 对应 stage_cards 管线 key
    script: str                 # 仓库相对脚本路径
    desc: str
    stages: Tuple[str, ...]     # --stage 的真实 choices (该编排器无 --stage 时为空)
    required: Tuple[str, ...] = ()
    any_required: Tuple[Tuple[str, ...], ...] = ()   # 每组内全要, 组间满足其一
    optional: Tuple[str, ...] = ()
    dangerous: Tuple[str, ...] = ()
    cond_required: Tuple[Tuple[str, Tuple[str, ...]], ...] = ()  # (stage集合, flags)
    hint: str = ""
    stage_required: bool = False                     # --stage 本身是否 argparse 必填 (required=True)
    subcommands: Tuple[SubCommand, ...] = ()         # 顶层子命令 (仅 submission)
    aliases: Tuple[Tuple[str, str], ...] = ()        # 短选项 → 长选项 (canonical 名)
    # 编排器已被归档时, 填其**真实**仓库相对路径 (如 "archive/2026.../x.py")。
    # 非空 ⟹ 该 spec 只用于「历史命令仍可校验」, 不可运行:
    #   · `script` 保留历史路径 (validate/_match_spec 要认得老命令 & 老文档里的写法);
    #   · `archived_to` 指向真实位置, 两条路径都算命中;
    #   · validate() 做存在性判定, 命中归档 spec → 明确警告 + 强制人工确认门。
    archived_to: str = ""


# ── ② data_preprocessing.py (--help 逐条核对, 2026-10) ─────────────────
_PREP_FLAGS = ("--input_reads", "--output_dir", "-o", "--io-layout", "--host_db", "--kraken2_db",
               "--host_align_db", "--host_fasta", "--host_species", "--host_taxid",
               "--aligner", "--seq_type", "--rrna", "--rrna_tool", "--silva_index",
               "--kraken2_confidence", "--skip_clumpify", "--dedup", "--clumpify_memory",
               "--no_compress", "--threads", "-t", "--jobs", "-j", "--force",
               "--keep_rrna", "--tmp_dir", "--warn_retained")

# ── ③ virome_pipeline.py (156 个 flag 全量登记, 按阶段分组) ────────────
_DISC_IO_FLAGS = ("--input_reads", "--io-layout", "--cluster_input", "--input_assembly",
                  "--config", "--profile", "--dump-config")
# filter 阶段的 graph 子步骤 (默认关闭, --skip-graph 显式跳过)
_DISC_GRAPH_FLAGS = ("--skip-graph", "--graph-enabled", "--graph-min-link-pairs", "--graph-gc-diff",
                     "--graph-depth-ratio", "--graph-require-proper", "--graph-resolution",
                     "--graph-seed-min-votes", "--graph-nonviral-min-len",
                     "--graph-majority-ratio", "--graph-max-reads", "--graph-map-min-len",
                     "--graph-jobs", "--graph-keep-bam", "--graph-no-candidates",
                     "--no-graph-demote", "--no-graph-rescue")
# filter / 复核 (SDT + CT3-HMM)
_DISC_FILTER_FLAGS = ("--host-filter", "--skip-sdt", "--sdt-min-ratio", "--sdt-threads",
                      "--sdt-max-seqs", "--tree-min-ratio", "--keep-min-len", "--ct3-hmm-dir",
                      "--skip-ct3-hmm", "--no-hmm-rescue", "--no-hmm-in-score",
                      "--hmm-min-specific", "--hmm-promote-keep", "--no-hmm-promote-keep",
                      "--hmm-promote-min", "--hmm-require-family", "--hmm-count-ambiguous")
# 流程控制 + clean/deplete 细化
_DISC_RUN_FLAGS = ("--skip_clean", "--skip_depletion", "--skip_clumpify", "--dry-run",
                   "--log-level", "--stop-on-error", "--dedup", "--clumpify_memory",
                   "--no_compress", "--clean_debug", "--deplete_tmp", "--kraken2_confidence",
                   "--keep_rrna", "--keep_taxids", "--vmtouch", "--rrna_chunk_size",
                   "--rrna_report", "--deplete_steps", "--align_config", "--deplete_debug",
                   "--refinec_threads", "--refinec_frag_min_len", "--refinec_min_id",
                   "--refinec_min_cov", "--n-total-samples")
# assembly / identification / cobra / cluster
_DISC_ASM_FLAGS = ("--asm_tmp_dir", "--asm_keep_temp", "--nr_db", "--skip_uniprot_filter",
                   "--skip_nr_filter", "--skip_id_plots", "--clean_failed", "--ident_ext",
                   "--filter-mode", "--no-merge", "--cobra_mink", "--cobra_maxk",
                   "--cobra_jobs", "--cobra_linkage_mismatch", "--cobra_verbose",
                   "--skip_vclust", "--vclust_cluster_file", "--skip-rmdup", "--rmdup-length",
                   "--skip-flye", "--flye-min-overlap", "--flye-read-error")
# taxonomy / host
_DISC_TAX_FLAGS = ("--tax_tools", "--tax_jobs", "--tax_ext", "--tax_sample_name",
                   "--tax_remove_suffix", "--skip_rnavirhost", "--skip_phabox", "--skip_ictv")
# 数据库与资源 (含工具选择/阈值/线程)
_DISC_DB_FLAGS = ("--host_db", "--kraken2_db", "--host_align_db", "--virus_db", "--checkv_db",
                  "--blast-db", "--aligner", "--seq_type", "--rrna", "--rrna_tool",
                  "--silva_index", "--coassembly", "--assembler", "--contig-length", "-l",
                  "--identify_tools", "--virus_mode", "--genomad_db", "--cat_db", "--cat_tax",
                  "--mmseqs_db", "--vitap_db", "--acvirus_db", "--vcontact3_db", "--viroids_db",
                  "--virsorter_db", "--viralverify_hmm", "--metabuli_db", "--virus_protein_db",
                  "--uniprot_db", "--virus_taxid", "--virhunter_path", "--virhunter_weights",
                  "--virbot_path", "--viralm_path", "--virsorter_group", "--blast_mode",
                  "--blast_evalue", "--blast_top_n", "--phabox-db", "--host-mode", "--prob-dir",
                  "--ref-genomes", "--threads", "-t", "--memory", "-m", "--jobs", "-j",
                  "--min-length", "--ani", "--qcov", "--cdhit_ani", "--cdhit_qcov",
                  "--virseqimprover-path", "--salmon-bin", "--max_vsi_samples", "--min_vsi_len",
                  "--checkv_threshold")

# ── ④ auto_known_virus.py (87 个 flag 全量登记) ───────────────────────
_ANA_IO_FLAGS = ("--reads_dir", "--output_dir", "-o", "--io-layout", "--ref_info",
                 "--reference", "--no-resume", "--force", "--dry_run", "--log_level",
                 "--threads", "--jobs", "--align_threads", "--profile", "--config")
_ANA_DETECT_FLAGS = ("--tool", "--batch_size", "--coverage", "--ratio", "--meandepth",
                     "--min_tpm", "--min_uniq_reads", "--sp_thresh", "--genes_cov",
                     "--min_gene_total_cov", "--min_gene_avr_cov", "--taxid_clusters",
                     "--use_coverm", "--min_aln_len", "--min_aln_prop", "--min_pid",
                     "--single_end", "--keep_tmp", "--verbose")
_ANA_FILTER_FLAGS = ("--filter", "--no-filter", "--filter_cov", "--filter_depth",
                     "--filter_reads", "--filter_keyword", "--filter_tpm", "--filter_poisson")
_ANA_VARIANT_FLAGS = ("--variant_caller", "--snpeff", "--snpeff_jar", "--snpeff_config",
                      "--snpeff_mem", "--snpgenie", "--no_call_variants", "--bam",
                      "--disable_dynamic_vcf", "--vc_qual", "--vc_depth", "--vc_freq",
                      "--vc_ambig")
_ANA_FULL_FLAGS = ("--no_extract_reads", "--no_consensus", "--assembly_tools", "--min_covered",
                   "--extra_args", "--virus_full_script", "--gb", "--extract_target",
                   "--max_n_genome", "--min_length")
_ANA_POST_FLAGS = ("--post_min_dp", "--post_min_af", "--skip_vcf_viz", "--skip_vcf_merge",
                   "--meta", "--metadata_tsv", "--skip_metadata_assoc", "--skip_snpeff_macro",
                   "--skip_maftools", "--skip_snpgenie")
_ANA_SIM_DVG_FLAGS = ("--sim_ref", "--sim_mode", "--sim_cdhit", "--virema_script",
                      "--dvg_seed", "--dvg_mindel", "--dvg_min_cov", "--dvg_shm",
                      "--dvg_reads")
_ANA_REPORT_FLAGS = ("--ai_api_key", "--ai_model")

# ── ① public_data_pipeline.py (40 个 flag 全量登记) ───────────────────
_PUB_FLAGS = ("--deepseek-api", "--deepseek-model", "--ncbi-api", "--work-dir", "--io-layout",
              "--bin-dir", "--source-type", "--detailed", "--no-detailed", "--prescreen-filter",
              "--prescreen-jobs", "--prescreen-gsa-jobs", "--prescreen-virus-pct-ident",
              "--prescreen-virus-pct", "--prescreen-no-gsa", "--virus-profile",
              "--no-virus-profile", "--virus-profile-plot", "--no-virus-profile-plot",
              "--virus-profile-min-samples", "--virus-profile-top-n", "--virus-profile-by",
              "--taxdump", "--ngdc-method", "--ngdc-concurrency", "--prefetch-concurrency",
              "--skip-list", "--threads", "--convert-jobs", "--convert-threads",
              "--host-species", "--host-taxid", "--host-fasta", "--host-download-source",
              "--hostdb-tools", "--force", "--dry-run")

# ── ① build_host_pipeline.py (15 个 flag 全量登记) ────────────────────
_HOSTREF_BUILD_FLAGS = ("--ncbi-api", "--work-dir", "--bin-dir", "--genome-fasta",
                        "--download-source", "--taxonomy-dir", "--hostdb-tools",
                        "--seq-types", "--k2-libs", "--threads", "--dry-run")

# ── ⑤ phylo_pipeline.py: 只保留真实存在的公共 flag (-j/--jobs、--resume、
#    --log_level 该编排器没有, 已从旧 _COMMON 里剔除) ────────────────
_PHYLO_COMMON_FLAGS = ("-t", "--threads", "--io-layout", "--config", "--force")

# ── ⑦ submission_pipeline.py: 顶层是 subparsers (known/discovery), 无 --stage ──
_SUB_COMMON_FLAGS = ("--suvtk-db", "--run-title", "-t", "--threads")

COMMAND_SPECS: Dict[str, CommandSpec] = {
    "preprocessing": CommandSpec(
        key="preprocessing",
        script="data_preprocessing_pipeline/data_preprocessing.py",
        desc="② 质控+去宿主 (clean/deplete/bbnorm/report)",
        stages=("all", "clean", "deplete", "bbnorm", "report"),
        required=("--output_dir",),
        any_required=(("--input_reads",), ("--output_dir",)),  # report/deplete 可缺 reads
        optional=_PREP_FLAGS,
        aliases=(("-o", "--output_dir"), ("-t", "--threads"), ("-j", "--jobs")),
        dangerous=("--force",),
        cond_required=((("clean", "all"), ("--input_reads",)),),
        hint="python data_preprocessing_pipeline/data_preprocessing.py --stage clean deplete "
             "--input_reads raw/ --output_dir <PROJECT_ROOT>/02_Preprocessing --host_db <HOST_DB> -t 40 -j 10",
    ),
    "discovery": CommandSpec(
        key="discovery",
        script="virome_discovery_pipeline/virome_pipeline.py",
        desc="③ 病毒发现 15 阶段 (clean→…→report)",
        stages=("all", "clean", "deplete", "assembly", "identification", "filter",
                "cobra", "merge", "cluster", "taxonomy", "host", "checkv",
                "rescue", "analysis", "analysis_verify", "report"),
        required=("--output_dir",),
        optional=(_DISC_IO_FLAGS + _DISC_GRAPH_FLAGS + _DISC_FILTER_FLAGS
                  + _DISC_RUN_FLAGS + _DISC_ASM_FLAGS + _DISC_TAX_FLAGS + _DISC_DB_FLAGS),
        aliases=(("-t", "--threads"), ("-j", "--jobs"), ("-l", "--contig-length"),
                 ("-m", "--memory")),
        dangerous=("--force",),
        cond_required=((("clean", "deplete", "assembly", "all"), ("--input_reads",)),),
        hint="python virome_discovery_pipeline/virome_pipeline.py --stage all "
             "--output_dir <PROJECT_ROOT>/03_Discovery --input_reads <DEPLETED_READS_DIR> "
             "--host_db <HOST_DB> --virus_db <VIRUS_DB> --checkv_db <CHECKV_DB> -t 120 -j 20",
    ),
    "analysis": CommandSpec(
        key="analysis",
        script="virome_analysis_pipeline/auto_known_virus.py",
        desc="④ 已知病毒深度分析 (detect/variants/full/dvg/report…)",
        stages=("all", "detect", "filter", "variants", "full", "extract", "post",
                "similarity", "dvg", "report"),
        # 该编排器手工校验 reads_dir/output_dir/ref_info/reference 四件套
        # (--config 可从 profile 里补, 见 auto_known_virus.py 参数校验段)
        any_required=(("--reads_dir", "--reference", "--ref_info", "--output_dir"),
                      ("--config",)),
        optional=(_ANA_IO_FLAGS + _ANA_DETECT_FLAGS + _ANA_FILTER_FLAGS
                  + _ANA_VARIANT_FLAGS + _ANA_FULL_FLAGS + _ANA_POST_FLAGS
                  + _ANA_SIM_DVG_FLAGS + _ANA_REPORT_FLAGS),
        aliases=(("-o", "--output_dir"),),
        dangerous=("--force",),
        hint="python virome_analysis_pipeline/auto_known_virus.py --stage all "
             "--reads_dir <DEPLETED_READS_DIR> --reference <ref.fasta> "
             "--ref_info <ref_info.tsv> --output_dir <PROJECT_ROOT>/04_Analysis "
             "--threads 40 --jobs 4",
    ),
    "public": CommandSpec(
        key="public",
        script="public_metadata_pipeline/public_data_pipeline.py",
        desc="① 公共数据检索/元数据/下载/可视化",
        stages=("search", "info", "prescreen", "down", "convert", "plot",
                "hostref", "hostdb", "report", "all"),
        required=("--species", "--taxid"),
        stage_required=True,          # --stage 是 argparse 必填 (见 --help usage)
        optional=_PUB_FLAGS,
        dangerous=("--deepseek-api", "--ncbi-api", "--force"),  # 内联 key → 强提醒走环境变量
        hint="python public_metadata_pipeline/public_data_pipeline.py "
             "--species 'Solanum lycopersicum' --taxid 4081 --stage search info down plot",
    ),
    # build_host_pipeline.py 独立编排器 (dict 键 "hostref"; 别名见下方 hostref_build)。
    # ① 公共管道 public_data_pipeline.py 自身的 hostref/hostdb 两个 --stage 值
    #   归上面 "public" spec (同脚本不重复登记, 否则 _match_spec 歧义)。
    "hostref": CommandSpec(
        key="hostref",
        script="public_metadata_pipeline/build_host_pipeline.py",
        desc="① 宿主参考基因组下载+四索引建库",
        stages=("genome-down", "hostdb", "all"),
        required=("--species", "--taxid"),
        stage_required=True,          # --stage 是 argparse 必填 (见 --help usage)
        optional=_HOSTREF_BUILD_FLAGS,
        dangerous=("--ncbi-api", "--force"),
        hint="python public_metadata_pipeline/build_host_pipeline.py "
             "--species 'Solanum lycopersicum' --taxid 4081 --stage all --threads 30",
    ),
    # ⑥ EVE 已归档 2026-10-11 → archive/20261011_endogenous_virus_pipeline/。
    # spec 保留 (test_agent_tools/test_ai_driven 的 validate 用例依赖; 历史命令仍可校验)。
    # `script` 故意保留归档前的历史路径 —— 老日志/老文档里的命令写法必须继续认得;
    # 真实位置写在 `archived_to`, _match_spec 两条路径都命中, validate 会做存在性判定
    # 并给出「已归档, 真实路径在此」的警告 + 强制人工确认门 (永不自动执行)。
    "eve": CommandSpec(
        key="eve",
        script="endogenous_virus_pipeline/eve_screen.py",
        archived_to="archive/20261011_endogenous_virus_pipeline/eve_screen.py",
        desc="⑥ EVE 内源病毒筛查【已归档, 不可运行】(discover/verdict/annotate)",
        stages=("all", "discover", "verdict", "annotate"),
        required=("-o",),
        any_required=(("-g",), ("-B",)),
        optional=("-t", "-J", "-n", "-o", "-g", "-B", "--fast", "--merge", "--no-merge",
                  "--overlap", "--cmd-timeout", "--io-layout", "--resume", "--force"),
        dangerous=("--force",),
        hint="python endogenous_virus_pipeline/eve_screen.py -g <HOST_GENOME_FASTA> "
             "-n <SPECIES_NAME> -o <PROJECT_ROOT>/06_EVE -t 40",
    ),
    "phylo": CommandSpec(
        key="phylo",
        script="virome_phylo_pipeline/phylo_pipeline.py",
        desc="⑤ 系统发育/进化 18 stage 8 模块",
        stages=("all", "prep", "metadata", "online", "clean", "align", "splitstree",
                "rdp5", "tree", "popgen", "host", "capheine", "clock", "beast",
                "gene_dating", "phylogeo", "geo", "geo_paths", "report",
                # 8 大模块分组名 (--stage 顶层名)
                "data", "phylogeny", "recomb", "select", "time", "geography"),
        required=(),
        any_required=(("--config",), ("--extract_dir", "--metadata")),
        optional=_PHYLO_COMMON_FLAGS + ("--extract_dir", "--metadata", "--meta"),
        aliases=(("-t", "--threads"),),
        dangerous=("--force",),
        hint="python virome_phylo_pipeline/phylo_pipeline.py "
             "--extract_dir <④06_Extraction> --metadata <metadata.tsv> --stage all -t 20",
    ),
    "submission": CommandSpec(
        key="submission",
        script="virome_submission_pipeline/submission_pipeline.py",
        desc="⑦ GenBank/CNCB 提交 (.sqn) — 子命令 known / discovery",
        stages=(),                      # 该编排器无 --stage, 顶层是子命令
        optional=_SUB_COMMON_FLAGS,
        aliases=(("-t", "--threads"),),
        dangerous=(),
        hint="python virome_submission_pipeline/submission_pipeline.py discovery "
             "--work-dir <WORK_DIR> --run-title <RUN_TITLE> --mode both -t 40",
        subcommands=(
            SubCommand(
                name="known", desc="已知病毒提交 (INFO → suvtk)",
                # 序列来源二选一: all_summary 表 或 ④06_Extraction 全长产物
                any_required=(("--all-summary",), ("--extraction-dir",)),
                optional=("--all-summary", "--extraction-dir", "--fasta-dir", "--metadata-dir",
                          "-o", "--output", "--assembler", "--sequencer", "--enrichment",
                          "--pipeline", "--source-uvig", "--metagenome-source", "--authors"),
                aliases=(("-o", "--output"),),
            ),
            SubCommand(
                name="discovery", desc="新病毒提交 (discovery pipeline 的拯救/Rescue 产物)",
                required=("--work-dir",),
                optional=("--work-dir", "--fasta", "--rvdb-db", "--mode", "--fetch-metadata",
                          "--metadata-dir", "--config", "--interactive", "--pipeline-type",
                          "--ref-info", "--force-taxonomy", "--force-features", "--force-hypo",
                          "--export-metadata", "--skip-gate"),
            ),
        ),
    ),
}
# 命名别名: "hostref" 与 "hostref_build" 指向同一 spec (build_host_pipeline.py)。
# 两个键共享一个 CommandSpec 对象 —— _match_spec 按脚本匹配结果一致, 不会歧义。
COMMAND_SPECS["hostref_build"] = COMMAND_SPECS["hostref"]

SHELL_METACHARS = re.compile(r"[;|&`><\n]|\$\(")


# ═══════════════════════════════════════════════════════════════════
# 解析与校验
# ═══════════════════════════════════════════════════════════════════

@dataclass
class Plan:
    spec: CommandSpec
    tokens: List[str]
    stages: List[str] = field(default_factory=list)
    flags: Dict[str, str] = field(default_factory=dict)   # flag → value("" for store_true)
    subcommand: str = ""                                  # 命中的顶层子命令 (无则空)


@dataclass
class ValidationResult:
    ok: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    requires_confirmation: bool = False
    preview: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {"ok": self.ok, "errors": self.errors, "warnings": self.warnings,
                "requires_confirmation": self.requires_confirmation,
                "preview": self.preview}


def _match_spec(script_token: str) -> Optional[CommandSpec]:
    t = script_token.replace("\\", "/").strip()
    for spec in COMMAND_SPECS.values():
        if t.endswith(spec.script) or t == spec.script:
            return spec
        # 已归档 spec: 归档后的真实路径同样命中 (两条路径都认, 老命令不失效)
        if spec.archived_to and (t.endswith(spec.archived_to) or t == spec.archived_to):
            return spec
    return None


def _script_state(spec: CommandSpec) -> str:
    """脚本落盘状态: "live" | "archived" | "missing"。

    仅做存在性判定, 不执行任何东西。用 `archived_to` 区分「本来就没写」和
    「曾经有、现在被归档了」—— 后者要给出可操作的真实路径, 不能只说文件不存在。
    """
    if (_REPO_ROOT / spec.script).is_file():
        return "live"
    if spec.archived_to and (_REPO_ROOT / spec.archived_to).is_file():
        return "archived"
    return "missing"


def _archived_warning(spec: CommandSpec) -> str:
    """归档 spec 的警告文案 —— 必须给出真实路径, 否则等于换个地方断链。"""
    return (f"脚本已归档: {spec.script} 不在仓库现状里; 真实路径 "
            f"{spec.archived_to} (归档 ≠ 可运行: 该 spec 仅用于历史命令校验, "
            f"不建议直接执行归档副本)")


# 允许出现在脚本之前的解释器/启动器词元 (python / pixi run 等)
_INTERPRETER_TOKENS = {"python", "python3", "pixi", "run", "uv", "uvicorn"}
# 解释器 basename 白名单: python / python3 / python3.11 / python.exe … 或 py (Windows 启动器)。
# 只认"python[数字[.数字]][.exe]" —— 收紧了旧判据 base.startswith("python")
# (pythonXxx / pythonevil 之类 basename 会被当解释器剥掉, 配合 agent 直执行是真风险)。
_PY_LAUNCHER_RE = re.compile(r"python\d*(\.\d+)*(\.exe)?|py(\.exe)?")


def _split_interpreter(tokens: List[str]) -> Tuple[List[str], Optional[str]]:
    """剥掉开头的解释器词元, 返回 (剩余词元, 命中的脚本词元或 None)。"""
    i = 0
    while i < len(tokens):
        bare = tokens[i].strip().strip('"\'')
        base = bare.replace("\\", "/").rsplit("/", 1)[-1].lower()
        if base in _INTERPRETER_TOKENS or _PY_LAUNCHER_RE.fullmatch(base):
            i += 1
            continue
        return tokens[i:], bare
    return [], None


def parse_tokens(command: str) -> List[str]:
    """切词 (posix=False 保住 Windows 反斜杠), 去包裹引号; 元字符先行拒绝。"""
    if SHELL_METACHARS.search(command):
        raise ValueError("命令含 shell 元字符 (; | & ` > < $( ) 换行), 已拒绝")
    toks = shlex.split(command, posix=False)
    return [t[1:-1] if len(t) >= 2 and t[0] == t[-1] and t[0] in "\"'" else t
            for t in toks]


def validate(command: str) -> ValidationResult:
    """白名单校验。永远不执行 —— 只裁决。"""
    try:
        tokens = parse_tokens(command)
    except ValueError as e:
        return ValidationResult(ok=False, errors=[str(e)])
    if not tokens:
        return ValidationResult(ok=False, errors=["空命令"])
    rest, script_tok = _split_interpreter(tokens)
    spec = _match_spec(script_tok or "")
    if spec is None:
        return ValidationResult(
            ok=False,
            errors=[f"脚本不在白名单: {script_tok or '(空)'} (允许: "
                    f"{', '.join(dict.fromkeys(s.script for s in COMMAND_SPECS.values()))})"])

    errors: List[str] = []
    warnings: List[str] = []
    dangerous_hit: List[str] = []
    stages: List[str] = []
    flags: Dict[str, str] = {}
    accepts_stage = bool(spec.stages)     # 该编排器是否有 --stage (choices 非空即有)

    full_tokens = list(tokens)   # Plan 保留完整命令 (含解释器+脚本) 用于 dry-run 预览
    body = rest[1:]              # rest[0] 是已匹配白名单的脚本词元

    # ── 子命令 (⑦ submission {known,discovery}): 首个非 flag 词元必须是其一 ──
    sub: Optional[SubCommand] = None
    if spec.subcommands:
        names = [s.name for s in spec.subcommands]
        if body and body[0] in names:
            sub = next(s for s in spec.subcommands if s.name == body[0])
            body = body[1:]
        elif body and not body[0].startswith("-"):
            errors.append(f"子命令不在 {spec.key} 允许集合: {body[0]} (允许: {', '.join(names)})")
            body = body[1:]
        else:
            errors.append(f"缺子命令 (允许: {', '.join(names)})")

    # flag 全量面: 通用面 + 各子命令面 (子命令未给对时报真实错, 不级联报一串未知 flag)
    known = set(spec.required) | set(spec.optional) | set(spec.dangerous)
    for s in spec.subcommands:
        known |= set(s.required) | set(s.optional)
    if accepts_stage:
        known.add("--stage")
    # 短选项 → 长选项: -o/-t/-j/-l/-m 归一到 canonical 名, 必填判定与取值才对得上
    aliases = dict(spec.aliases)
    if sub is not None:
        aliases.update(sub.aliases)

    i = 0
    while i < len(body):
        tok = body[i]
        if tok in ("-h", "--help"):
            flags[tok] = ""
            i += 1
            continue
        if tok == "--stage" and accepts_stage:
            # --stage 后的词元全是 stage 名 (nargs='+'), 直到下一个 - 开头 flag
            i += 1
            while i < len(body) and not body[i].startswith("-"):
                s = body[i]
                if s not in spec.stages:
                    errors.append(f"--stage 词元不在 {spec.key} choices 里: {s}")
                else:
                    stages.append(s)
                i += 1
            continue
        if tok.startswith("-"):
            flag = aliases.get(tok, tok)
            if flag not in known:
                errors.append(f"未知 flag: {tok} ({spec.key} 允许: {', '.join(sorted(known))})")
                i += 1
                continue
            if flag in spec.dangerous:
                dangerous_hit.append(flag)
            # 取值: 下一词元若不以 - 开头则视为该 flag 的值
            if i + 1 < len(body) and not body[i + 1].startswith("-"):
                flags[flag] = body[i + 1]
                i += 2
            else:
                flags[flag] = ""
                i += 1
            continue
        errors.append(f"多余位置参数: {tok}")
        i += 1

    # ── 必填: 通用 required/any_required + 命中的子命令各自的 required/any_required ──
    required = list(spec.required) + (list(sub.required) if sub is not None else [])
    for req in required:
        if req not in flags:
            errors.append(f"缺必填参数 {req}")
    groups = list(spec.any_required) + (list(sub.any_required) if sub is not None else [])
    if groups and not any(all(f in flags for f in g) for g in groups):
        errors.append("缺必填参数组: 满足其一即可 — " +
                      " 或 ".join("+".join(g) for g in groups))
    if spec.stage_required and not stages:
        errors.append("缺必填参数 --stage (该编排器 --stage 为 argparse 必填)")
    stage_set = set(stages)
    for stage_names, req_flags in spec.cond_required:
        if stage_set & set(stage_names):
            for rf in req_flags:
                if rf not in flags:
                    errors.append(f"stage 含 {stage_names[0]} 时缺 {rf}")
    if accepts_stage and not stages and not spec.stage_required:
        warnings.append("未指定 --stage, 编排器将按其默认 (通常为 all) 执行")
    if "--deepseek-api" in flags or "--ncbi-api" in flags:
        warnings.append("检测到内联 API key 参数: 建议改用环境变量 "
                        "(DEEPSEEK_API_KEY / NCBI_API_KEY), 避免密钥进 shell 历史/日志")
    for f_ in dangerous_hit:
        warnings.append(f"危险参数 {f_}: 将强制重跑已完成阶段, 请确认磁盘与时长预算")

    # ── 归档断链防护: spec 的脚本已不在仓库现状 (被归档) 时, 说清楚真实路径 ──
    # 这类 spec 保留只为「历史命令仍可校验」, 绝不能被当成可跑命令放行:
    # ok 时强制人工确认门 (agent 收到 requires_confirmation 且无 confirm_fn 会直接跳过)。
    archived = _script_state(spec) == "archived"
    if archived:
        warnings.append(_archived_warning(spec))

    plan = Plan(spec=spec, tokens=full_tokens, stages=stages, flags=flags,
                subcommand=sub.name if sub is not None else "")
    preview = dry_run_preview(plan)
    # 危险 flag 才需要"先确认再跑"; 校验失败就是拒绝 (ok=False), 不是"确认后可跑"。
    # 归档 spec 校验通过时同样置确认门 (它压根不该被执行)。
    return ValidationResult(ok=not errors, errors=errors, warnings=warnings,
                            requires_confirmation=bool(dangerous_hit) or (archived and not errors),
                            preview=preview)


def dry_run_preview(plan: Plan) -> str:
    lines = [
        f"[dry-run 预览] 管线: {plan.spec.key} — {plan.spec.desc}",
        f"脚本: {plan.spec.script}",
    ]
    if _script_state(plan.spec) == "archived":
        # 预览里也要显式标出, 否则人只看到「脚本: <旧路径>」会以为能直接粘终端跑。
        # 用 [归档] 而非 emoji: Windows GBK 控制台 print 会 UnicodeEncodeError。
        lines.append(f"[归档] 不可运行; 真实路径: {plan.spec.archived_to}")
    if plan.subcommand:
        lines.append(f"子命令: {plan.subcommand}")
    lines.append(f"stage: {' '.join(plan.stages) if plan.stages else '(编排器默认)'}")
    for k, v in plan.flags.items():
        lines.append(f"  {k}" + (f" {v}" if v else ""))
    lines.append("完整命令:")
    lines.append("  " + " ".join(plan.tokens))
    lines.append("(本副驾驶不执行命令; 请复制到终端人工执行)")
    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════
# 意图规划: 规则路径 (离线)
# ═══════════════════════════════════════════════════════════════════

_KEYWORDS: Tuple[Tuple[Tuple[str, ...], str], ...] = (
    (("sra", "gsa", "检索", "下载", "公共数据", "下载测序"), "public"),
    (("质控", "清洗", "去宿主", "clean", "deplete", "预处理"), "preprocessing"),
    (("发现", "denovo", "de novo", "组装", "鉴定", "聚类", "votu", "拯救", "rescue"), "discovery"),
    (("定量", "变异", "snp", "全长", "dvg", "已知病毒"), "analysis"),
    (("进化", "发育", "系统", "树", "地理", "定年", "选择压力", "phylo"), "phylo"),
    (("eve", "内源", "化石"), "eve"),
    (("提交", "genbank", "sqn", "sequin", "ncbi 提交"), "submission"),
    (("宿主参考", "建库", "索引", "hostdb", "kraken2 库"), "hostref"),
)


def suggest_intents(text: str, top_k: int = 3) -> List[Dict[str, str]]:
    """关键词 → 管线建议 + 示例命令。纯规则, 离线。"""
    low = text.lower()
    scored: List[Tuple[int, str]] = []
    for kws, key in _KEYWORDS:
        score = sum(1 for k in kws if k in low)
        if score:
            scored.append((score, key))
    scored.sort(reverse=True)
    out = []
    for _, key in scored[:top_k]:
        spec = COMMAND_SPECS.get(key)
        if spec:
            out.append({"pipeline": key, "desc": spec.desc, "example": spec.hint})
    return out


# ═══════════════════════════════════════════════════════════════════
# 意图规划: AI 路径 (输出仍过 validate())
# ═══════════════════════════════════════════════════════════════════

def _stages_text() -> str:
    """各管线合法 stage / 子命令清单 (喂 AI 意图提示词)。

    "hostref"/"hostref_build" 是同脚本别名 → 只列一次, 免得提示词里出现重复行。
    """
    lines = []
    seen_scripts = set()
    for spec in COMMAND_SPECS.values():
        if spec.script in seen_scripts:
            continue
        seen_scripts.add(spec.script)
        if spec.subcommands:
            choices = "子命令: " + " ".join(s.name for s in spec.subcommands)
        else:
            choices = " ".join(spec.stages)
        tag = " [已归档, 不可运行]" if spec.archived_to else ""
        lines.append(f"{spec.key} ({spec.script}){tag}: {choices}")
    return "\n".join(lines)


def plan_with_ai(text: str, *,
                 state_dir=None, config: Any = None,
                 client: Any = None) -> Dict[str, Any]:
    """NL → AI 候选命令 → 白名单校验。返回 {status, plan|reason, validation}。

    AI 的输出从不被直接执行: 只作为草稿进 validate(), 校验不过就是不过。
    """
    pipelines_text = stage_cards.summary_text(max_stage_desc_chars=40)
    prompt_text, meta = render_prompt("chat_intent",
                                      pipelines_text=pipelines_text,
                                      stages_text=_stages_text())
    res = chat([
        {"role": "system", "content": prompt_text},
        {"role": "user", "content": text},
    ], config=config, client=client, json_mode=True, temperature=0.0)
    if res.status != "ok":
        return {"status": res.status, "reason": res.reason,
                "fallback": suggest_intents(text)}
    try:
        intent = parse_json(res.text)
    except ValueError as e:
        return {"status": "error", "reason": f"AI 意图不可解析: {e}",
                "fallback": suggest_intents(text)}
    if not isinstance(intent, dict) or intent.get("pipeline") not in COMMAND_SPECS:
        return {"status": "error", "reason": "AI 给出的管线不在白名单",
                "fallback": suggest_intents(text)}
    spec = COMMAND_SPECS[intent["pipeline"]]
    tokens = ["python", spec.script]
    # 带子命令的编排器 (⑦ submission): AI 选对才拼; 选错/漏给交给 validate 报错
    sub = str(intent.get("subcommand") or "")
    if spec.subcommands and sub in [s.name for s in spec.subcommands]:
        tokens.append(sub)
    st = intent.get("stages") or []
    if st:
        tokens += ["--stage"] + [str(s) for s in st]
    known = set(spec.required) | set(spec.optional) | set(spec.dangerous)
    for s_ in spec.subcommands:
        known |= set(s_.required) | set(s_.optional)
    accepted_flags: List[str] = []
    for flag, val in (intent.get("args") or {}).items():
        if flag in known:
            accepted_flags.append(flag)
            tokens += [flag] + ([str(val)] if str(val) else [])
    validation = validate(" ".join(tokens))
    if state_dir is not None:
        # 只记结构信息 (管线名 + 子命令 + 通过白名单的 flag 键名), 不落用户自由文本/
        # AI 原始输出 —— ai_log.jsonl 是论文溯源凭证, 不该沉淀用户原话与模型胡写
        log_ai_event(state_dir, {"kind": "call", "op": "plan_with_ai",
                                 "model": res.model, "prompt": meta.to_dict(),
                                 "pipeline": intent.get("pipeline"),
                                 "subcommand": sub if spec.subcommands else "",
                                 "flags": sorted(accepted_flags),
                                 "ok": validation.ok})
    return {"status": "ok", "intent": intent,
            "missing": intent.get("missing", []),
            "explain": intent.get("explain", ""),
            "command": " ".join(tokens),
            "validation": validation.to_dict()}
