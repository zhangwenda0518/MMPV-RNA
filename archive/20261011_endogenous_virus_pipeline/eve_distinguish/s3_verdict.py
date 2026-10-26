#!/usr/bin/env python3
"""S3 判别 [判定逻辑最后变更: 模块 v7]: 单一位点架构 (S2/S2b) x 密码子退化 (S1) -> verdict.

v3 的规则把"完整"当 EVE 证据: ORF 完整 + loci>=2 或 span>=1500 -> review_extend;
stops>=3 -> ancient_EVE. 两个方向都不对:
  * RNA-seq 里真 DNA 病毒本来就拼不全, "不全" 不能推出 EVE;
  * stops>=3 是绝对计数, 被组装 indel 污染 (v4 S1 已改用同批读框基线 + 分布剖面).
v4 判定表 (只对 Caulimoviridae / 面板可用的 DNA 候选):

  密码子退化 \\ 单一位点架构   齐全(>=4/5 组件, 或 S2b 位点齐全)   部分/单组件
  -------------------------  ----------------------------------  ---------------
  分布式退化 (化石签名)       EVE_STRONG_provirus (整合前病毒)     ancient_EVE / EVE_suspect
  完好 (编码受保留)           virus_candidate (感染中的真病毒)     virus_fragment_review
  组装断裂 (单点 indel)       assembly_breakpoint_review          assembly_breakpoint_review
  组成本底噪声                compositional_noise_review         compositional_noise_review

  注意"部分 + 退化"只给 EVE_suspect: 它既可能是化石碎片, 也可能是有 indel 的真病毒碎片,
  contig 级无法定案, 必须靠 S2b 位点或寄主基因组.
  非 Cauli (面板不适用): 只走结构通道, 退而不判 EVE.

v5 修订 (2026-09-23):
  * 未寄主映射候选 (s2b 的 locus_arch=no_host_locus) **不参与**"单位点齐全"判定。
    旧版它们被 s2b 合并成伪位点, 若干碎片凑出"齐全"并伪造 EVE_STRONG_provirus。
  * 新增 host_scope 列 (own_species / cross_species_guess / none): 同一条否决里
    "自身物种"和"跨物种旁证"是两种强度的证据, 分开写, 下游可按此列过滤。
  * provirus_scale 列写 TRUE/FALSE (旧版写 Python 的 True/False, 按 "TRUE" 过滤的
    下游脚本会静默漏掉这些行)。

v6 修订 (2026-09-23): **寄主否决不再能一票定生死**。旧规则
"自身寄主 >=90%/80% 且无 MP/CP/AP -> host_contamination_likely" 在 1027 条真实候选上
被调用 48 次, 每一次都落在带上游病毒科注释的 contig 上 (Rhabdoviridae 12 / Geminiviridae
12 / Betaflexiviridae 6 / Aspiviridae 3 / Virgaviridae 3 / ...), 没有一条是寄主基因 ——
有一条 885 nt 的候选 nt_pident 93.5% 对着 Betaflexiviridae 参考, 仍被判成"寄主污染"。
三个缺陷叠在一起: 跨物种命中也走同一条否决; `not has_mp_cp_ap` 对非 Cauli 候选是空转的
(面板 222 条全是 Caulimoviridae, 压根评估不了 Geminiviridae 的结构基因, "没测到"不等于
"没有"); 上游病毒注释与寄主同源冲突时旧版无条件采信后者。现在改成三道闸:

  * host_scope=cross_species_guess -> host_homology_cross_species (REVIEW, 不动数据):
    样本在映射表里没有条目, 命中来自"哪个物种比得上就拿哪个"的启发式, 不当自身寄主证据;
  * own_species 但上游给了病毒科, 或蛋白/核苷酸通道有对已知序列的命中
    -> host_conflict_review (REVIEW): 两个通道给出互相矛盾的答案, contig 级定不了案,
    不由 s3 仲裁;
  * 只有 own_species 且三条同源通道全空, 才给 host_contamination_likely (REMOVE)。

另两处 v6 改动: `locus_full` 不再把 s2 的 provirus_scale 算进来 —— 旧版那样写, 19 条
EVE_STRONG_provirus 里有 6 条 arch 是 no_host_locus (根本没有寄主位点), 旗标的
ERR2040732_NODE_11 就在其中; verdict 表的 locus_ncomp 改取 s2b 的位点级值, 不再用 s2
的单 contig 值 (旧版因此出现 arch=locus_full 而 ncomp=2 的自相矛盾行)。

v6.1 修订 (2026-09-23): 补上核苷酸通道。上面那道"上游有病毒信号"原先只读了病毒科和
蛋白一致两条, 于是 48 条旧 REMOVE 里有 2 条只带 blastn 命中 (71.1% / 72.1% 对着
PP728250.1 / HQ633072.1) 的候选被漏判成删除。现在 nt_pident 非空且带物种名同样算一条
病毒同源记录, 不设百分比下限 —— 这是直接序列同一性, 上游 blastn 既已写下物种名, 记录
本身就是证据; 低于多少算噪声该由上游定, 不由 s3 例代。

**"上游有病毒信号"只认这三条通道**, 不把 CDD 域级低一致命中算进去: 实测剩下的候选里,
CDD 命中的是 dUTPase / RdRp / RNase H 这类广谱域, 寄主基因同样带, 驳不倒"这是寄主序列"。
同理 aa_pident 30–50% 的命中也不算。所以 REMOVE 分支查的不是"上游无任何病毒信号"这么
绝对的一句话, 而是"三条同源通道都没有记录", 判词里也照这个写。

输入: s1_decay.tsv s2_domains.tsv [s2b_locus_architecture.tsv] rescue_evidence_scored.tsv <out.tsv>
      [refgenome_evidence.tsv]

v6.3 修订: **参考基因组通道 (S2c, 可选第 6 参)**。候选比对回自身物种的染色体级参考
基因组后, 有两类整合信号是 contig 级通道 (S1/S2/S2b) 拿不到的:
  * EVE_flank_confirmed —— contig 两端各一段高一致宿主比对, 且落在同一条染色体同链、
    间隔有限 = 基因组片段跨着整合位点被组装出来 (前病毒铁证);
  * EVE_splice_chimera —— 宿主外显子嵌合/剪接样式 (contig 级代理, read 级待验证).
调整是**单向**的: 只把 review 档往 EVE 方向升、或把 REMOVE/KEEP 与整合证据的矛盾
转成 host_conflict_review 人工复核; 不比对上 (no_hit / no_ref_genome / 无结构) **不
改任何 verdict** —— 参考基因组缺失是常态, "没比对上"绝不是"真病毒"的证据。
host_dominant (参考基因组聚合覆盖占优且无 MP/CP/AP) 走与 S2b 同款的三道闸:
上游有病毒信号 -> host_conflict_review, 否则 host_contamination_likely。
rg_call / rg_scope 两列无论调整与否都随表输出, 下游可按列过滤。

v7 修订 (P0-2, 2026-10-11 审查实测): **evidence 覆盖率闸门**。旧版 `ev.get(cid, {})`
对缺行候选静默取空 —— fam 落空串 -> upstream_viral=False, 若该候选又恰有自身物种
寄主同源, 就落进 host_contamination_likely, 即判别链里唯一会删数据的分支; 实测同一
候选只换一份缺行的 evidence 表, verdict 就从 host_conflict_review(REVIEW) 翻成
REMOVE。v7 起缺行候选一律判 **missing_evidence** (REVIEW 类, s4 映射为 REVIEW),
三条上游通道没评估过的候选不得自动 REMOVE/MOVE_EVE/KEEP; s2c 调整对它同样不生效。
开跑时打印覆盖率, 并把 query/evidence 的路径 + sha256 + 行数 + 覆盖清单写进
`<out>.provenance.tsv` 溯源。run_all.sh 另有前置 ID 重叠预检 (零重叠直接中止)。

v7.1 修订 (2026-10-17, 审稿意见落实):
  * **no_panel 不再允许 REMOVE**: 面板只含 Caulimoviridae, 对非 Cauli 候选零命中是
    "面板不适用"而非"确认无结构基因"; 旧版只加注脚仍判 host_contamination_likely
    (唯一删数据分支)。现在该情形判 **host_homology_panel_na** (REVIEW), 并写明
    缺的是哪类证据。
  * **evidence_tier 列**: A_genomic_flank (s2c 侧翼/嵌合) > B_locus_architecture
    (s2b 位点齐全) > C_decay_only (仅退化通道) > D_none。退化通道与整合位点证据
    强度不同, 不许合并表述; s4 --strict-eve 用它把 C/D 档 MOVE_EVE 降为 REVIEW。
  * 组件判定改为**整词集合** (旧版子串匹配对拼接串的语义依赖拼写)。

v7.5 修订 (2026-10-20, 文献对齐 P0): **coding_intact + locus_full 不再单独放行 KEEP**。
vead088 实测 32.3% 的基因组 EVE ORF 完整, Rozo-Lopez 2023 的 5 个全长组成性表达转录本
全是 EVE —— "编码完好 + 位点齐全"是外源的必要条件而非充分条件。第 8 参 (可选) 给
eve_ref_classify 的双库分类表 (ref_classify.tsv) 时, virus_candidate 分支先复查双库:
ref_class 为 endogenous_supported (EVE 元件库 bits 领先) 或 both_families (两侧并立)
的候选转 **virus_eve_conflict_review** (REVIEW 类, s4 人工定案), 不自动 KEEP。不给
该表时行为与 v7.4 完全一致 (向后兼容)。
"""
import sys
import os
from collections import Counter

HOST_PIDENT, HOST_COV = 90.0, 0.80   # 与自身寄主物种 >=90% 一致且区间覆盖 >=80% -> 寄主序列
CANON = ("MP", "CP", "AP", "RT", "RH")   # s2 的规范组件集 (判定按整词集合, 不按子串)
POL_COMPS = {"RT", "RH"}   # 纯 pol 区 (无 MP/CP/AP) 更像转座子/退化 pol
AA_VIRAL_PID = 95.0        # 上游最佳蛋白命中一致到这个程度 -> "这是某个已知蛋白的序列"
UPSTREAM_NA = {"", "-", "NA", "N/A", "未知", "unknown", "None"}   # 上游没给出科的写法

# v6.3: S2c 的整合证据可以把这些 review 档 verdict 往 EVE 方向升 (单向表)。
# 不在表里的 (EVE_STRONG_provirus / virus_candidate / host_*) 各有专属处理或不动。
RG_EVE_UPGRADE_FROM = {"review", "virus_fragment_review", "compositional_noise_review",
                       "structure_intact_review", "structure_decay_review",
                       "assembly_breakpoint_review", "EVE_suspect", "ancient_EVE",
                       "EVE_LTR_TE", "host_conflict_review",
                       "host_homology_cross_species"}


def rg_adjust(verdict, rg_call, upstream_viral, upstream_hit):
    """S2c 参考基因组证据 -> verdict 的单向调整。返回调整后的 verdict。

    冲突一律转 host_conflict_review (两通道矛盾, 人工定案, 不由本函数仲裁):
      * 整合证据 (侧翼/嵌合) 撞上 virus_candidate (判成感染中的真病毒) 或
        host_contamination_likely (判成寄主序列) —— 两种都是矛盾;
      * host_dominant 撞上上游病毒信号 (科注释 / 高一致蛋白或核苷酸命中)。
    """
    if rg_call == "host_dominant":
        if upstream_viral or upstream_hit:
            return "host_conflict_review"
        return "host_contamination_likely"
    if rg_call in ("EVE_flank_confirmed", "EVE_splice_chimera"):
        if verdict in ("virus_candidate", "host_contamination_likely"):
            return "host_conflict_review"
        if verdict in RG_EVE_UPGRADE_FROM:
            return "EVE_STRONG_provirus" if rg_call == "EVE_flank_confirmed" \
                else "EVE_suspect"
    return verdict

# 输出表头: 单独抽成常量, 测试直接 from s3_verdict import VERDICT_HDR 引用,
# 免得两头各抄一份、加列时只有一边改 (旧版测过 s2b 表头就是这么漏的)
VERDICT_HDR = ("contig_id\ttax_family\tcategory\tcheckv_completeness\tlocus_ncomp\tlocus_completeness\t"
               "provirus_scale\tlocus_arch\tlocus_members\tdecay_class\tstop_enrichment\t"
               "premature_stops_region\tframe_switches\torf_max_fraction\tte_flag\thost_wpid\t"
               "host_cov\thost_scope\trg_call\trg_scope\tevidence_tier\t"
               "recon_stop\trecon_frameshift\trecon_status\tverdict\tverdict_reason\n")
VERDICT_COLS = VERDICT_HDR.rstrip("\n").split("\t")

# v7.3: 蛋白重建层的交叉验证 (可选第 7 参, eve_reconstruct.py 产物)
#   recon_stop / recon_frameshift 是 GeneWise 在重建比对串上核验过的终止/移码数;
#   recon_status: ok / filtered / miniprot_only / no_hit / tool_missing
# 用途: s1 的 decay_class 是在**原始序列**上数的终止, 会被组装 indel 污染;
# 重建层的计数是独立通道。两者矛盾时降级为 REVIEW, 而不是各自坚持。
RECON_TIER_UPGRADE = {"ok"}

# 证据层级 (v7.1, 审稿意见落实): 不同强度的证据不许合并表述。
#   A_genomic_flank       基因组侧翼/嵌合 (s2c): 真正的整合位点证据
#   B_locus_architecture  宿主位点级组件齐全 (s2b locus_full)
#   C_decay_only          仅密码子退化通道 (无位点/侧翼证据)
#   D_none                以上都没有 (证据不足)
# 下游按此列过滤/分层报告; s4 --strict-eve 用它把 C/D 档的 MOVE_EVE 降为 REVIEW。
TIER_FLANK = "A_genomic_flank"
TIER_LOCUS = "B_locus_architecture"
TIER_DECAY = "C_decay_only"
TIER_NONE = "D_none"


def evidence_tier(rg_call, locus_full, decay):
    """该候选可用的**最强**证据层级 (与 verdict 无关, 只描述证据形态)."""
    if rg_call in ("EVE_flank_confirmed", "EVE_splice_chimera"):
        return TIER_FLANK
    if locus_full:
        return TIER_LOCUS
    if decay == "distributed_decay":
        return TIER_DECAY
    return TIER_NONE


def read_tsv(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        d = {}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p:
                d[p[0]] = dict(zip(hdr, p))
        return d


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _sha256(path):
    """文件内容指纹 (溯源用); 文件不存在返回空串 (read_tsv 同样容错)."""
    if not path or not os.path.exists(path):
        return ""
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_provenance(out_tsv, s1_tsv, ev_tsv, n_ev_rows, ids, missing_ev, extra_ev):
    """P0-2 验收: 记录 query/evidence 的来源标识与内容指纹.

    写 <out>.provenance.tsv (field\\tvalue): 双方路径 + sha256 + 行数 + 覆盖清单。
    误传旧批次/不同候选集时, 复核者拿这份文件能直接对出"这份 verdict 是拿哪两份
    输入算的", 而不是只看到一列安静的 verdict。
    """
    prov = f"{out_tsv}.provenance.tsv"
    with open(prov, "w", encoding="utf-8", newline="\n") as fo:
        fo.write("field\tvalue\n")

        def clip(lst):
            s = ";".join(lst[:50])
            return s + (" ...(+%d)" % (len(lst) - 50) if len(lst) > 50 else "")

        for k, v in (
            ("query_tsv", s1_tsv),
            ("query_sha256", _sha256(s1_tsv)),
            ("query_rows", len(ids)),
            ("evidence_tsv", ev_tsv),
            ("evidence_sha256", _sha256(ev_tsv)),
            ("evidence_rows", n_ev_rows),
            ("covered_candidates", len(ids) - len(missing_ev)),
            ("missing_evidence_ids", clip(missing_ev)),
            ("extra_evidence_ids", clip(extra_ev)),
        ):
            fo.write(f"{k}\t{v}\n")
    return prov


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if len(args) < 5:
        print(__doc__)
        sys.exit(1)
    s1_tsv, s2_tsv, s2b_tsv, ev_tsv, out_tsv = args[0], args[1], args[2], args[3], args[4]
    rg_tsv = args[5] if len(args) > 5 else ""   # v6.3: 可选, s2c 的 refgenome_evidence.tsv
    recon_tsv = args[6] if len(args) > 6 else ""  # v7.3: 可选, eve_reconstruct 的 recon.tsv
    refclass_tsv = args[7] if len(args) > 7 else ""  # v7.5: 可选, eve_ref_classify 的 ref_classify.tsv
    s2b = read_tsv(s2b_tsv)
    s1 = read_tsv(s1_tsv)
    s2 = read_tsv(s2_tsv)

    ev = read_tsv(ev_tsv)
    rg = read_tsv(rg_tsv)
    recon = read_tsv(recon_tsv)
    rc = read_tsv(refclass_tsv)   # v7.5: 双库分类 (P0 复查用), 缺省空表 = 行为不变

    # ---- v7 P0-2 闸门: evidence 表覆盖率 (逐行闸门在循环里, 这里只报账) ----
    ids = list(s1)
    missing_ev = [cid for cid in ids if cid not in ev]
    extra_ev = sorted(cid for cid in ev if cid not in s1)
    print(f"[gate] evidence 覆盖 {len(ids) - len(missing_ev)}/{len(ids)} 候选"
          + (f"; 缺行 {len(missing_ev)} 条 -> missing_evidence (REVIEW), "
             "不进入 REMOVE/MOVE_EVE/KEEP" if missing_ev else "")
          + (f"; evidence 多出 {len(extra_ev)} 行不在候选集里 (疑似批次错配): "
             f"{extra_ev[:3]}" if extra_ev else ""))
    if missing_ev:
        print("[gate] WARN: evidence 缺行 = 上游三条通道从未评估这些候选, "
              "s3 拒绝为它们下发任何可执行动作", file=sys.stderr)
    prov = write_provenance(out_tsv, s1_tsv, ev_tsv, len(ev), ids,
                            missing_ev, extra_ev)
    print(f"[gate] 溯源指纹 -> {prov}")

    out = open(out_tsv, "w")
    out.write(VERDICT_HDR)
    verdicts = Counter()
    tiers_seen = []
    for cid, d1 in s1.items():
        d2 = s2.get(cid, {})
        db = s2b.get(cid, {})
        e = ev.get(cid)
        has_ev = e is not None
        if not has_ev:
            e = {}
        fam = (e.get("tax_family", "") or "").strip()
        checkv = _f(e.get("checkv_completeness"))
        # v6: locus_ncomp 取 s2b 的位点级值。旧版取 s2 的 best_locus_ncomp(单 contig
        # 级), 于是出现 arch=locus_full 而 ncomp=2 的自相矛盾行; s2b 缺失时才退回 s2。
        _ln = _f(db.get("locus_ncomp"))
        best_n = int(_ln) if _ln is not None else int(_f(d2.get("best_locus_ncomp")) or 0)
        comps_str = d2.get("components", "-") or "-"
        provirus = d2.get("provirus_scale", "FALSE") == "TRUE"
        locus_arch = db.get("locus_arch", "")
        locus_hint = db.get("locus_hint", "")
        # v6: 只认真位点。旧版 `or provirus` 让 s2 的 provirus_scale 顶掉 s2b 的
        # no_host_locus —— 那 6 条 EVE_STRONG 根本没有寄主位点, 是这么来的。
        locus_full = locus_arch == "locus_full"
        decay = d1.get("decay_class", "") or "no_hsp"
        enrich = _f(d1.get("stop_enrichment"))
        stops = _f(d1.get("premature_stops_region"))
        stops = int(stops) if stops is not None else None
        switches = _f(d1.get("frame_switches"))
        switches = int(switches) if switches is not None else None
        te = d2.get("te_flag", "FALSE")
        te_bits = _f(d2.get("te_bitscore")) or 0
        cauli_bits = _f(d2.get("cauli_bitscore")) or 0
        # v7.1: 按整词集合判定, 不再子串匹配 —— `"MP" in comps_str` 对 "-" 之外
        # 的任何拼接串都成立与否取决于拼写, 而组件名是两位/三位缩写; 集合判定与
        # s2 的 CANON 序输出同一口径 (s2 历史 bug 正是按字符迭代把 AP 撕成 A+P)。
        comp_set = {c for c in comps_str.split("+") if c in CANON}
        has_mp_cp_ap = bool(comp_set & {"MP", "CP", "AP"})
        no_panel = not comp_set           # 面板零命中: MP/CP/AP 缺席是否不可评估
        pol_only = bool(comp_set) and not has_mp_cp_ap
        hpid, hcov = _f(db.get("host_wpid")) or 0.0, _f(db.get("host_cov")) or 0.0
        host_scope = db.get("host_scope", "") or ("none" if not db else "")
        hcode = db.get("host_code", "") or "-"
        aa_pid = _f(e.get("aa_pident")) or 0.0
        aa_spec = (e.get("aa_species", "") or "").strip()
        nt_pid = _f(e.get("nt_pident")) or 0.0
        nt_spec = (e.get("nt_species", "") or "").strip()
        # 上游已经把这条候选归到某个病毒科 / 或高一致命中某个已知序列 —— 无论哪种,
        # "它是寄主序列"这个结论都不能只由 s2b 的同源强度单方面下达
        upstream_viral = fam not in UPSTREAM_NA
        upstream_aa = aa_pid >= AA_VIRAL_PID and bool(aa_spec)
        # 核苷酸通道是蛋白通道的对等口径: 上游 blastn 拿得到一致度和物种名, 就是一条
        # 病毒同源记录。v6.1 补上它 —— 48 条旧 REMOVE 里有 2 条只带 nt 命中 (71.1% /
        # 72.1% 对着 PP728250.1 / HQ633072.1), 只读科和蛋白通道就把它们漏判成删除了。
        # 不设百分比下限: 这是直接序列同一性, 上游 blastn 既然写下了物种名, 记录本身
        # 就说明它在病毒库里比中过东西。低一致的自有阈值只该由上游定, 不该由 s3 例代。
        upstream_nt = nt_pid > 0.0 and bool(nt_spec)
        upstream_hit = upstream_aa or upstream_nt
        reason = []
        # v6.3: emit 不再立即写行 —— verdict 要等参考基因组通道 (可选) 单向调整后
        # 才落定, reason 也可能被 rg 块续写 (如 virus_candidate 分支的 locus_hint)。
        # 写行挪到循环体末尾, verdict 计数也改在定稿后统计。
        verdict = ""

        def emit(v):
            nonlocal verdict
            verdict = v

        # ---- v7 P0-2 闸门: evidence 缺行的候选不得进入任何自动判定分支 ----
        # 三条上游通道 (病毒科 / 蛋白 / 核苷酸) 对它从未评估, "没评估"绝不是
        # "评估过且无病毒信号"; 删数据的 host_contamination_likely 与 MOVE_EVE /
        # KEEP_virus 都必须有 evidence 在场才允许下达。
        if not has_ev:
            reason.append("evidence 表缺该候选行: 上游三条同源通道(病毒科/蛋白/核苷酸)"
                          "从未评估, 不得自动 REMOVE/MOVE_EVE/KEEP, 转人工复核 (P0-2 闸门)")
            emit("missing_evidence")

        # ---- 寄主同源否决: 候选就是寄主的序列 (如 HSP70 假阳性) ----
        # OneKP 转录组组装是碎的, 所以用物种级区间并集覆盖而非单条比对覆盖.
        # 但"同源强度够"只是必要条件, 不是充分条件: 一条广宿主范围的植物病毒在混合
        # 转录组里同样能对自己的寄主打到 95-100% 一致。v6 起按证据口径分流, 见 docstring。
        elif hpid >= HOST_PIDENT and hcov >= HOST_COV and not has_mp_cp_ap:
            if host_scope != "own_species":
                reason.append(f"跨物种旁证 (样本无自身物种映射, 命中物种 {hcode}): 聚合 "
                              f"{hpid:.1f}% 一致性 / {hcov:.2f} 区间覆盖; 口径不是自身物种, "
                              "不能据此判寄主序列")
                emit("host_homology_cross_species")
            elif upstream_viral or upstream_hit:
                why = []
                if upstream_viral:
                    why.append(f"上游归为 {fam}")
                if upstream_aa:
                    why.append(f"最佳蛋白命中 {aa_pid:.1f}% 一致到 {aa_spec}")
                if upstream_nt:
                    why.append(f"核苷酸命中 {nt_pid:.1f}% 一致到 {nt_spec}")
                reason.append(f"自身寄主物种 {hcode} 聚合 {hpid:.1f}% 一致性 / {hcov:.2f} "
                              "区间覆盖, 但同时存在病毒侧证据(" + "; ".join(why) + "): "
                              "两通道结论矛盾, 需人工定案, 不据此删数据")
                emit("host_conflict_review")
            else:
                # v7.1 (审稿意见): no_panel = 面板对这条候选**不适用**, "无 MP/CP/AP"
                # 是"没测到"而不是"没有"。旧版只加一句注脚, verdict 仍落到唯一会删
                # 数据的 host_contamination_likely —— 用只含 Caulimoviridae 的面板
                # 给非 Cauli 候选下"确证寄主序列"的结论, 逻辑上不成立。
                # 现在拆成两个 verdict: 面板适用时才允许 REMOVE; 不适用时只到 REVIEW,
                # 并明确写出"缺的是结构基因阴性证据, 需补该科的蛋白面板或基因组比对"。
                if no_panel:
                    reason.append(
                        f"自身寄主物种 {hcode} 聚合 {hpid:.1f}% 一致性 / {hcov:.2f} "
                        "区间覆盖; 但本候选在 Caulimoviridae 专用面板上零命中, "
                        "MP/CP/AP 阴性无法评估 (面板不适用), 上游三条同源通道也无记录 "
                        "-> 不足以判寄主序列, 转人工 (补该科面板或基因组比对后复判)")
                    emit("host_homology_panel_na")
                else:
                    reason.append(f"自身寄主物种 {hcode} 聚合 {hpid:.1f}% 一致性 / {hcov:.2f} "
                                  "区间覆盖, 无 MP/CP/AP 病毒结构基因, 上游三条同源通道"
                                  "(病毒科 / 蛋白打到已知序列 / 核苷酸命中)均无记录 "
                                  "-> 判寄主序列")
                    emit("host_contamination_likely")

        # ---- 非 Cauli / 无面板命中: 只走结构通道, 退而不判 EVE ----
        elif fam != "Caulimoviridae" and not comp_set:
            if te == "TRUE" and te_bits >= 100:
                reason.append("TE bitscore 领先且无病毒结构基因")
                emit("EVE_LTR_TE")
            elif decay == "coding_intact":
                reason.append("结构完整/编码保留, 面板不适用")
                emit("structure_intact_review")
            elif decay == "distributed_decay":
                reason.append("分布式退化但无 Cauli 面板证据, 不下 EVE 结论")
                emit("structure_decay_review")
            else:
                reason.append("证据不足")
                emit("review")

        # ---- TE 否决 (仅适用于无 MP/CP/AP 的纯 pol 区) ----
        elif te == "TRUE" and te_bits > cauli_bits and pol_only:
            reason.append("纯 pol 区且 TE bitscore 领先 -> 转座子")
            emit("EVE_LTR_TE")
        elif fam == "Metaviridae":
            reason.append("Metaviridae 即 LTR 逆转录转座子科")
            emit("EVE_LTR_TE")

        else:
            # ---- 二维判定: 单一位点架构 x 密码子退化 ----
            if decay == "distributed_decay":
                if locus_full:
                    reason.append("单一位点架构齐全(>=4/5 组件)且密码子分布式退化 = 整合前病毒")
                    emit("EVE_STRONG_provirus")
                elif pol_only:
                    reason.append("纯 pol 区分布式退化")
                    emit("ancient_EVE")
                else:
                    reason.append("架构不全 + 分布式退化; contig 级无法定案, 需位点/寄主基因组")
                    emit("EVE_suspect")
            elif decay == "coding_intact":
                if locus_full:
                    cv = "" if checkv is None else f", CheckV 完整度 {checkv:.0f}%"
                    # v7.5 (P0, 文献对齐): vead088 实测 32.3% 基因组 EVE ORF 完整,
                    # Rozo-Lopez 2023 的 5 个全长组成性转录本全是 EVE —— "编码完好 +
                    # 位点齐全"是外源的必要条件而非充分条件。双库分类 (可选第 8 参)
                    # 可用时必须复查: EVE 元件库 bits 领先 (endogenous_supported) 或
                    # 两侧并立 (both_families) 时, KEEP 结论与同源证据矛盾, 转 REVIEW
                    # 人工定案。表缺席 (ref_class 空) 时维持 v7.4 行为。
                    rc_call = ((rc.get(cid) or {}).get("ref_class", "") or "").strip()
                    if rc_call in ("endogenous_supported", "both_families"):
                        _rc = rc.get(cid) or {}
                        reason.append(
                            "单一位点架构齐全且密码子完好, 但双库分类 %s (EVE bits %s / "
                            "病毒 bits %s): 文献口径下编码完整不是外源证据, "
                            "KEEP 与同源证据矛盾 -> 人工复核 (v7.5 P0)"
                            % (rc_call, _rc.get("eve_bits", "?"), _rc.get("virus_bits", "?")))
                        if locus_hint == "locus_full_but_assembled_breakpoint":
                            reason.append("位点齐全但存在组装断裂, 复拼后可升级")
                        emit("virus_eve_conflict_review")
                    else:
                        reason.append("单一位点架构齐全且密码子完好 = 感染中的真病毒 (非 EVE)" + cv)
                        if locus_hint == "locus_full_but_assembled_breakpoint":
                            reason.append("位点齐全但存在组装断裂, 复拼后可升级")
                        emit("virus_candidate")
                else:
                    reason.append("架构不全 + 编码完好: RNA-seq 里真 DNA 病毒本就拼不全, 不判 EVE")
                    emit("virus_fragment_review")
            elif decay == "assembly_breakpoint":
                reason.append("终止集中在单个断裂点 (同批读框基线已排除组成本底), "
                              "更像组装 indel, 建议 rescue 延伸后复判")
                emit("assembly_breakpoint_review")
            elif decay == "compositional_noise":
                reason.append("病毒读框不比其他读框更脏, 终止来自组成/组装噪声, 无退化证据")
                emit("compositional_noise_review")
            else:
                reason.append("同源区过短或模式不明")
                emit("review")

        # ---- v6.3 参考基因组通道 (可选第 6 参, s2c 产物): 证据列 + 单向调整 ----
        # rg_call 为空 (没给 s2c / 该候选无参考基因组 / 比对零命中) 时不动任何 verdict:
        # "没比对上"不是真病毒的证据。调整只发生两种情况: review 档往 EVE 方向升,
        # 或整合证据与 REMOVE/KEEP 的既有结论矛盾时转 host_conflict_review 人工定案。
        rg_row = rg.get(cid) or {}
        rg_call = (rg_row.get("rg_call", "") or "").strip()
        rg_scope = (rg_row.get("rg_scope", "") or "").strip()
        if rg_call and verdict == "missing_evidence":
            # P0-2: evidence 缺行的候选连基础判定都没做, s2c 的单向调整更不能
            # 替它做决定 (host_dominant 若放行会把它推回 host_contamination_likely)
            reason.append("s2c 有 %s 但 evidence 缺行, 调整不生效 (先补 evidence)"
                          % rg_call)
        elif rg_call:
            reason.append("s2c %s: %s"
                          % (rg_call, (rg_row.get("call_reason", "") or "").strip()))
            new_v = rg_adjust(verdict, rg_call, upstream_viral, upstream_hit)
            if new_v != verdict:
                reason.append("s2c 调整 %s -> %s" % (verdict, new_v))
                verdict = new_v
        # ---- v7.3 蛋白重建层交叉验证 (可选第 7 参) ----
        # s1 的 decay_class 在**原始序列**上数终止, 会被组装 indel 污染;
        # 重建层 (miniprot 定位 + GeneWise 计数) 是独立通道, 每个 X 都核对过密码子。
        # 两者矛盾 -> 转 REVIEW (不各自坚持); 重建层确认有事件 -> 允许升级 tier。
        tier = evidence_tier(rg_call, locus_full, decay)
        tiers_seen.append(tier)
        rc_row = recon.get(cid) or {}
        recon_stop = (rc_row.get("n_stop", "") or "").strip()
        recon_fs = (rc_row.get("n_frameshift", "") or "").strip()
        recon_status = (rc_row.get("status", "") or "").strip()
        recon_events = None
        if recon_stop.isdigit() and recon_fs.isdigit():
            recon_events = int(recon_stop) + int(recon_fs)
        if recon_status and recon_status not in ("tool_missing", "no_hit"):
            if recon_status == "filtered":
                reason.append("重建层: 该位点的突变被质量过滤器判为不可靠 "
                              "(晚期突变/单事件低局部一致性/多事件质量不足)")
                if verdict in ("EVE_STRONG_provirus", "ancient_EVE"):
                    reason.append("重建层过滤 -> %s 降为 REVIEW" % verdict)
                    verdict = "review"
            elif recon_events == 0 and decay == "distributed_decay":
                reason.append("重建层 0 终止/0 移码, 与 s1 的 distributed_decay 矛盾 "
                              "(疑组装 indel 造成的错读框) -> 转 REVIEW")
                verdict = "review"
            elif recon_events and recon_events > 0:
                reason.append("重建层确认 %d 个事件 (stop %s / fs %s), 支持退化判定"
                              % (recon_events, recon_stop, recon_fs))
                if tier == TIER_DECAY:
                    tier = TIER_LOCUS if locus_full else tier

        verdicts[verdict] += 1
        out.write(f"{cid}\t{fam}\t{e.get('category','')}\t{e.get('checkv_completeness','')}\t"
                  f"{best_n}\t{d2.get('locus_completeness','')}\t"
                  f"{'TRUE' if provirus else 'FALSE'}\t{locus_arch}\t"
                  f"{db.get('locus_members','')}\t{decay}\t{enrich}\t{stops}\t{switches}\t"
                  f"{d1.get('orf_max_fraction','')}\t{te}\t{hpid}\t{hcov}\t{host_scope}\t"
                  f"{rg_call}\t{rg_scope}\t{tier}\t"
                  f"{recon_stop}\t{recon_fs}\t{recon_status}\t"
                  f"{verdict}\t{'|'.join(reason)}\n")
    out.close()
    print(f"[done] s3(logic v7, evidence-gated) -> {out_tsv}")
    print("verdict 分布:")
    for v, n in verdicts.most_common():
        print(f"  {v:30s} {n}")
    # v7.4: 证据层级覆盖统计 —— 让"有多少候选真有基因组侧翼证据"一眼可见。
    # 这是审稿口径的关键数字: 若 A 档为 0, 说明 s2c 没跑或没有参考基因组,
    # 全部 EVE 判定都只建立在转录本位点/退化通道上。
    tier_counts = Counter(t for t in tiers_seen)
    print("evidence_tier 分布:")
    for t in (TIER_FLANK, TIER_LOCUS, TIER_DECAY, TIER_NONE):
        print(f"  {t:26s} {tier_counts.get(t, 0)}")
    if not tier_counts.get(TIER_FLANK):
        print("[gate] WARN: A 档 (基因组侧翼证据) 为 0 —— s2c 未运行或无参考基因组。"
              "此时所有 MOVE_EVE 都只依赖转录本位点/退化通道, 不能声称'整合位点已确认'。",
              file=sys.stderr)


if __name__ == "__main__":
    main()
