#!/usr/bin/env python3
"""S3 判别 [判定逻辑最后变更: 模块 v6.3]: 单一位点架构 (S2/S2b) x 密码子退化 (S1) -> verdict.

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
"""
import sys
import os
from collections import Counter

HOST_PIDENT, HOST_COV = 90.0, 0.80   # 与自身寄主物种 >=90% 一致且区间覆盖 >=80% -> 寄主序列
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
               "host_cov\thost_scope\trg_call\trg_scope\tverdict\tverdict_reason\n")
VERDICT_COLS = VERDICT_HDR.rstrip("\n").split("\t")


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


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if len(args) < 5:
        print(__doc__)
        sys.exit(1)
    s1_tsv, s2_tsv, s2b_tsv, ev_tsv, out_tsv = args[0], args[1], args[2], args[3], args[4]
    rg_tsv = args[5] if len(args) > 5 else ""   # v6.3: 可选, s2c 的 refgenome_evidence.tsv
    s2b = read_tsv(s2b_tsv)
    s1 = read_tsv(s1_tsv)
    s2 = read_tsv(s2_tsv)

    ev = read_tsv(ev_tsv)
    rg = read_tsv(rg_tsv)

    out = open(out_tsv, "w")
    out.write(VERDICT_HDR)
    verdicts = Counter()
    for cid, d1 in s1.items():
        d2 = s2.get(cid, {})
        db = s2b.get(cid, {})
        e = ev.get(cid, {})
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
        has_mp_cp_ap = any(c in comps_str for c in ("MP", "CP", "AP"))
        no_panel = comps_str in ("-", "")   # 面板零命中: MP/CP/AP 缺席是否不可评估
        pol_only = bool(comps_str not in ("-", "")) and not has_mp_cp_ap
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

        # ---- 寄主同源否决: 候选就是寄主的序列 (如 HSP70 假阳性) ----
        # OneKP 转录组组装是碎的, 所以用物种级区间并集覆盖而非单条比对覆盖.
        # 但"同源强度够"只是必要条件, 不是充分条件: 一条广宿主范围的植物病毒在混合
        # 转录组里同样能对自己的寄主打到 95-100% 一致。v6 起按证据口径分流, 见 docstring。
        if hpid >= HOST_PIDENT and hcov >= HOST_COV and not has_mp_cp_ap:
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
                note = ""
                if no_panel:
                    # 面板 222 条全是 Caulimoviridae: 非 Cauli 候选零命中是"没法测",
                    # 不是"测了没有"。这里只在无任何上游病毒信号时才走到, 说明一下口径。
                    note = (" (注: 面板为 Cauli 专用, 本候选零命中属'面板不适用'而非"
                            "'确认无病毒结构基因')")
                reason.append(f"自身寄主物种 {hcode} 聚合 {hpid:.1f}% 一致性 / {hcov:.2f} "
                              "区间覆盖, 无 MP/CP/AP 病毒结构基因, 上游三条同源通道"
                              "(病毒科 / 蛋白打到已知序列 / 核苷酸命中)均无记录"
                              + note + " -> 判寄主序列")
                emit("host_contamination_likely")

        # ---- 非 Cauli / 无面板命中: 只走结构通道, 退而不判 EVE ----
        elif fam != "Caulimoviridae" and comps_str in ("-", ""):
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
        if rg_call:
            reason.append("s2c %s: %s"
                          % (rg_call, (rg_row.get("call_reason", "") or "").strip()))
            new_v = rg_adjust(verdict, rg_call, upstream_viral, upstream_hit)
            if new_v != verdict:
                reason.append("s2c 调整 %s -> %s" % (verdict, new_v))
                verdict = new_v
        verdicts[verdict] += 1
        out.write(f"{cid}\t{fam}\t{e.get('category','')}\t{e.get('checkv_completeness','')}\t"
                  f"{best_n}\t{d2.get('locus_completeness','')}\t"
                  f"{'TRUE' if provirus else 'FALSE'}\t{locus_arch}\t"
                  f"{db.get('locus_members','')}\t{decay}\t{enrich}\t{stops}\t{switches}\t"
                  f"{d1.get('orf_max_fraction','')}\t{te}\t{hpid}\t{hcov}\t{host_scope}\t"
                  f"{rg_call}\t{rg_scope}\t{verdict}\t{'|'.join(reason)}\n")
    out.close()
    print(f"[done] s3(logic v6.3) -> {out_tsv}")
    print("verdict 分布:")
    for v, n in verdicts.most_common():
        print(f"  {v:30s} {n}")


if __name__ == "__main__":
    main()
