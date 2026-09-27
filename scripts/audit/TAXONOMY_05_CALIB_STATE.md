# 05 分类选择层：科属不一致现状

核查时间：2026-09-15（服务器实测，非推断）
范围：**仅** `05_Taxonomy` 的分类选择与共识投票层。不涉及宿主判定。

---

## 一、05 的选型链路（实测）

```
04_CLUSTER/4_centroids/final_centroids.fasta
  └─ virus_classifier.py (md5 7cbe6915fdfc93666d4922e3a58f93c4, mtime 2026-07-25 12:07, 未动)
       7 工具并跑: genomad / metabuli / CAT / diamond_lca / VITAP / mmseqs / ACVirus
       合并 → Votus.classed/Votus_combined_taxonomy.tsv   (枸杞 62,373 条)
  └─ virus_classifier_analysis.R (md5 eee0957279a794bdc2bd2ba8ae85fa90, mtime 2026-09-14 22:31)
       逐 rank 加权投票 (RANK_DEPTH_WEIGHTS Realm=1 … Species=128; TOOL_BIAS ACVirus 1.2 … diamond_lca 0.8)
       → Votus.integrated/final_integrated_classification.tsv (枸杞 20,892 行 / md5 12e8f1715922ee967db9b10cb70e523a)
```

科属自洽的闸门**在 05 这一层**，不在下游任何一步。

## 二、机制：09-14 往 R 共识脚本里补了两道

四级备份链（可逐级回退）：

| 备份 | 时间 | 大小 | 对应改动 |
|---|---|---|---|
| `.bak` | 2026-08-04 21:26 | 37,642 B | 原始 |
| `.bak_harmonize_20260904` | 2026-09-04 02:00 | 38,827 B | 属-种协调 |
| `.bak_famgenus_20260914` | 2026-09-04 02:00 | 38,875 B | 科-属校准 |
| `.bak_rankcontain_20260914` | 2026-09-14 20:51 | 42,175 B | 逐级相容性约束 |
| 现版 | 2026-09-14 22:31 | 47,419 B | 上述两者入库 |

调用点在共识流程末端（先清属，再收尾）：

```r
442:  wide <- harmonize_family_genus(wide, stacked, tool_weights)
445:  wide <- enforce_rank_containment(wide, GENUS_FAMILY_REF)
```

- `harmonize_family_genus`（:246）：Genus 必须能追溯到「自报 Family == 共识 Family」的工具所报的属；不在候选集则换票首属，无候选则置空。
- `enforce_rank_containment`（:306）：行内 Genus 在**两把尺子**（NCBI `rankedlineage.dmp` 首次出现科 / VMR MSL41 `acvirus_db/taxa.txt`）都收录、且**两者给出的科都与行内 Family 不同** → 置空 Genus + Species。单侧判据一律不动。
- fail-safe：参照表缺失即跳过；`clean_all_ranks`（:108-109）另把非 `viridae$` 结尾的 Family、`viridae|virinae$` 结尾的 Genus 置 NA。
- 参照表：`~/database/taxonomy/genus_family_ref.tsv`，17,009,255 B，234,804 属，mtime 2026-09-14 22:29。

两道闸门自 09-14 起**下次跑 05 自动生效**；已有输出未重跑，改走后处理补丁。

## 三、结果：枸杞已按同一口径补过一次（只补错行，原表零改动）

产物 `Votus.integrated/calibration_20260914/`（09-14 22:18）：

```
final_integrated_classification.calibrated.tsv   8,934,718 B   (+6 列 calib_*)
calibration_changes.tsv                            590,190 B   1,524 数据行
calibration_report.md                               35,143 B
```

原表 `final_integrated_classification.tsv` md5 未变（`12e8f171…`），回滚 = 不用新文件。

### 参照表口径复算（脚本 `bl_05state.py`）

| 类别 | 亲代表（08-12） | 校准表（09-14） |
|---|---|---|
| 相容 | 13,202 | 13,202 |
| **双侧不相容（可置空）** | **1,246** | **0** |
| 单侧不相容 | 304 | 304 |
| Genus 为空 | 4,870 | **6,116** |
| Family 空、Genus 有值 | 947 | 947 |
| 属无跨参照记录（无法判） | 323 | 323 |
| 合计 | 20,892 | 20,892 |

6,116 − 4,870 = **1,246**，与置空数精确对应，双侧硬冲突清零。

### 304 单侧的拆分（`calib_action` / `calib_flag_single_ref` 实测）

```
review              278   (ncbi_only 27 / vmr_only 251)
conflict_one_side    26   (两参照都有该属，一相容一不相容)
blank             1,246
(空)             19,342
```

278 + 26 = 304，与单侧数闭合。278 review 与 26 conflict 的 Genus **全部保留**（275+3 与 22+4）。

### 置空行的实际形态

`blank` 1,246 行 = Genus 与 Species 清空、**Family 保留**。原属值留在 `calib_prev_genus`：

```
Fam=Mimiviridae       Genus=''  prev=Nitunavirus
Fam=Pithoviridae      Genus=''  prev=Chlorovirus
Fam=Schizomimiviridae Genus=''  prev=Fadolivirus
Fam=Hepaciviridae     Genus=''  prev=Fadolivirus
```

表内无 `prev_family` 列，校准**从不改写 Family**。

### 251 vmr_only 的性质

MSL41 把原 Mimiviridae 拆出新科 Hydriviridae，各工具的库还停在老科名。行内 Family 是过时科名，属级主张（如 Alphahydrivirus）未必错，判为**版本漂移**而非嵌合，故只标记不置空。典型行：

```
Family=Mimiviridae     Genus=Alphahydrivirus   VMR 定型科=Hydriviridae
Family=Phycodnaviridae Genus=Plazymidvirus     VMR 定型科=Peduoviridae
```

对侧 27 `ncbi_only` 含名称层级混入，如 `Human papillomavirus`、`uncultured partitivirus`。

## 四、覆盖缺口（实测）

| 项 | 状态 |
|---|---|
| 枸杞 05 结果层 | 已按口径 A 校准（独立 calibrated 表） |
| OneKP 05 | **未跑**。`05_Taxonomy` 无任何 `*calib*` / `*changes*`；`final_integrated_classification.tsv` 停在 2026-09-04 02:13（179,810,597 B），无 `calib_*` 列 |
| 下游植物表 | **未接**。枸杞 Lycium barbarum `All_plant.viruses_info.tsv` 2026-08-20 16:20（281,204 B）；OneKP 同名表 2026-09-11 16:22（5,161,437 B）。两者都还是校准前的分类 |
| 补丁脚本位置 | `apply_calib_A.py` / `verify_calib_A.py` 只在本地工作区 `scripts/audit/`；服务器上不存在（服务器 `~/MMPV-RNA/scripts/` 是另一套编号脚本）。`apply_calib_A.py` 的 `BASE` 硬编码枸杞路径，指到 OneKP 需改路径 |

## 五、结论

- 双侧硬冲突（1,246 行）：**已解决**，机制 + 结果两层都到位（枸杞）。
- 单侧冲突 304 行：**按设计保留**，其中 251 行属分类学版本漂移、27 行属名称层级混入、26 行两参照互相矛盾。字面上仍科属不一致，性质不是嵌合。
- Family 字段的过时科名（Mimiviridae 等）：**未处理**，本次口径不含科名映射。
- OneKP 与下游表：**未覆盖**。

## 六、待定夺

1. 是否把同一口径指向 OneKP 的 05 输出（改 `BASE` 即可，不重跑管线）。
2. 是否补一步科名映射（老科 → MSL41 新科），处理 251 行版本漂移与 1,246 行保留的老科名。
3. 947 行「有属无科」、323 行「属无跨参照记录」是否要单独定性。

---

## 七、7 个非 MSL41 科的宿主归属与影响面（2026-09-15 补，脚本 `bl_05plantfam4.py` / `bl_05ourmia.py`）

### 7.1 它们是不是植物科

三重口径（是否在 `acvirus_db/taxa.txt` / 是否被 `Plant.tsv` 的 `Virus_lineage` 收录为科 / 是否进宿主黑名单 `NON_PLANT_FAMILIES_FALLBACK`）：

| 科 | 行数 | 在 taxa.txt | 植物库认 | 宿主黑名单 | 实际宿主 |
|---|---|---|---|---|---|
| Hepaciviridae | 120 | 否 | 否 | 在 | 脊椎动物（Orthohepacivirus 76 / Pegivirus 23 / 空 21）|
| Pestiviridae | 25 | 否 | 否 | 不在 | 脊椎动物（Orthopestivirus 20 等）|
| Zimmerviridae | 19 | 否 | 否 | 不在 | Duplodnaviria，噬菌体类（Jouvirus 7 / Nesevirus 6 / Glaedevirus 4 / Bievrevirus 2）|
| Ambiguiviridae | 7 | 否 | 否 | 在 | 非植物（Alphambiguivirus 6 / Gammambiguivirus 1）|
| Autographiviridae | 5 | 否 | 否 | 不在 | 噬菌体（5 个属各 1 行）|
| **Ourmiaviridae** | **3** | 否 | 否 | 在 | **植物（Ourmiavirus 2 / Iotaourmiavirus 1）** |
| Parahypoviridae | 1 | 否 | 否 | 不在 | 真菌（Betahypovirus）|

**180 行里 178 行的属全部是非植物属，只有 Ourmiavirus 2 行（+ Iotaourmiavirus 1 行）是植物相关。**

`Plant.tsv` 无 Family 列，判植物科只能按 `Virus_lineage` 分号切分后**精确比对 token**（子串比对会把 `botourmiaviridae` 误判成 `ourmiaviridae`）。植物病毒库共认 50 个科、7,924 个属。

### 7.2 「不在 taxa.txt」不等于「科名过时」

`acvirus_db/taxa.txt`（21,541 行 / 367 科 / 2,929 属）**不是全量现行 ICTV**：现行噬菌体大科 Autographiviridae 也不在里面。它更像「ACVirus 收得到基因组的那部分分类单元」。

所以「7 个科不在 taxa.txt」要拆成三种不同性质：

- **172 行**（Hepaciviridae 120 + Pestiviridae 25 + Zimmerviridae 19 + Ambiguiviridae 7 + Parahypoviridae 1）：属在参照表里只有 NCBI 一侧有定型科（`VMR_n=0`），即**快照比 NCBI 旧**。
- **5 行** Autographiviridae：现行 ICTV 有、快照没收录，属**快照不全**。
- **3 行** Ourmiaviridae：**两套体系的真实分歧**。参照表实测 `Ourmiavirus` → `NCBI_Family=Ourmiaviridae` / `VMR_Family=Botourmiaviridae`；`Iotaourmiavirus` → NCBI 侧 `Ourmiaviridae`、VMR 侧空。

### 7.3 为什么 Ourmiaviridae 的行没被闸门置空

实测 `calib_action = conflict_one_side`，`calib_evidence = Family=Ourmiaviridae(3/4: ACVirus,CAT,diamond_lca)|Genus=Ourmiavirus(1/3: ACVirus)|NCBI定型科=Ourmiaviridae|VMR定型科=Botourmiaviridae`。

判据要求 `NCBI 与 VMR 都判不相容`。这里 NCBI 侧与行内科名一致（Ourmiaviridae），只算单侧冲突，按设计保留属与科。**不是漏判。**

### 7.4 对植物结果的实际影响

| 面 | 数字 |
|---|---|
| 7 个科的 180 行里，植物相关 | **2～3 行**（Ourmiavirus 2 + Iotaourmiavirus 1）|
| 其中落进枸杞下游植物表 | **2 行**（`Ourmiaviridae` + `Ourmiavirus`，落在 `10_Reports/All_plant.viruses_info.tsv`）|
| 1,246 blank 行里属是植物病毒属 | **25 行**（其余 1,221 行属非植物/未知）|
| 251 vmr_only 漂移行里属是植物病毒属 | **0 行**（全部非植物）|
| 现存活行里「科是病毒科但不被植物库收录 + 属是植物病毒属」 | **2 行**（就是那 2 行 Ourmiaviridae）|

版本漂移那 251 行（Mimiviridae→Hydriviridae 等）**属全部非植物**，对植物病毒结果零影响。

### 7.5 植物入库前的实际形态（校准前下游表示例）

`10_Reports/All_plant.viruses_info.tsv`（1,061 行）里科属打架的行长这样：

```
Family=Mimiviridae       Genus=Potyvirus        Species=Potyvirus rapae      tool=metabuli
Family=Marseilleviridae  Genus=Betanucleorhabdovirus  Species=Coptis betanucleorhabdovirus 1
Family=Pithoviridae      Genus=Potexvirus       Species=NA
Family=Ourmiaviridae     Genus=Ourmiavirus      Species=Ourmiavirus cucurbitae
```

前三种是**科字段错、属种正确**（属种都是植物病毒），这正是 1,246 行被置空的形态：闸门清掉了对的属，留下了错的科。若改用「按属的定型科改写科名」（两参照一致时）会更贴事实。

### 7.6 结论（对「影响大不大」的回答）

6/7 个科不是植物科，占 178/180 行，且宿主黑名单已覆盖其中 3 个。**这部分对植物病毒论文的影响可以忽略。**

唯一需要处理的是 `Ourmiaviridae`：它是植物病毒科（Ourmiavirus 类），3 行，2 行进了枸杞植物表。它不是错，是 NCBI 与 ICTV 的科名分歧（NCBI 用 Ourmiaviridae，ICTV 归 Botourmiaviridae）。

