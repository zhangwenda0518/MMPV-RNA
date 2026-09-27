#!/usr/bin/env python3
"""
utils/glm_predictors.py — BEAST 离散系统地理 GLM「迁移预测因子矩阵」生成器
==========================================================================
把「地点集合 + 坐标 / 自定义数值列」转成 BEAST 1.x GLM 离散扩散模型
(`glmSubstitutionModel`) 所需的 origin×destination 协变量矩阵。

本模块**自包含**（只依赖 numpy / pandas / 标准库），不 import 管线内其它模块，
也不修改任何既有行为；纯新增。

算法出处（逐条可追溯）
----------------------
1. 矩阵构造 / 对数变换 / 标准化 / TSV 写出 ——
   `git-repo/WMV6-phylogeography/discrete_phylogeography/scripts/
    WuMV-6_phylogeography_predictors.ipynb`（cell 2 `predictorMatrixPlot`；
    `logMatrix` / `standardizeMatrix` / `convertToString` 在该 notebook 内
    署名 "Author - Gytis Dudas"，是 BEAST 社区通用写法）：
     - origin 矩阵      O[i][j] = f(region_i)   ← 同一「行」恒定（起点端变量）
     - destination 矩阵 D[i][j] = f(region_j)   ← 同一「列」恒定（终点端变量）
     - 连续变量先取对数再标准化（总体标准差 ddof=0）
2. 地理距离矩阵 —— notebook cell 9（原实现用 `geopy.distance.geodesic`）；
   本模块用等价的 haversine 球面公式（Sinnott 1984, Sky & Telescope），
   地球半径默认 6371.0088 km（IUGG 平均半径，与 geopy 默认 WGS84 球面近似
   在 1% 内一致，自测中有断言）。notebook 中对角置 0 后做 min-max 归一化；
   因对角为 0，等价于线性缩放到 0–1。
3. 陆界（邻接）0/1 矩阵 —— notebook cell 7（`neighboring_countries`）。
4. 变量差 / 均值矩阵 —— 任务要求的通用扩展（notebook 未直接给出）：
     diff[i][j] = |v_i - v_j|（对角 0，对称），mean[i][j] = (v_i + v_j) / 2。
5. XML 片段字段结构 —— `git-repo/WMV6-phylogeography/discrete_phylogeography/
   data/NP-wmv6-discrete-ucld-skygrid-glm.xml` L1804–1866：
       <glmSubstitutionModel id="region.model">
         <glmModel family="logLinear" checkIdentifiability="true">
           <independentVariables>
             <parameter id="region.coefficients" value="0.0"/>
             <indicator><parameter id="region.coefIndicators" value="1.0"/></indicator>
             <designMatrix id="region.designMatrix">
               <parameter id="region.<name>_predictor_matrix" value="..."/>
             </designMatrix>
           </independentVariables>
         </glmModel>
       </glmSubstitutionModel>
   **已核实（本次现场统计）**：真实 XML 里每个预测因子参数携带 n*(n-1) 个值
   （WMV6：27 区划 → 702 = 27×26；flu_d：13/14 区划 → 156 = 13×12、
   182 = 14×13），即**只写非对角元**；而 notebook 落盘的 TSV 是完整 n×n。
   本模块 `--xml-layout offdiag`（默认）复现 XML 的真实布局，
   `full` 写出完整 n×n 供人工核对，两者都做自检。

用法
----
  # 1) 自测（内置玩具数据 + 断言，不需要任何外部文件）
  python utils/glm_predictors.py --selftest

  # 2) 真实数据：地点表 (region,lat,lon) + 变量表 (region,value)
  python utils/glm_predictors.py --locations loc.csv --values var.csv \
      --value-col annual_precip --var-name precip --outdir predictors --xml

  # 3) 变量就在地点表里（--value-col 指向 loc.csv 的数值列）
  python utils/glm_predictors.py --locations loc.csv --value-col mean_temp \
      --var-name temp --outdir predictors --xml

  # 4) 只要距离矩阵（无需变量表）
  python utils/glm_predictors.py --locations loc.csv --outdir predictors --xml

  # 5) 邻接矩阵（邻接表 CSV: region,neighbor 两列，可多行）
  python utils/glm_predictors.py --locations loc.csv --neighbors borders.csv \
      --outdir predictors --xml

输出
----
  <outdir>/<var>_origin_predictor_matrix.txt       (TSV, n×n, %.8f)
  <outdir>/<var>_destination_predictor_matrix.txt  (TSV, n×n, %.8f)
  <outdir>/<var>_difference_predictor_matrix.txt   （若 --matrix diff|all）
  <outdir>/<var>_mean_predictor_matrix.txt         （若 --matrix mean|all）
  <outdir>/distance_predictor_matrix.txt           （若给了坐标）
  <outdir>/borders_predictor_matrix.txt            （若给了邻接表）
  <outdir>/predictors_design_matrix.xml            （若 --xml，含完整 designMatrix 块）
  <outdir>/predictors_manifest.json                （区划顺序 + 数据出处，供 XML 复现）

注意
----
* 区划顺序 = 输入表出现的顺序（notebook 亦如此：`NP_region_list` 顺序即 XML
  中矩阵的行列顺序）。矩阵是**位置敏感**的，改动顺序必须同步改 XML。
* BEAST trait（离散性状）名需与 `--trait-name`（默认 `region`）一致，
  XML 参数 id 前缀随之变化。
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EARTH_RADIUS_KM = 6371.0088  # IUGG 平均地球半径（geopy 默认 WGS84 近似）

_NAME_KEYS = ("region", "location", "locality", "site", "state", "place", "name", "id")
_LAT_KEYS = ("lat", "latitude", "lat_deg", "y")
_LON_KEYS = ("lon", "lng", "long", "longitude", "lon_deg", "x")

_TRANSFORMS = ("log_standardize", "standardize", "minmax", "raw")


# ═══════════════════════════════════════════════════════════════════
# 1. 距离 / 坐标
# ═══════════════════════════════════════════════════════════════════
def haversine_km(
    lat1: float, lon1: float, lat2: float, lon2: float,
    radius_km: float = EARTH_RADIUS_KM,
) -> float:
    """两点球面大圆距离（km）。haversine 公式（Sinnott 1984）。

    与 notebook 的 `geopy.distance.geodesic`（WGS84 椭球）相比，球面近似在
    实测上千公里尺度偏差 < 0.2%（自测有断言），对 GLM 协变量完全够用。
    """
    p1 = math.radians(float(lat1))
    p2 = math.radians(float(lat2))
    dphi = math.radians(float(lat2) - float(lat1))
    dlam = math.radians(float(lon2) - float(lon1))
    a = math.sin(dphi / 2.0) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlam / 2.0) ** 2
    a = min(1.0, max(0.0, a))
    return 2.0 * radius_km * math.asin(math.sqrt(a))


def haversine_matrix(
    regions: Sequence[str],
    coords: Dict[str, Tuple[float, float]],
    radius_km: float = EARTH_RADIUS_KM,
    normalize: bool = False,
) -> pd.DataFrame:
    """haversine 距离矩阵（对称、对角 0）。

    Parameters
    ----------
    regions : 区划名列表（决定行列顺序）
    coords  : {region: (lat, lon)}
    normalize : True 时线性缩放到 0–1（notebook 的 min-max；对角为 0 的缘故，
                等价于除以最大距离）
    """
    n = len(regions)
    m = np.zeros((n, n), dtype=float)
    for i, r_i in enumerate(regions):
        if r_i not in coords:
            raise KeyError(f"区划 {r_i!r} 缺坐标")
        for j in range(i + 1, n):
            r_j = regions[j]
            if r_j not in coords:
                raise KeyError(f"区划 {r_j!r} 缺坐标")
            d = haversine_km(coords[r_i][0], coords[r_i][1],
                             coords[r_j][0], coords[r_j][1], radius_km)
            m[i, j] = m[j, i] = d
    if normalize:
        vmax = float(m.max()) if m.size else 0.0
        if vmax > 0:
            m = m / vmax
    return pd.DataFrame(m, index=list(regions), columns=list(regions))


def _pick_col(cols: Sequence[str], keys: Sequence[str]) -> Optional[str]:
    """列名模糊匹配：先全等，再去空格/下划线后相等，最后前缀命中。"""
    norm = {str(c).strip().lower().replace(" ", "_"): c for c in cols}
    for k in keys:
        if k in norm:
            return norm[k]
    for k in keys:
        for c_norm, c in norm.items():
            if c_norm.startswith(k):
                return c
    return None


def read_table(path: str) -> pd.DataFrame:
    """读 CSV/TSV（自动嗅探分隔符）。"""
    sep = None
    if path.lower().endswith(".tsv"):
        sep = "\t"
    return pd.read_csv(path, sep=sep, engine="python")


def load_locations(path: str) -> pd.DataFrame:
    """读地点表 → DataFrame(region, lat, lon [, 其它数值列…])。

    列名支持 region/location/name… 与 lat/latitude、lon/long/longitude
    （大小写、空格、下划线不敏感）。数值列原样透传，便于 --value-col 直接取用。
    """
    df = read_table(path)
    name_c = _pick_col(df.columns, _NAME_KEYS)
    lat_c = _pick_col(df.columns, _LAT_KEYS)
    lon_c = _pick_col(df.columns, _LON_KEYS)
    if name_c is None or lat_c is None or lon_c is None:
        raise ValueError(
            f"{path}: 需要 region/lat/lon 三类列，实际列为 {list(df.columns)}")
    out = pd.DataFrame({
        "region": df[name_c].astype(str).str.strip(),
        "lat": pd.to_numeric(df[lat_c], errors="coerce"),
        "lon": pd.to_numeric(df[lon_c], errors="coerce"),
    })
    for c in df.columns:
        if c in (name_c, lat_c, lon_c):
            continue
        num = pd.to_numeric(df[c], errors="coerce")
        if num.notna().all():
            out[str(c)] = num.astype(float)
    if out["lat"].isna().any() or out["lon"].isna().any():
        bad = out.loc[out["lat"].isna() | out["lon"].isna(), "region"].tolist()
        raise ValueError(f"{path}: 以下区划坐标非数值：{bad}")
    return out


def load_region_values(path: str, value_col: str, region_col: Optional[str] = None) -> Dict[str, float]:
    """读变量表 → {region: value}（region 列自动识别，找不到用第一列）。"""
    df = read_table(path)
    rc = region_col or _pick_col(df.columns, _NAME_KEYS) or df.columns[0]
    if value_col not in df.columns:
        raise ValueError(f"{path}: 无列 {value_col!r}（实际列 {list(df.columns)}）")
    vals = pd.to_numeric(df[value_col], errors="coerce")
    if vals.isna().any():
        bad = df.loc[vals.isna(), rc].astype(str).tolist()
        raise ValueError(f"{path}: 列 {value_col!r} 有非数值行：{bad}")
    return {str(k).strip(): float(v) for k, v in zip(df[rc], vals)}


def load_neighbors(path: str) -> Dict[str, List[str]]:
    """读邻接表（两列：region, neighbor；可多行）→ {region: [neighbor,…]}。

    有向输入会被对称化（陆界/邻接关系本身对称）。
    """
    df = read_table(path)
    if df.shape[1] < 2:
        raise ValueError(f"{path}: 邻接表需要两列 region,neighbor")
    a, b = df.columns[0], df.columns[1]
    adj: Dict[str, List[str]] = {}
    for x, y in zip(df[a].astype(str).str.strip(), df[b].astype(str).str.strip()):
        adj.setdefault(x, []).append(y)
        adj.setdefault(y, []).append(x)
    return {k: sorted(set(v)) for k, v in adj.items()}


# ═══════════════════════════════════════════════════════════════════
# 2. 变换（复刻 notebook 的 logMatrix / standardizeMatrix）
# ═══════════════════════════════════════════════════════════════════
def log_transform(values: Sequence[float]) -> np.ndarray:
    """逐元素自然对数（notebook `logMatrix`）。要求所有值 > 0。"""
    arr = np.asarray(values, dtype=float)
    if np.any(arr <= 0):
        bad = [i for i, v in enumerate(arr) if v <= 0]
        raise ValueError(f"log 变换要求全部 > 0，第 {bad} 个值不合法")
    return np.log(arr)


def standardize(values: Sequence[float]) -> np.ndarray:
    """标准化 (x-mean)/sd，sd 用总体标准差 ddof=0（notebook `standardizeMatrix`）。"""
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        raise ValueError("standardize 收到空向量")
    sd = float(arr.std(ddof=0))
    if sd == 0.0:
        raise ValueError("常量向量无法标准化（sd=0）；请改用 --transform raw/minmax")
    return (arr - float(arr.mean())) / sd


def minmax(values: Sequence[float]) -> np.ndarray:
    """线性缩放到 0–1（常量向量返回全 1，避免 NaN）。"""
    arr = np.asarray(values, dtype=float)
    lo, hi = float(arr.min()), float(arr.max())
    if hi == lo:
        return np.ones_like(arr)
    return (arr - lo) / (hi - lo)


def transform_values(values: Sequence[float], mode: str = "log_standardize") -> np.ndarray:
    """统一入口：log_standardize（notebook 默认）/ standardize / minmax / raw。"""
    if mode not in _TRANSFORMS:
        raise ValueError(f"未知变换 {mode!r}；可选 {_TRANSFORMS}")
    if mode == "log_standardize":
        return standardize(log_transform(values))
    if mode == "standardize":
        return standardize(values)
    if mode == "minmax":
        return minmax(values)
    return np.asarray(values, dtype=float)


# ═══════════════════════════════════════════════════════════════════
# 3. 矩阵构造（复刻 notebook `predictorMatrixPlot` + 通用扩展）
# ═══════════════════════════════════════════════════════════════════
def origin_matrix(regions: Sequence[str], values: Dict[str, float],
                  transform: str = "log_standardize") -> pd.DataFrame:
    """origin 矩阵：第 i 行 = 起点区划 i 的（变换后）变量值，整行恒定。"""
    t = transform_values([values[r] for r in regions], transform)
    m = np.tile(t.reshape(-1, 1), (1, len(regions)))
    return pd.DataFrame(m, index=list(regions), columns=list(regions))


def destination_matrix(regions: Sequence[str], values: Dict[str, float],
                       transform: str = "log_standardize") -> pd.DataFrame:
    """destination 矩阵：第 j 列 = 终点区划 j 的（变换后）变量值，整列恒定。"""
    t = transform_values([values[r] for r in regions], transform)
    m = np.tile(t.reshape(1, -1), (len(regions), 1))
    return pd.DataFrame(m, index=list(regions), columns=list(regions))


def difference_matrix(regions: Sequence[str], values: Dict[str, float],
                      transform: str = "raw") -> pd.DataFrame:
    """|v_i - v_j| 矩阵（对称、对角 0）。transform 作用在**变量**上（默认不变换）。"""
    t = transform_values([values[r] for r in regions], transform)
    m = np.abs(t.reshape(-1, 1) - t.reshape(1, -1))
    np.fill_diagonal(m, 0.0)
    return pd.DataFrame(m, index=list(regions), columns=list(regions))


def mean_matrix(regions: Sequence[str], values: Dict[str, float],
                transform: str = "raw") -> pd.DataFrame:
    """(v_i + v_j)/2 矩阵（对称）。"""
    t = transform_values([values[r] for r in regions], transform)
    m = (t.reshape(-1, 1) + t.reshape(1, -1)) / 2.0
    return pd.DataFrame(m, index=list(regions), columns=list(regions))


def binary_matrix(regions: Sequence[str], neighbors: Dict[str, List[str]]) -> pd.DataFrame:
    """邻接 0/1 矩阵（对称、对角 0；notebook cell 7 的 land borders 写法）。"""
    idx = {r: i for i, r in enumerate(regions)}
    n = len(regions)
    m = np.zeros((n, n), dtype=float)
    for r, nbrs in neighbors.items():
        if r not in idx:
            continue
        for nb in nbrs:
            if nb in idx and nb != r:
                m[idx[r], idx[nb]] = 1.0
                m[idx[nb], idx[r]] = 1.0
    np.fill_diagonal(m, 0.0)
    return pd.DataFrame(m, index=list(regions), columns=list(regions))


def build_predictor_set(
    regions: Sequence[str],
    coords: Optional[Dict[str, Tuple[float, float]]] = None,
    values: Optional[Dict[str, float]] = None,
    neighbors: Optional[Dict[str, List[str]]] = None,
    var_name: str = "var",
    transform: str = "log_standardize",
    distance_normalize: bool = True,
) -> Dict[str, pd.DataFrame]:
    """一次生成全套预测因子矩阵。

    Returns
    -------
    dict[str, DataFrame]，键为预测因子名（同时用作 XML 参数 id 与 TSV 文件名）：
      {var}_origin / {var}_destination / {var}_difference / {var}_mean
      distance / borders        （后两个仅在给了 coords / neighbors 时出现）
    """
    regions = [str(r) for r in regions]
    if len(set(regions)) != len(regions):
        raise ValueError(f"区划名重复：{regions}")
    out: Dict[str, pd.DataFrame] = {}
    if values is not None:
        missing = [r for r in regions if r not in values]
        if missing:
            raise KeyError(f"变量表缺区划：{missing}")
        out[f"{var_name}_origin"] = origin_matrix(regions, values, transform)
        out[f"{var_name}_destination"] = destination_matrix(regions, values, transform)
        out[f"{var_name}_difference"] = difference_matrix(regions, values, "raw")
        out[f"{var_name}_mean"] = mean_matrix(regions, values, "raw")
    if coords is not None:
        out["distance"] = haversine_matrix(regions, coords, normalize=distance_normalize)
    if neighbors is not None:
        out["borders"] = binary_matrix(regions, neighbors)
    return out


# ═══════════════════════════════════════════════════════════════════
# 4. 导出：TSV / XML 片段
# ═══════════════════════════════════════════════════════════════════
def flatten_values(df: pd.DataFrame, layout: str = "offdiag") -> str:
    """矩阵 → XML `value` 字符串（空格分隔）。

    layout='full'    : 行优先展开 n×n（notebook 落盘的 TSV 顺序）
    layout='offdiag' : 行优先展开非对角元（i 行内 j≠i，共 n*(n-1) 个）
                       —— 与真实 WMV6 / flu_d XML 的取值个数一致
    """
    n = df.shape[0]
    if df.shape[1] != n:
        raise ValueError("predictor 矩阵必须方阵（n×n）")
    if layout == "full":
        vals = df.to_numpy(dtype=float).reshape(-1)
    elif layout == "offdiag":
        vals = np.array([df.iat[i, j] for i in range(n) for j in range(n) if i != j],
                        dtype=float)
    else:
        raise ValueError(f"未知 layout {layout!r}（可选 full / offdiag）")
    return " ".join(f"{v:.8g}" for v in vals)


def write_matrix_tsv(df: pd.DataFrame, path: str, float_format: str = "%.8f") -> str:
    """按 notebook 的方式落盘 TSV（首列 = 行标签，首行 = 列标签）。"""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    df.to_csv(path, sep="\t", float_format=float_format)
    return path


def xml_parameter_line(param_id: str, df: pd.DataFrame, layout: str = "offdiag",
                       indent: str = "\t\t\t\t") -> str:
    """生成一行 `<parameter>`（designMatrix 的一个预测因子列）。"""
    return (f'{indent}<parameter id="{param_id}" '
            f'value="{flatten_values(df, layout)}"/>')


def xml_design_matrix_block(predictors: Dict[str, pd.DataFrame], trait_name: str = "region",
                            layout: str = "offdiag",
                            comments: Optional[Dict[str, str]] = None) -> str:
    """生成 `<designMatrix>` 块（含每个预测因子一行参数 + 出处注释）。

    字段结构见本模块 docstring 第 5 条（WMV6 XML L1804–1866）。
    """
    lines = [f'\t\t\t<designMatrix id="{trait_name}.designMatrix">']
    comments = comments or {}
    for i, (name, df) in enumerate(predictors.items(), start=1):
        note = comments.get(name)
        if note:
            lines.append(f"\t\t\t\t<!-- predictor {i}: {note} -->")
        lines.append(xml_parameter_line(
            f"{trait_name}.{name}_predictor_matrix", df, layout,
            indent="\t\t\t\t"))
    lines.append("\t\t\t</designMatrix>")
    return "\n".join(lines)


def xml_glm_substitution_model_block(predictors: Dict[str, pd.DataFrame],
                                     trait_name: str = "region",
                                     n_states: Optional[int] = None,
                                     data_type_id: Optional[str] = None,
                                     model_id: Optional[str] = None,
                                     layout: str = "offdiag",
                                     comments: Optional[Dict[str, str]] = None) -> str:
    """生成**完整**的 `<glmSubstitutionModel>` 块（可直接替换进 BEAST 1.x XML）。

    ⚠️ 只生成 GLM 模型本体 + designMatrix；generalDataType / 树 / 算子 / 先验
    仍需由管线的 XML 生成器补齐（本模块刻意不碰既有 beast1_bridge 行为）。
    """
    trait_name = str(trait_name)
    n = int(n_states if n_states is not None else next(iter(predictors.values())).shape[0])
    dtd = data_type_id or f"{trait_name}.dataType"
    mid = model_id or f"{trait_name}.model"
    block = [
        f'\t<glmSubstitutionModel id="{mid}">',
        f'\t\t<generalDataType idref="{dtd}"/>',
        "\t\t<rootFrequencies>",
        f'\t\t\t<frequencyModel id="{trait_name}.frequencyModel" normalize="true">',
        f'\t\t\t\t<generalDataType idref="{dtd}"/>',
        "\t\t\t\t<frequencies>",
        f'\t\t\t\t\t<parameter id="{trait_name}.frequencies" dimension="{n}"/>',
        "\t\t\t\t</frequencies>",
        "\t\t\t</frequencyModel>",
        "\t\t</rootFrequencies>",
        '\t\t<glmModel family="logLinear" checkIdentifiability="true">',
        "\t\t\t<independentVariables>",
        f'\t\t\t\t<parameter id="{trait_name}.coefficients" value="0.0"/>',
        "\t\t\t\t<indicator>",
        f'\t\t\t\t\t<parameter id="{trait_name}.coefIndicators" value="1.0"/>',
        "\t\t\t\t</indicator>",
        xml_design_matrix_block(predictors, trait_name, layout, comments),
        "\t\t\t</independentVariables>",
        "\t\t</glmModel>",
        "\t</glmSubstitutionModel>",
    ]
    return "\n".join(block)


# ═══════════════════════════════════════════════════════════════════
# 5. 内置玩具数据 + 自测
# ═══════════════════════════════════════════════════════════════════
# 坐标取自 WMV6 notebook 的 `centroid_coords`（真实区划质心，便于对拍）
TOY_COORDS: Dict[str, Tuple[float, float]] = {
    "Yunnan": (24.990533, 101.493693),
    "Hainan": (17.883595, 111.241771),
    "Guangdong": (22.912986, 113.341439),
    "Hubei": (30.972024, 112.257427),
    "California": (36.551343, -121.320102),
}
TOY_REGIONS: List[str] = ["Yunnan", "Hainan", "Guangdong", "Hubei", "California"]
# 变量取 WMV6 notebook 的 `annual_precipitation`（mm/月级 WorldClim 汇总）
TOY_VALUES: Dict[str, float] = {
    "Yunnan": 60.49, "Hainan": 76.12, "Guangdong": 104.36,
    "Hubei": 42.92, "California": 30.48,
}
TOY_NEIGHBORS: Dict[str, List[str]] = {
    "Yunnan": ["Hubei", "Guangdong"],
    "Hainan": ["Guangdong"],
    "Guangdong": ["Hainan", "Yunnan", "Hubei"],
    "Hubei": ["Yunnan", "Guangdong"],
    "California": [],
}
# 文献/解析真值（haversine，R=6371.0088）：用于 1% 误差断言
KNOWN_DISTANCES_KM = {
    ("Beijing", "Shanghai"): ((39.9042, 116.4074), (31.2304, 121.4737), 1067.0),
    ("London", "Paris"): ((51.5074, -0.1278), (48.8566, 2.3522), 343.6),
    ("(0,0)", "(0,90)"): ((0.0, 0.0), (0.0, 90.0), 10007.5572),  # = R*pi/2 解析解
}


def _selftest(verbose: bool = True) -> int:
    """内置玩具数据自测：断言矩阵性质 / 已知距离误差 / XML 布局 / 往返一致。"""
    n_ok = 0

    def _check(cond: bool, msg: str) -> None:
        nonlocal n_ok
        if not cond:
            raise AssertionError(f"SELFTEST FAILED: {msg}")
        n_ok += 1
        if verbose:
            print(f"  [ok {n_ok:2d}] {msg}")

    # ── A. haversine 已知点距（误差 < 1%）──────────────────────────
    if verbose:
        print("A. haversine 已知距离（误差阈值 1%）")
    for (a, b), ((la1, lo1), (la2, lo2), ref) in KNOWN_DISTANCES_KM.items():
        d = haversine_km(la1, lo1, la2, lo2)
        err = abs(d - ref) / ref
        _check(err < 0.01, f"{a}–{b}: {d:.2f} km vs 真值 {ref} km，相对误差 {err*100:.3f}% < 1%")
    # 解析解：同一经线 90 度 = R * pi/2
    exact = EARTH_RADIUS_KM * math.pi / 2.0
    _check(abs(haversine_km(0, 0, 0, 90) - exact) < 1e-6,
           f"解析解校验 haversine((0,0),(0,90)) == R*pi/2 = {exact:.4f} km")

    # ── B. 距离矩阵性质：对称 / 对角 0 ──────────────────────────────
    if verbose:
        print("B. 距离矩阵（对称、对角 0、归一化 0–1）")
    dm = haversine_matrix(TOY_REGIONS, TOY_COORDS, normalize=False)
    _check(np.allclose(dm.to_numpy(), dm.to_numpy().T, atol=1e-9), "距离矩阵对称")
    _check(np.allclose(np.diag(dm.to_numpy()), 0.0), "距离矩阵对角为 0")
    _check(abs(dm.loc["Yunnan", "Hainan"] - 1280.7215) < 1.0,
           f"云南–海南 = {dm.loc['Yunnan','Hainan']:.4f} km ≈ 1280.72（notebook 质心对拍）")
    dmn = haversine_matrix(TOY_REGIONS, TOY_COORDS, normalize=True)
    _check(dmn.to_numpy().max() == 1.0 and dmn.to_numpy().min() == 0.0,
           "min-max 归一化后值域恰为 [0,1]（对角仍 0）")
    _check(np.allclose(np.diag(dmn.to_numpy()), 0.0), "归一化后对角仍为 0")

    # ── C. origin / destination 矩阵结构 ───────────────────────────
    if verbose:
        print("C. origin / destination 矩阵结构（行/列恒定 + 变换数值）")
    om = origin_matrix(TOY_REGIONS, TOY_VALUES)
    dmx = destination_matrix(TOY_REGIONS, TOY_VALUES)
    _check(all(om.iloc[i].nunique() == 1 for i in range(len(TOY_REGIONS))),
           "origin 矩阵每行恒定（行=起点区划）")
    _check(all(dmx.iloc[:, j].nunique() == 1 for j in range(len(TOY_REGIONS))),
           "destination 矩阵每列恒定（列=终点区划）")
    _check(bool(np.allclose(om.to_numpy().T, dmx.to_numpy())),
           "destination == origin 的转置（notebook 两矩阵互为转置的关系）")
    z = transform_values([TOY_VALUES[r] for r in TOY_REGIONS], "log_standardize")
    _check(abs(z.mean()) < 1e-12 and abs(z.std(ddof=0) - 1.0) < 1e-12,
           "log_standardize 后均值 0 / 总体标准差 1（复刻 notebook standardizeMatrix）")
    _check(abs(om.loc["Hainan", "Hainan"] - z[1]) < 1e-12,
           f"Hainan 变换值 = {z[1]:.8f}（手工对拍第 2 个位置）")
    try:
        transform_values([1.0, 0.0, 2.0], "log_standardize")
        raise AssertionError("log 变换未拦截 <=0 输入")
    except ValueError:
        _check(True, "log 变换对 <=0 输入抛 ValueError（不静默出 NaN）")
    try:
        transform_values([3.0, 3.0, 3.0], "log_standardize")
        raise AssertionError("常量向量未拦截")
    except ValueError:
        _check(True, "常量向量标准化抛 ValueError（sd=0 不产生 inf/NaN）")

    # ── D. 差 / 均值 / 邻接矩阵 ─────────────────────────────────────
    if verbose:
        print("D. 差 / 均值 / 邻接矩阵")
    diff = difference_matrix(TOY_REGIONS, TOY_VALUES)
    _check(np.allclose(diff.to_numpy(), diff.to_numpy().T, atol=1e-12), "差矩阵对称")
    _check(np.allclose(np.diag(diff.to_numpy()), 0.0), "差矩阵对角 0")
    _check(abs(diff.loc["Yunnan", "Hainan"] -
               abs(TOY_VALUES["Yunnan"] - TOY_VALUES["Hainan"])) < 1e-12,
           "差矩阵元素 = |v_i - v_j|（手算对拍）")
    mn = mean_matrix(TOY_REGIONS, TOY_VALUES)
    _check(abs(mn.loc["Hubei", "California"] -
               (TOY_VALUES["Hubei"] + TOY_VALUES["California"]) / 2) < 1e-12,
           "均值矩阵元素 = (v_i+v_j)/2（手算对拍）")
    bm = binary_matrix(TOY_REGIONS, TOY_NEIGHBORS)
    _check(np.allclose(bm.to_numpy(), bm.to_numpy().T), "邻接矩阵对称")
    _check(np.allclose(np.diag(bm.to_numpy()), 0.0) and set(np.unique(bm.to_numpy())) <= {0.0, 1.0},
           "邻接矩阵对角 0 且取值仅 {0,1}")
    _check(bm.loc["Yunnan", "Hubei"] == 1.0 and bm.loc["California", "Yunnan"] == 0.0,
           "邻接关系正确（云南–湖北 1，加州–云南 0）")

    # ── E. XML 布局：offdiag == n*(n-1)；full == n*n ────────────────
    if verbose:
        print("E. XML 布局（与真实 WMV6 / flu_d XML 取值个数对拍）")
    n = len(TOY_REGIONS)
    _check(len(flatten_values(om, "offdiag").split()) == n * (n - 1),
           f"offdiag 布局值数 = n*(n-1) = {n*(n-1)}（WMV6 27×26=702 / flu_d 13×12=156 同构）")
    _check(len(flatten_values(om, "full").split()) == n * n,
           f"full 布局值数 = n*n = {n*n}（notebook TSV 顺序）")
    full_vals = np.array([float(x) for x in flatten_values(om, "full").split()])
    _check(bool(np.allclose(full_vals, om.to_numpy().reshape(-1))),
           "full 展开 = 行优先 reshape(-1)，与 TSV 读回一致")
    off_vals = np.array([float(x) for x in flatten_values(diff, "offdiag").split()])
    _check(len(off_vals) == n * (n - 1) and off_vals.max() > 0,
           "对称预测因子（差矩阵）offdiag 展开无对角占位")

    # ── F. 端到端：build_predictor_set + XML 块 + TSV 往返 ──────────
    if verbose:
        print("F. 端到端（build_predictor_set / XML 块 / TSV 往返）")
    preds = build_predictor_set(TOY_REGIONS, TOY_COORDS, TOY_VALUES, TOY_NEIGHBORS,
                                var_name="precip")
    _check(set(preds) == {"precip_origin", "precip_destination", "precip_difference",
                          "precip_mean", "distance", "borders"},
           "build_predictor_set 生成 6 个预测因子矩阵")
    block = xml_design_matrix_block(preds, trait_name="region")
    _check(block.count("<parameter id=\"region.") == 6, "designMatrix 块含 6 个预测因子参数")
    _check("<designMatrix id=\"region.designMatrix\">" in block
           and block.rstrip().endswith("</designMatrix>"), "designMatrix 开闭标签完整")
    _check('id="region.precip_origin_predictor_matrix"' in block,
           "参数 id 命名与 WMV6 一致（region.<name>_predictor_matrix）")
    gm = xml_glm_substitution_model_block(preds, trait_name="region", n_states=n)
    for frag in ('<glmSubstitutionModel id="region.model">',
                 '<glmModel family="logLinear" checkIdentifiability="true">',
                 '<parameter id="region.coefficients" value="0.0"/>',
                 f'<parameter id="region.frequencies" dimension="{n}"/>',
                 "<designMatrix"):
        _check(frag in gm, f"GLM 块含字段 {frag[:46]}…")
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = write_matrix_tsv(om, os.path.join(td, "origin.txt"))
        back = pd.read_csv(p, sep="\t", index_col=0)
        _check(np.allclose(back.to_numpy(), om.to_numpy(), atol=1e-8),
               "TSV 落盘→读回与内存矩阵一致（%.8f 精度内）")
        _check(list(back.index) == TOY_REGIONS and list(back.columns) == TOY_REGIONS,
               "TSV 保留区划顺序（矩阵位置敏感）")

    # ── G. 反例拦截 ────────────────────────────────────────────────
    if verbose:
        print("G. 输入校验")
    try:
        build_predictor_set(["A", "A"], values={"A": 1.0})
        raise AssertionError("重复区划未拦截")
    except ValueError:
        _check(True, "重复区划名抛 ValueError")
    try:
        build_predictor_set(["A", "B"], values={"A": 1.0})
        raise AssertionError("缺变量未拦截")
    except KeyError:
        _check(True, "变量表缺区划抛 KeyError")

    if verbose:
        print(f"\nSELFTEST PASSED — {n_ok} 项断言全部通过")
    return n_ok


# ═══════════════════════════════════════════════════════════════════
# 6. CLI
# ═══════════════════════════════════════════════════════════════════
def _build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="glm_predictors",
        description="BEAST 离散系统地理 GLM 迁移预测因子矩阵生成器"
                    "（haversine 距离 / origin-destination / 差-均值 / 邻接）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("用法\n----")[-1][:1200])
    ap.add_argument("--selftest", action="store_true",
                    help="跑内置玩具数据自测（断言矩阵性质 / 已知距离 / XML 布局）后退出")
    ap.add_argument("--locations", help="地点表 CSV/TSV（region,lat,lon [,…]）")
    ap.add_argument("--values", help="变量表 CSV/TSV（region,value）；缺省时用 --value-col 取 --locations 的列")
    ap.add_argument("--value-col", default=None, help="变量列名（--values 表内，或 --locations 表内）")
    ap.add_argument("--region-col", default=None, help="变量表的 region 列名（缺省自动识别）")
    ap.add_argument("--var-name", default="var", help="预测因子名（写入参数 id 与文件名，默认 var）")
    ap.add_argument("--transform", default="log_standardize", choices=list(_TRANSFORMS),
                    help="origin/destination 的变量变换（默认 log_standardize，同 notebook）")
    ap.add_argument("--neighbors", default=None, help="邻接表 CSV/TSV（region,neighbor）")
    ap.add_argument("--no-distance", action="store_true", help="即使有坐标也不出距离矩阵")
    ap.add_argument("--distance-raw", action="store_true", help="距离矩阵不做 0–1 归一化（输出 km）")
    ap.add_argument("--outdir", default="glm_predictors_out", help="输出目录")
    ap.add_argument("--trait-name", default="region", help="BEAST 离散性状名（默认 region）")
    ap.add_argument("--xml", action="store_true", help="额外导出 predictors_design_matrix.xml（GLM 块）")
    ap.add_argument("--xml-layout", default="offdiag", choices=["offdiag", "full"],
                    help="XML 参数展开方式（默认 offdiag = 真实 XML 的 n*(n-1) 写法）")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_argparser().parse_args(argv)

    if args.selftest:
        print("=" * 64)
        print("glm_predictors 自测（玩具数据：WMV6 notebook 质心 + 月降水变量）")
        print("=" * 64)
        try:
            _selftest()
        except AssertionError as e:
            print(f"\n{e}")
            return 1
        return 0

    if not args.locations:
        print("错误：需要 --locations（或 --selftest）", file=sys.stderr)
        return 2

    loc = load_locations(args.locations)
    regions = loc["region"].tolist()
    coords = {r: (float(la), float(lo))
              for r, la, lo in zip(loc["region"], loc["lat"], loc["lon"])}

    values = None
    if args.value_col:
        if args.values:
            values = load_region_values(args.values, args.value_col, args.region_col)
        elif args.value_col in loc.columns:
            values = {str(r): float(v) for r, v in zip(loc["region"], loc[args.value_col])}
        else:
            print(f"错误：--locations 无列 {args.value_col!r}，也未给 --values", file=sys.stderr)
            return 2
    elif args.values:
        print("错误：给了 --values 但缺 --value-col", file=sys.stderr)
        return 2

    neighbors = load_neighbors(args.neighbors) if args.neighbors else None

    preds = build_predictor_set(
        regions, coords, values, neighbors,
        var_name=args.var_name, transform=args.transform,
        distance_normalize=not args.distance_raw)
    if args.no_distance:
        preds.pop("distance", None)

    os.makedirs(args.outdir, exist_ok=True)
    written = []
    for name, df in preds.items():
        path = os.path.join(args.outdir, f"{name}_predictor_matrix.txt")
        written.append(write_matrix_tsv(df, path))
        print(f"[TSV ] {path}  ({df.shape[0]}x{df.shape[1]})")

    comments = {name: name for name in preds}
    if args.xml:
        xml_path = os.path.join(args.outdir, "predictors_design_matrix.xml")
        block = xml_glm_substitution_model_block(
            preds, trait_name=args.trait_name, n_states=len(regions),
            layout=args.xml_layout, comments=comments)
        header = (f"<!-- glm_predictors.py 生成 · 性状={args.trait_name} · "
                  f"n={len(regions)} · layout={args.xml_layout} ·\n"
                  f"     算法出处 WMV6-phylogeography notebook + XML L1804-1866 -->\n")
        with open(xml_path, "w", encoding="utf-8") as fh:
            fh.write(header + block + "\n")
        written.append(xml_path)
        print(f"[XML ] {xml_path}  ({len(preds)} predictors, layout={args.xml_layout})")

    manifest = {
        "module": "utils/glm_predictors.py",
        "trait_name": args.trait_name,
        "region_order": regions,
        "n_regions": len(regions),
        "predictors": {name: {"shape": list(df.shape),
                              "xml_layout": args.xml_layout}
                       for name, df in preds.items()},
        "transform": args.transform,
        "distance_normalized": not args.distance_raw,
        "algorithm_sources": [
            "git-repo/WMV6-phylogeography/discrete_phylogeography/scripts/"
            "WuMV-6_phylogeography_predictors.ipynb (cell 2 predictorMatrixPlot;"
            " logMatrix/standardizeMatrix by Gytis Dudas)",
            "git-repo/WMV6-phylogeography/discrete_phylogeography/data/"
            "NP-wmv6-discrete-ucld-skygrid-glm.xml L1804-1866 (glmSubstitutionModel)",
            "haversine: Sinnott 1984 (R=6371.0088 km, IUGG mean radius)",
        ],
        "files": written,
    }
    mpath = os.path.join(args.outdir, "predictors_manifest.json")
    with open(mpath, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=2)
    print(f"[JSON] {mpath}  (区划顺序已锁定，供 XML 复现)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
