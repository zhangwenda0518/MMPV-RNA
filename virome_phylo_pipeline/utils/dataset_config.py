"""dataset_config.py — virome_phylo_pipeline 通用数据集配置加载器

配置文件: datasets.yaml (与本包同级的 pipeline 根目录)
新数据集只需在 yaml 的 datasets 下加条目, 无需改任何代码。

路径占位符与环境变量（2026-09-22 增加, 目标是"换环境只改一处"）：
  {repo}       仓库根目录（本文件上溯 3 级），如 {repo}/biosoft/pypopart/src
               —— 本地与服务器都在仓库内, 用它可两边通用
  {base}       yaml 顶层 base 字段（数据根目录）
  {results}    yaml 顶层 results 字段（运行产物根目录）
  {dir}        数据集自身的 dir 字段（供 fasta/metadata/haplo.* 等复用, 避免重复长路径）
  ${VAR}       环境变量；未定义则展开为空串
  ${VAR:-默认} 环境变量；未定义时用默认值
               —— 换环境可 export MMPV_DATA_BASE=... 而不动 yaml

加载顺序：环境变量展开 → 占位符替换。故 ${MMPV_X:-{repo}/y} 可正常工作。
"""
import os
import re

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_CONFIG_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "datasets.yaml")

# 支持 ${VAR} / ${VAR:-默认}；默认值允许一层 {占位符} 嵌套，如
# ${MMPV_RESULTS:-{repo}/virome_phylo_pipeline/phylo_results}
_ENV_PAT = re.compile(r"\$\{([A-Za-z_]\w*)(?::-((?:[^{}]|\{[^{}]*\})*))?\}")

_cache = None


def _expand_env(s):
    """展开 ${VAR} 与 ${VAR:-默认}（未定义且无默认 → 空串）。

    迭代展开以支持嵌套（如 ${A:-${B:-x}/y}）：re.sub 单次替换不会再扫描
    替换结果，外层展开后残留的内层变量需再跑一轮。
    """

    def _sub(m):
        name, default = m.group(1), m.group(2)
        val = os.environ.get(name)
        if val:
            return val
        return default if default is not None else ""

    for _ in range(10):
        new = _ENV_PAT.sub(_sub, s)
        if new == s:
            break
        s = new
    return s


def _subst(s, base, results, repo, ddir=None):
    """先展开环境变量, 再替换占位符。"""
    s = _expand_env(s)
    if "{repo}" in s:
        s = s.replace("{repo}", repo)
    if "{base}" in s:
        s = s.replace("{base}", base)
    if "{results}" in s:
        s = s.replace("{results}", results)
    if ddir is not None and "{dir}" in s:
        s = s.replace("{dir}", ddir)
    return s


def _resolve(obj, base, results, repo, ddir=None):
    """递归替换字符串中的占位符与环境变量。"""
    if isinstance(obj, str):
        return _subst(obj, base, results, repo, ddir)
    if isinstance(obj, dict):
        return {k: _resolve(v, base, results, repo, ddir) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_resolve(v, base, results, repo, ddir) for v in obj]
    return obj


def load_config(path=None):
    """读取并缓存配置。path 缺省用 pipeline 根目录 datasets.yaml。

    顶层字段先解析（base/results 自身支持环境变量与 {repo}），
    再对每个 dataset 解析一次（此时可用 {dir}）。
    """
    global _cache
    if _cache is not None:
        return _cache
    if yaml is None:
        raise SystemExit("[dataset_config] 需要 pyyaml，请安装: pip install pyyaml")
    p = path or _CONFIG_PATH
    if not os.path.exists(p):
        raise SystemExit(f"[dataset_config] 配置文件不存在: {p}")
    with open(p, encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}

    # 顶层: base / results 支持 ${VAR} 与 {repo}
    base = _subst(str(raw.get("base", "")), "", "", _REPO_ROOT)
    results = _subst(str(raw.get("results", "")), base, "", _REPO_ROOT)

    cfg = dict(raw)
    cfg["base"] = base
    cfg["results"] = results
    cfg["repo"] = _REPO_ROOT
    cfg["env"] = _resolve(raw.get("env", {}), base, results, _REPO_ROOT)

    datasets = {}
    for name, entry in (raw.get("datasets") or {}).items():
        if not isinstance(entry, dict):
            continue
        # 先解析 dir, 再用它解析该条目其余字段（支持 {dir} 自引用）
        ddir = _subst(str(entry.get("dir", "")), base, results, _REPO_ROOT)
        datasets[name] = _resolve(entry, base, results, _REPO_ROOT, ddir)
    cfg["datasets"] = datasets

    _cache = cfg
    return cfg


def dataset(name, key=None, path=None):
    """按病毒名取数据集条目（自动大写匹配）。key 为 None 返回整个条目。"""
    d = load_config(path).get("datasets", {}).get(str(name).upper())
    if d is None:
        return None
    return d.get(key) if key else d


def known_dataset_names(path=None):
    """所有已配置的数据集名（大写）。"""
    return sorted(load_config(path).get("datasets", {}).keys())


def colors(kind="host", path=None):
    """着色方案表，kind ∈ {host, location}。"""
    return load_config(path).get("colors", {}).get(kind, {})


def provinces(path=None):
    """采样点经纬度表 {省: [lat, lon]}。"""
    return load_config(path).get("province", {})


def env(key=None, path=None):
    """运行环境配置（beast_cp/beast_jlp/pypopart_path）。"""
    e = load_config(path).get("env", {})
    return e.get(key) if key else e


def base(path=None):
    return load_config(path).get("base", "")


def results(path=None):
    """运行产物根目录（{results} 占位符的值）。"""
    return load_config(path).get("results", "")


def repo(path=None):
    """仓库根目录（{repo} 占位符的值）。"""
    return load_config(path).get("repo", _REPO_ROOT)
