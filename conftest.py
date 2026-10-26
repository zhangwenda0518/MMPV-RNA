# -*- coding: utf-8 -*-
"""pytest 全局隔离: 测试期间关闭图件四格式自动导出。

背景
----
2026-10-20 起所有出图脚本在导入时安装 `Figure.savefig` 钩子 (默认产出
SVG/PNG/PDF/hires 四格式)。这在生产环境是期望行为, 但会污染测试:

  · `test_plot_io.py` 的 `pdf_to_png` 用例先造 PDF 再断言「兜底渲染出 PNG」;
    钩子已在造 PDF 时顺手写了 PNG → 兜底逻辑判定「已存在」而跳过 → 假失败。
  · 全量 pytest 时若先导入了任一管线脚本, 钩子就被装上; 单跑该文件则不会,
    于是表现为「单跑绿、全量红」。

口径: 测试默认 `MMPV_FIG_EXPORT=0` (关), 需要验证四格式的用例自己显式打开
(见 `mmpv_common/tests/test_figure_export.py` 用 monkeypatch 控制环境变量)。
"""
import os
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

# 必须在任何管线模块被导入前设置 —— conftest 在收集阶段最早执行。
os.environ["MMPV_FIG_EXPORT"] = "0"


def pytest_configure(config):
    """双保险: 若有模块已在 conftest 之前装上钩子, 这里撤掉。"""
    try:
        from mmpv_common.figure_export import uninstall_savefig_hook
        uninstall_savefig_hook()
    except Exception:  # noqa: BLE001 — 测试隔离失败不应中断收集
        pass


# ── 同名顶层包 utils/ 隔离 ────────────────────────────────────────────────────
"""一条 pytest 命令跑多个管线目录时, 各自的 `utils/` 顶层包会互相顶掉。

问题
----
仓库有 4 个管线各自带一个**同名顶层包** `utils/`:

  · public_metadata_pipeline/utils/   (pipeline_utils.py …)
  · virome_phylo_pipeline/utils/      (ess.py, plot_io.py …)
  · virome_discovery_pipeline/utils/  (namespace 包, 无 __init__.py)
  · virome_analysis_pipeline/utils/   (namespace 包, 无 __init__.py)

测试代码写的是裸 `from utils.X import ...`, 靠把自己管线目录插进 sys.path 生效。
单目录跑没问题; 多目录一条命令时, 谁先把自己 import 进 `sys.modules["utils"]`,
后面所有目录就只能看到那一个 `utils`:

  收集 virome_phylo_pipeline/tests → `sys.modules["utils"]` = phylo/utils
  收集 public_metadata_pipeline/tests → `import utils.pipeline_utils`
      → 只在 phylo/utils 里找 → ModuleNotFoundError: No module named 'utils.pipeline_utils'

为什么要放在**根 conftest + 收集期钩子**, 而不是各管线的 tests/conftest.py
------------------------------------------------------------------------
直觉方案是「各目录 tests/conftest.py 里插 sys.path + 清 utils 缓存」, 理由是
「dir2 的 conftest 总在 dir2 的测试模块之前 import」。这句话本身没错, 但**不够**:
pytest 8 在 `pytest_load_initial_conftests` 阶段就把**所有命令行参数目录**的
conftest.py 一次性全部提前 import, 之后才开始收集。所以 ⑥ 的 conftest 清缓存发生
在 ⑤ 的测试模块 import `utils` **之前** —— 清了也白清: 随后收集 ⑤ 目录时 phylo 的
utils 又被装回 sys.modules, 等轮到 ⑥ 时早就不干净了。

实测 (2026-10-26, pytest 8.4.2) 在 public_metadata_pipeline/tests/conftest.py 里
`traceback.print_stack()` 得到的调用链:

    _pytest/config/__init__.py:1431  pytest_load_initial_conftests
    _pytest/config/__init__.py:583   _set_initial_conftests
    _pytest/config/__init__.py:621   _try_load_conftest   ← 所有 arg 的 conftest 都在这里
    _pytest/config/__init__.py:661   _loadconftestmodules
    _pytest/config/__init__.py:712   _importconftest

即目录级 conftest 的副作用早于**任何收集**, 因此同步不到「收集到某个目录」这一刻。
能对上这一刻的只有按节点触发的钩子:

  · `pytest_collectstart`  在 `Module.collect()` 之前触发, 正好在本目录测试模块
    被 import 之前 —— 与单目录跑时的时机一致;
  · `pytest_runtest_setup`  在每个用例的 fixture 之前触发, 覆盖函数体内的
    惰性 `from utils.X import ...` (例如 test_plot_io.py 里测试函数内部的导入)。

做法
----
按节点路径找到它所属的「带 utils/ 的管线目录」, 把该目录顶到 sys.path[0], 并把
`sys.modules` 里的 `utils*` 换成该管线的副本。切换时**不重新 import**, 而是从
收集期存下的快照里把原来那批模块对象原样装回去 —— 这样每个管线的 utils 全跑完
只被 import 一次, 模块 identity 与分开跑完全一致 (避免二次 import 的副作用,
也让 test_plot_io.py 那种靠 `__globals__` 打补丁的用例继续有效)。

不改任何被测源码、不改测试断言/用例数, 只让 `utils` 的解析与模块缓存按目录隔离。
"""

_UTILS_SNAPSHOTS = {}   # 管线目录 -> {模块名: 模块对象}  收集期留存, 运行期原样还原
_UTILS_PIPELINE_CACHE = {}  # 节点路径 -> 所属管线目录 (None 表示不属于任何管线)


def _utils_pipeline_for(path):
    """节点路径所属的、带同名顶层包 `utils/` 的管线目录; 没有则 None。"""
    if path is None:
        return None
    cache_key = str(path)
    if cache_key in _UTILS_PIPELINE_CACHE:
        return _UTILS_PIPELINE_CACHE[cache_key]
    result = None
    try:
        here = Path(path).resolve()
    except OSError:  # 已删除/不可访问的路径: 当作无管线
        here = None
    if here is not None:
        if here.is_file():
            here = here.parent
        for parent in [here, *here.parents]:
            if (parent / "utils").is_dir():
                result = parent
                break
    _UTILS_PIPELINE_CACHE[cache_key] = result
    return result


def _utils_owner():
    """当前 `sys.modules["utils"]` 绑定到哪个管线; 未加载 -> None, 跨管线 -> "<mixed>"。"""
    mod = sys.modules.get("utils")
    if mod is None:
        return None
    locations = []
    origin = getattr(mod, "__file__", None)
    if origin:
        locations.append(Path(origin))
    else:  # namespace 包 (discovery/analysis 的 utils 没有 __init__.py)
        locations.extend(Path(p) for p in (getattr(mod, "__path__", None) or []))
    owners = set()
    for loc in locations:
        utils_dir = loc.parent if loc.is_file() else loc
        if utils_dir.name != "utils":
            return None           # 不是本仓库管线的 utils, 不归我们管
        owners.add(str(utils_dir.parent))
    if not owners:
        return None
    return owners.pop() if len(owners) == 1 else "<mixed>"


def _snapshot_utils(pipeline):
    """把当前 utils* 模块对象按管线留存 (只在绑定的确是这条管线时调用)。"""
    _UTILS_SNAPSHOTS[pipeline] = {
        name: sys.modules[name]
        for name in list(sys.modules)
        if name == "utils" or name.startswith("utils.")
    }


def _hoist_utils_dir(pipeline):
    if pipeline in sys.path:
        sys.path.remove(pipeline)
    sys.path.insert(0, pipeline)


def _activate_utils_for(path):
    """让本节点的 `utils` 解析到它自己那条管线的 utils/ (幂等, 未切换时只调整 path 顺序)。"""
    pipeline = _utils_pipeline_for(path)
    if pipeline is None:
        return
    pipeline = str(pipeline)
    owner = _utils_owner()
    if owner == pipeline:
        _snapshot_utils(pipeline)   # 顺手刷新快照 (本管线又新 import 了子模块)
        _hoist_utils_dir(pipeline)
        return
    if owner is not None and owner != "<mixed>":
        _snapshot_utils(owner)      # 换人之前先把上一条管线的模块对象存下来
    for name in [n for n in list(sys.modules) if n == "utils" or n.startswith("utils.")]:
        del sys.modules[name]
    _hoist_utils_dir(pipeline)
    saved = _UTILS_SNAPSHOTS.get(pipeline)
    if saved:                       # 二次进入同一管线: 原样还原, 不再 import 一遍
        sys.modules.update(saved)


def pytest_collectstart(collector):
    """在收集每个节点 (含目录/测试模块) 之前切好 utils —— 早于测试模块的 import。"""
    _activate_utils_for(getattr(collector, "path", None))


def pytest_runtest_setup(item):
    """每个用例开工前再切一次, 覆盖测试函数体内的惰性 `from utils.X import ...`。"""
    _activate_utils_for(getattr(item, "path", getattr(item, "fspath", None)))
