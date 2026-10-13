#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ai.py — MMPV AI 基础设施 (四层递进架构 · L2 解读层核心)
========================================================

全平台唯一的 AI 出入口。设计纪律沿袭 gsa_sra.info.py 的元数据清洗经验:

  1. AI 输出一律是「建议」: 调用方必须落 AI 专属字段 (ai_* 前缀),
     永不覆盖任何数据列 / 任何统计结论 (triage.py 强制执行)。
  2. 全量溯源: 每次调用把 prompt 元数据 (模板名/版本/渲染体哈希)、
     模型、温度、token 用量、耗时、输入摘要哈希写进
     <项目>/90_Handoff/state/ai_log.jsonl —— 论文可复现性的凭证。
  3. 离线优雅降级: 未配 key / 无网络 / SDK 缺失时, 所有高层接口返回
     status="unavailable" 的结构化结果, 绝不抛异常打断管线。
  4. 提示词版本化: 模板登记在 PROMPTS, 渲染体带 sha256 ——
     同一版本同输入必得同 prompt, 跨项目可比。

配置 (环境变量, 优先级从高到低):
  MMPV_AI_API_KEY   > DEEPSEEK_API_KEY        API key
  MMPV_AI_BASE_URL  (默认 https://api.deepseek.com)   OpenAI 兼容端点
                    (Ollama: http://localhost:11434/v1, key 随便填 "ollama")
  MMPV_AI_MODEL     (默认 deepseek-chat)
  MMPV_AI_TIMEOUT   (默认 120 秒)

纯 OpenAI 兼容协议 (openai SDK 可选依赖): 无 SDK / 无 key 时 ai_status()
返回 (False, 原因), 调用方走规则回退。
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import re
import time
from dataclasses import dataclass, field
from html import escape as _html_escape
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# ═══════════════════════════════════════════════════════════════════
# 配置
# ═══════════════════════════════════════════════════════════════════

DEFAULT_BASE_URL = "https://api.deepseek.com"
DEFAULT_MODEL = "deepseek-chat"
DEFAULT_TIMEOUT = 120.0


@dataclass
class AIConfig:
    """一次 AI 会话的连接配置。从环境变量构建: from_env()。"""

    api_key: str = ""
    base_url: str = DEFAULT_BASE_URL
    model: str = DEFAULT_MODEL
    timeout: float = DEFAULT_TIMEOUT
    temperature: float = 0.2   # 解读类任务要稳, 默认低温; triage 建议 0.0

    @classmethod
    def from_env(cls, env: Optional[Dict[str, str]] = None) -> "AIConfig":
        e = env if env is not None else os.environ
        key = e.get("MMPV_AI_API_KEY") or e.get("DEEPSEEK_API_KEY") or ""
        base = e.get("MMPV_AI_BASE_URL") or DEFAULT_BASE_URL
        model = e.get("MMPV_AI_MODEL") or DEFAULT_MODEL
        try:
            timeout = float(e.get("MMPV_AI_TIMEOUT") or DEFAULT_TIMEOUT)
        except ValueError:
            timeout = DEFAULT_TIMEOUT
        return cls(api_key=key, base_url=base, model=model, timeout=timeout)

    def masked_key(self) -> str:
        if not self.api_key:
            return ""
        if len(self.api_key) <= 8:
            return "***"
        return self.api_key[:4] + "***" + self.api_key[-4:]


class AIUnavailable(RuntimeError):
    """AI 不可用 (未配 key / SDK 缺失 / 网络失败)。携带人类可读原因。"""


def ai_status(config: Optional[AIConfig] = None) -> Tuple[bool, str]:
    """返回 (是否可用, 原因说明)。不发起网络请求, 只查配置与 SDK。"""
    cfg = config or AIConfig.from_env()
    if not cfg.api_key:
        return False, ("API key not configured: export MMPV_AI_API_KEY=... "
                       "(or DEEPSEEK_API_KEY, same convention as gsa_sra.info.py)")
    try:
        import openai  # noqa: F401
    except ImportError:
        return False, "openai SDK missing (pip install openai)"
    return True, f"{cfg.model} @ {cfg.base_url} (key {cfg.masked_key()})"


# ═══════════════════════════════════════════════════════════════════
# 版本化提示词注册表
# ═══════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class PromptMeta:
    name: str
    version: str
    sha256: str          # 渲染后全文哈希 (可复现凭证)
    text: str = ""       # 仅内存传递, 不写入日志正文之外的地方

    def to_dict(self) -> Dict[str, str]:
        return {"name": self.name, "version": self.version, "sha256": self.sha256}


PROMPTS: Dict[str, Dict[str, str]] = {}


def register_prompt(name: str, version: str, template: str) -> None:
    """登记提示词模板。同 (name, version) 重复登记以内容一致为条件, 否则报错。"""
    old = PROMPTS.get(name, {}).get(version)
    if old is not None and old != template:
        raise ValueError(f"提示词 {name}@{version} 已存在且内容不同, 请升版本号")
    PROMPTS.setdefault(name, {})[version] = template


def render_prompt(name: str, version: Optional[str] = None, *,
                  strict: bool = True, **fields: Any) -> Tuple[str, PromptMeta]:
    """渲染提示词, 返回 (渲染文本, 元数据)。默认取该名字的最高版本。"""
    if name not in PROMPTS:
        raise KeyError(f"未登记的提示词: {name}")
    if version is None:
        version = sorted(PROMPTS[name])[-1]
    if version not in PROMPTS[name]:
        raise KeyError(f"提示词 {name} 无版本 {version} (现有: {sorted(PROMPTS[name])})")
    text = PROMPTS[name][version]
    if strict:
        text = text.format(**fields)
    else:
        text = text.format_map(_SafeDict(fields))
    sha = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return text, PromptMeta(name=name, version=version, sha256=sha, text=text)


class _SafeDict(dict):
    """format_map 宽松模式: 缺字段留空占位而不抛错。"""

    def __missing__(self, key):  # noqa: D105
        return f"<{key}?>"


# ═══════════════════════════════════════════════════════════════════
# 调用层
# ═══════════════════════════════════════════════════════════════════

@dataclass
class AIResult:
    status: str                        # ok | unavailable | error
    text: str = ""
    model: str = ""
    prompt: Optional[PromptMeta] = None
    usage: Dict[str, Any] = field(default_factory=dict)
    latency_s: float = 0.0
    reason: str = ""                   # status != ok 时的人类可读原因
    cached: bool = False               # 命中本地缓存 (同 prompt 同模型不重复计费)

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {"status": self.status, "model": self.model,
                             "latency_s": round(self.latency_s, 2),
                             "reason": self.reason, "cached": self.cached}
        if self.prompt:
            d["prompt"] = self.prompt.to_dict()
        if self.usage:
            d["usage"] = self.usage
        return d


def _build_client(cfg: AIConfig):
    try:
        from openai import OpenAI
    except ImportError as e:
        raise AIUnavailable(f"openai SDK missing: {e}") from e
    try:
        return OpenAI(api_key=cfg.api_key, base_url=cfg.base_url, timeout=cfg.timeout)
    except Exception as e:  # pragma: no cover - SDK 初始化失败视环境而定
        raise AIUnavailable(f"OpenAI client init failed: {e}") from e


# ── 本地缓存 (报告重渲染不重复计费; MMPV_AI_CACHE=0 关闭) ────────────

def _ai_cache_file(state_dir: Optional[Path]) -> Optional[Path]:
    if state_dir is None or os.environ.get("MMPV_AI_CACHE", "1") == "0":
        return None
    return Path(state_dir) / "ai_cache.json"


def _cache_key(cfg: AIConfig, messages: List[Dict[str, str]],
               temperature: float, json_mode: bool) -> str:
    payload = {"model": cfg.model, "base_url": cfg.base_url,
               "temperature": temperature, "json_mode": json_mode,
               "messages": messages,
               # key 盐: 只入 8 位哈希 (不落明文) —— 换 key/换账号不命中旧缓存
               "key_salt": _sha(cfg.api_key)[:8] if cfg.api_key else ""}
    return hashlib.sha256(
        json.dumps(payload, ensure_ascii=False, sort_keys=True).encode("utf-8")
    ).hexdigest()


def _cache_get(path: Path, key: str) -> Optional[Dict[str, Any]]:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        hit = data.get(key)
        return hit if isinstance(hit, dict) else None
    except (OSError, json.JSONDecodeError):
        return None


def _cache_put(path: Path, key: str, res: "AIResult") -> None:
    try:
        data = {}
        p = Path(path)
        if p.exists():
            data = json.loads(p.read_text(encoding="utf-8"))
        data[key] = {"text": res.text, "model": res.model, "usage": res.usage,
                     "cached_at": time.strftime("%Y-%m-%d %H:%M:%S")}
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
    except (OSError, json.JSONDecodeError, TypeError):
        pass  # 缓存失败绝不影响主流程


def chat(messages: List[Dict[str, str]], *,
         config: Optional[AIConfig] = None,
         json_mode: bool = False,
         temperature: Optional[float] = None,
         client: Any = None,
         cache_path: Optional[Path] = None) -> AIResult:
    """一次 chat completion。client 可注入 (测试桩)。失败一律 AIResult(error)。

    cache_path 给定时: 同 (模型, 端点, 温度, json_mode, messages) 命中本地
    缓存直接返回 (res.cached=True, latency=0) —— 报告重渲染零计费,
    且同一输入跨次渲染文案稳定。
    """
    cfg = config or AIConfig.from_env()
    ok, why = ai_status(cfg)
    if not ok and client is None:
        return AIResult(status="unavailable", reason=why)
    temp = cfg.temperature if temperature is None else temperature
    if cache_path is not None:
        hit = _cache_get(Path(cache_path), _cache_key(cfg, messages, temp, json_mode))
        if hit is not None:
            return AIResult(status="ok", text=hit.get("text", ""),
                            model=hit.get("model", cfg.model),
                            usage=dict(hit.get("usage") or {}),
                            latency_s=0.0, reason="cache hit", cached=True)
    started = time.time()
    try:
        cli = client if client is not None else _build_client(cfg)
        kwargs: Dict[str, Any] = dict(
            model=cfg.model,
            messages=messages,
            temperature=temp,
        )
        if json_mode:
            kwargs["response_format"] = {"type": "json_object"}
        resp = cli.chat.completions.create(**kwargs)
        text = resp.choices[0].message.content or ""
        usage = {}
        u = getattr(resp, "usage", None)
        if u is not None:
            usage = {"prompt_tokens": getattr(u, "prompt_tokens", None),
                     "completion_tokens": getattr(u, "completion_tokens", None),
                     "total_tokens": getattr(u, "total_tokens", None)}
        res = AIResult(status="ok", text=text, model=cfg.model,
                       usage=usage, latency_s=time.time() - started)
        if cache_path is not None:
            _cache_put(Path(cache_path), _cache_key(cfg, messages, temp, json_mode), res)
        return res
    except AIUnavailable:
        raise
    except Exception as e:
        return AIResult(status="error", reason=f"{type(e).__name__}: {e}",
                        latency_s=time.time() - started)


def _strip_code_fence(text: str) -> str:
    t = text.strip()
    m = re.match(r"^```[a-zA-Z0-9]*\s*\n(.*?)\n?```\s*$", t, re.S)
    return m.group(1).strip() if m else t


def parse_json(text: str) -> Any:
    """宽容解析模型输出中的 JSON (剥代码围栏, 截取首尾大括号)。失败抛 ValueError。"""
    t = _strip_code_fence(text)
    try:
        return json.loads(t)
    except json.JSONDecodeError:
        pass
    for a, b in (("{", "}"), ("[", "]")):
        i, j = t.find(a), t.rfind(b)
        if i != -1 and j > i:
            try:
                return json.loads(t[i:j + 1])
            except json.JSONDecodeError:
                continue
    raise ValueError(f"AI output is not parseable as JSON: {text[:200]!r}")


# ═══════════════════════════════════════════════════════════════════
# 溯源日志 (90_Handoff/state/ai_log.jsonl)
# ═══════════════════════════════════════════════════════════════════

def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def log_ai_event(state_dir: Optional[Path], record: Dict[str, Any]) -> Optional[Path]:
    """追加一条 JSONL 溯源记录。state_dir 为 None 时跳过 (纯内存调用)。"""
    if state_dir is None:
        return None
    state_dir = Path(state_dir)
    state_dir.mkdir(parents=True, exist_ok=True)
    path = state_dir / "ai_log.jsonl"
    body = {"ts": time.strftime("%Y-%m-%d %H:%M:%S"), **record}
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(body, ensure_ascii=False) + "\n")
    return path


# ═══════════════════════════════════════════════════════════════════
# 内置提示词 (版本化; 修改措辞必须升版本号)
# ═══════════════════════════════════════════════════════════════════

register_prompt("run_summary", "1.0.0", """你是 MMPV (大规模植物病毒挖掘平台) 的分析解读助手。
下面是一次运行的 run digest (机器提取的阶段状态与关键计数 JSON)。

【纪律】你只解读, 不捏造: 数字必须出自 digest, digest 里没有的写"digest 未提供"。
不要给出超出 digest 证据范围的确凿结论; 不确定就用"建议核实"。

run digest:
{digest_json}

请用英文输出解读, 依次三段 (用 markdown 小标题):
1. 运行状态概览 (哪些管线已跑/未跑, 用 digest 数字)
2. 关键结果解读 (病毒发现/分类/拯救等关键计数说明了什么)
3. 建议的下一步 (与 advisor 规则建议相衔接, 最多 3 条)
""")

register_prompt("advisor_narrative", "1.0.0", """你是 MMPV 的分析规划助手。规则引擎已产出如下建议列表 (JSON),
以及本次运行的 run digest。

规则建议:
{rules_json}

run digest 摘要:
{digest_json}

请用英文写一段不超过 200 字的总体叙事: 把规则建议串成一条有优先顺序的执行路径,
点明为什么。不要新增规则里没有的建议; 数字只引用 digest 与规则建议里出现过的。
""")

register_prompt("row_triage", "1.0.0", """你是 MMPV 的候选 triage 助理。给你一批记录 (JSON 数组),
每条含 id 与若干证据列。请对每条给出:
  action: "KEEP" (证据一致支持保留) 或 "REVIEW" (证据冲突/信息不足)
  reason: 不超过 40 字的英文理由, 必须引用该条的证据值

【纪律】你只做标注, 不是裁判: 拿不准一律 REVIEW; 绝不输出删除类动作。
只输出 JSON 数组: [{{"id": "...", "action": "...", "reason": "..."}}]

记录:
{rows_json}
""")

register_prompt("process_narrative", "1.0.0", """你是 MMPV 平台的"过程叙事"助手。下面是一个项目的完整过程记录 JSON: 各管线推进状态、时间线快照、AI 参与清单 (哪些环节调过 AI、做了什么)、agent 受控运行记录。

【纪律】数字与事实只准出自给定 JSON; JSON 里没有 AI 调用记录就如实写"无 AI 调用记录 (纯规则运行)"; 不评价未在记录中的内容。

过程记录:
{process_json}

请用英文写不超过 300 字的过程叙事, 三个自然段:
1. 项目推进程度 (主链 ①→⑧ 走到哪, 关键计数)
2. 过程中 AI 参与了什么 (解读/建议/标注/规划各多少次, 在哪些环节)
3. 从过程看, 当前最值得做的下一步 (与建议清单衔接, 最多 2 条)
只输出正文, 不要标题。
""")

register_prompt("chat_intent", "1.1.0", """你是 MMPV 平台的命令规划助理。用户想用自然语言描述分析意图,
你需要把它映射到 MMPV 已有编排器的一条命令草稿。

可用管线 (key: 说明):
{pipelines_text}

各管线的合法 stage (标"子命令"的编排器没有 --stage, 只能选子命令):
{stages_text}

【纪律】只能从上面的管线/stage/子命令里选; 缺参数就如实列出 missing, 不要编造路径。
只输出 JSON:
{{"pipeline": "...", "subcommand": "无子命令时留空字符串", "stages": ["..."],
 "args": {{"--flag": "value"}}, "missing": ["缺什么参数的说明"],
 "explain": "一句话英文说明"}}
""")

register_prompt("ai_interpretation", "1.0.0", """你是 MMPV 的报告解读助手。下面是管线阶段统计 (JSON)。

{stage_stats_json}

请输出不超过 250 字的英文解读: 运行到哪一步、关键数字说明什么、下一步建议。
数字必须出自给定 JSON, 没有的信息写"报告未提供"。
""")


# ═══════════════════════════════════════════════════════════════════
# 高层封装 (供报告/CLI/MCP 直接使用, 全部离线安全)
# ═══════════════════════════════════════════════════════════════════

def summarize_run(digest: Dict[str, Any], *,
                  state_dir: Optional[Path] = None,
                  config: Optional[AIConfig] = None,
                  client: Any = None) -> Dict[str, Any]:
    """run digest → 中文运行解读。返回 {status, summary, meta?}。"""
    prompt_text, meta = render_prompt("run_summary",
                                      digest_json=json.dumps(digest, ensure_ascii=False)[:30000])
    res = chat([{"role": "user", "content": prompt_text}],
               config=config, client=client, cache_path=_ai_cache_file(state_dir))
    if res.status != "ok":
        return {"status": res.status, "summary": "", "reason": res.reason}
    log_ai_event(state_dir, {
        "kind": "call", "op": "summarize_run", "model": res.model, "cached": res.cached,
        "prompt": meta.to_dict(), "usage": res.usage,
        "latency_s": round(res.latency_s, 2),
        "digest_sha": _sha(json.dumps(digest, ensure_ascii=False)),
    })
    return {"status": "ok", "summary": res.text.strip(),
            "meta": {**res.to_dict(), "prompt": meta.to_dict()}}


def ai_interpretation_section(stage_stats: Dict[str, Any], *,
                              state_dir: Optional[Path] = None,
                              config: Optional[AIConfig] = None,
                              client: Any = None) -> Dict[str, Any]:
    """报告生成器接线点: stage_stats → {html, status}。

    status != ok 时返回空 html, 报告端原样跳过即可 —— 离线渲染不受影响。
    """
    prompt_text, meta = render_prompt("ai_interpretation",
                                      stage_stats_json=json.dumps(stage_stats, ensure_ascii=False)[:20000])
    res = chat([{"role": "user", "content": prompt_text}], config=config,
               client=client, cache_path=_ai_cache_file(state_dir))
    if res.status != "ok":
        return {"status": res.status, "html": "", "reason": res.reason}
    log_ai_event(state_dir, {
        "kind": "call", "op": "ai_interpretation", "model": res.model, "cached": res.cached,
        "prompt": meta.to_dict(), "usage": res.usage,
        "latency_s": round(res.latency_s, 2),
    })
    from mmpv_common.report_style import PALETTE
    # LLM 输出是不可信文本 → 一律 html 转义后插入 (与 note 参数同做法),
    # 防止模型输出 <script>/标签破坏报告标记
    body_html = _html_escape(res.text.strip())
    html = (
        f'<section id="mmpv-ai-section" style="margin:18px 0;padding:14px 18px;'
        f'background:{PALETTE["primary-light"]};border:1px solid {PALETTE["border"]};'
        f'border-radius:8px">'
        f'<h2 style="color:{PALETTE["primary"]};font-size:16px;margin-bottom:8px">'
        f'AI Interpretation <span style="font-size:11px;color:{PALETTE["ink-secondary"]}">'
        f'{res.model} · prompt {meta.name}@{meta.version}</span></h2>'
        f'<div style="font-size:13px;line-height:1.7">{body_html}</div></section>'
    )
    return {"status": "ok", "html": html,
            "meta": {**res.to_dict(), "prompt": meta.to_dict()}}


def ask(question: str, context: str = "", *,
        state_dir: Optional[Path] = None,
        config: Optional[AIConfig] = None,
        client: Any = None) -> Dict[str, Any]:
    """自由问答 (context 通常为 compact digest 或卡片摘要)。离线返回 unavailable。"""
    content = question if not context else (
        f"背景资料 (只引用其中出现的事实):\n{context[:25000]}\n\n问题: {question}")
    res = chat([{"role": "user", "content": content}], config=config,
               client=client, cache_path=_ai_cache_file(state_dir))
    if res.status != "ok":
        return {"status": res.status, "answer": "", "reason": res.reason}
    log_ai_event(state_dir, {"kind": "call", "op": "ask", "model": res.model, "cached": res.cached,
                             "usage": res.usage, "latency_s": round(res.latency_s, 2),
                             "question_sha": _sha(question)})
    return {"status": "ok", "answer": res.text.strip(),
            "meta": res.to_dict()}


AI_SECTION_ID = "mmpv-ai-section"


def append_ai_section_to_html(html_path: Path, stage_stats: Dict[str, Any], *,
                              state_dir: Optional[Path] = None,
                              config: Optional[AIConfig] = None,
                              client: Any = None,
                              note: str = "") -> Dict[str, Any]:
    """五报告统一 AI 接线点: 在已落盘 HTML 的 </body> 前插入 AI 解读段。

    三条契约:
      1. 幂等: 已含 AI 段 (id=mmpv-ai-section) 直接返回 already, 不调 AI 不改文件
         —— 事后补解读 (mmpv_ai.py enrich) 可反复跑;
      2. 纯增量: 离线/失败 (status != ok) 时【不改动文件】;
      3. 先查文件再调 AI: 无 body/已插入的文件不浪费一次调用。
    note: 追加在段内的小字来源说明 (如事后补解读的数据来源)。
    """
    html_path = Path(html_path)
    try:
        html = html_path.read_text(encoding="utf-8")
    except OSError as e:
        return {"status": "error", "written": False, "reason": f"读不回 HTML: {e}"}
    if AI_SECTION_ID in html:
        return {"status": "already", "written": False, "reason": "AI 解读段已存在"}
    if "</body>" not in html:
        return {"status": "error", "written": False,
                "reason": "HTML 无 </body>, 放弃插入 (不改动文件)"}
    res = ai_interpretation_section(stage_stats, state_dir=state_dir,
                                    config=config, client=client)
    res.setdefault("written", False)
    if res["status"] != "ok":
        return res
    section = res["html"]
    if note:
        from mmpv_common.report_style import PALETTE
        _note = (str(note).replace("&", "&amp;").replace("<", "&lt;")
                 .replace(">", "&gt;"))
        badge = (f'<div style="font-size:11px;color:{PALETTE["ink-secondary"]};'
                 f'margin-top:8px">{_note}</div>')
        section = section.replace("</section>", badge + "</section>", 1)
    html = html.replace("</body>", section + "\n</body>", 1)
    try:
        html_path.write_text(html, encoding="utf-8")
    except OSError as e:
        res["reason"] = f"写回失败: {e}"
        res["status"] = "error"
        return res
    res["written"] = True
    return res
