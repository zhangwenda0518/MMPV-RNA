"""MMPV 跨管线共享模块。

当前内容:
  io_layout     -- 五管线统一 I/O 目录布局 (legacy / standard 双布局)
  report_style  -- 七（八）管线 HTML 报告共享样式框架
  ai            -- AI 基础设施 (L2): provider 抽象/版本化提示词/溯源日志/离线降级
  stage_cards   -- 阶段卡片单一真相源 (L1 知识层)
  run_digest    -- 项目运行状态数字化 (L1 状态层)
  advisor       -- 规则顾问 + AI 叙事 (L3)
  triage        -- 候选三级标注, 只加 ai_* 列 (L3)
  copilot       -- 自然语言 → 白名单校验命令草稿, 不执行 (L4)
  mcp_server    -- MCP server: 把上述能力开放给 AI 客户端 (L4)
  process_report -- 过程报告层: digest 时间线 + AI 参与清单 + agent 留痕聚合
                    → 确定性过程面板 (嵌任意报告) + 项目级主报告 (L2)
  report_template -- 报告层 Jinja2 模板统一出口 (templates/*.html.j2,
                    autoescape + StrictUndefined, 零外部资源)
  provenance_export -- AI 参与溯源一键导出 (JSON + 方法学 MD, 零 AI 调用)
  figure_export -- 图件四格式统一导出 (SVG/PNG/PDF/PNG@1200dpi + savefig 钩子)

仓库根 mmpv_ai.py 为统一 CLI 入口; 设计论证见 doc/AI_DRIVEN_DESIGN.md。
图件导出规范见 doc/FIGURE_EXPORT.md。
"""

__version__ = "0.2.1"
