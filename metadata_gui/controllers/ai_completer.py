"""AI metadata completion — wraps gsa_sra.info.py sanitizer prompts for GUI use."""

import json
import re
import os

from PySide6.QtCore import QThread, Signal

# Same prompt as gsa_sra.info.py PROMPT_DATA_SANITIZER
AI_SYSTEM_PROMPT = """你是一个极其严谨的生物信息元数据清理程序。你的任务是【去污染、归位、格式化】，绝对不是【凭空捏造】。
我将给你一条 JSON 数据。请严格按以下规则逐字段处理：

1. **全面净化**：遇到 'not collected', 'missing', 'N/A', 'Not_Provided', 空字符串等无意义占位符，替换为 'Not_Provided'。

2. **ScientificName**：如果混入了部位（如 leaf）或来源词（如 wild），将其剔除，只保留纯物种名。

3. **Tissue**：强制小写单数（leaves→leaf, fruits→fruit, flowers→flower, roots→root）。如果原值为空但 Source 恰好是真正的部位词（如 leaf, root），才移动过来。Source 是品种名（如 Ningqi No.1）则不移动。

4. **Source**：如果被填成了物种名，清空为 Not_Provided。如果是 "wild", "cultivated" 或品种名（如 Ningqi No.1），必须 100% 保持原样。

5. **Location**：根据 CenterName 推断机构地理位置，输出严格三级格式【国家, 省/州, 市/县_AI】。参照示例：
   "Beijing Forestry University" → "China, Beijing, Beijing_AI"
   "Ningxia University" → "China, Ningxia, Yinchuan_AI"
   "North Minzu University" → "China, Ningxia, Yinchuan_AI"
   "University of Tokyo" → "Japan, Tokyo, Tokyo_AI"
   "USDA-ARS" → "USA, Maryland, Beltsville_AI"
   "Royal Botanic Gardens Kew" → "United Kingdom, England, London_AI"
   "Northwest A&F University" → "China, Shaanxi, Yangling_AI"
   "Gansu Agricultural University" → "China, Gansu, Lanzhou_AI"
   "Xinjiang University" → "China, Xinjiang, Urumqi_AI"
   "Qinghai University" → "China, Qinghai, Xining_AI"
   "Inner Mongolia University" → "China, Inner Mongolia, Hohhot_AI"
   "Henan Agricultural University" → "China, Henan, Zhengzhou_AI"
   不认识的机构保持 Not_Provided，绝对不要输出 Unknown。

6. **Age_GrowthStage**：标准化年龄和发育阶段文本。规则：
   - "3 year" → "3 years"（补全复数）
   - "2 years" → 保持原样（已规范）
   - "Archeocyte stage" → "archeocyte stage"（统一为首字母不大写的全小写，除非是专有名词）
   - 多个阶段用 "|" 分隔，如 "3 years | mature fruit stage"
   - 将中文风格描述转为英文，如 "flowering stage"、"ripening stage"
   - 如果为空则保持 Not_Provided

7. **LibrarySource**：如果传入长文本，浓缩为一个词：TRANSCRIPTOMIC / GENOMIC / METAGENOMIC。

8. **CollectionDate, CenterName, BioProject**：必须 100% 保持输入原样，一个字都不许改！

【核心纪律】：只做减法（去污染）、归位（移动）、格式标准化，绝对禁止凭空捏造！Location 没线索就 Not_Provided！

直接输出 9 个键的合法 JSON，勿带 ```json 标记：
{
  "CollectionDate": "...", "Location": "...", "Source": "...", "Tissue": "...",
  "Age_GrowthStage": "...", "ScientificName": "...", "LibrarySource": "...", "CenterName": "...", "BioProject": "..."
}"""

FILLABLE_COLS = [
    "CollectionDate", "Location", "Source", "Tissue",
    "Age_GrowthStage", "ScientificName", "LibrarySource",
    "CenterName", "BioProject",
]

EMPTY_VALS = {"", " ", "NA", "N/A", "Not_Provided", "not_provided",
              "not collected", "missing", "none", "unknown", "nan"}


class AICompleteWorker(QThread):
    """Background thread: call DeepSeek API to fill missing metadata fields."""
    progress = Signal(int, int)   # current, total
    finished = Signal(int)        # filled count
    error = Signal(str)

    def __init__(self, records: list, api_key: str, api_base: str,
                 model: str = "deepseek-chat"):
        super().__init__()
        self._records = records
        self._api_key = api_key
        self._api_base = api_base
        self._model = model
        self._filled = 0

    def run(self):
        try:
            from openai import OpenAI
            client = OpenAI(api_key=self._api_key, base_url=self._api_base)
        except Exception as e:
            self.error.emit(f"OpenAI client init failed: {e}")
            return

        total = len(self._records)
        results = []

        for i, (row_idx, record) in enumerate(self._records):
            target = {}
            for col in FILLABLE_COLS:
                val = record.get(col, "")
                target[col] = str(val).strip() if val else "Not_Provided"

            try:
                kwargs = self._build_kwargs(
                    self._model, AI_SYSTEM_PROMPT,
                    json.dumps(target, ensure_ascii=False))
                response = client.chat.completions.create(**kwargs)
                text = response.choices[0].message.content or ""
                text = self._clean_json(text)

                result = json.loads(text)
                applied = {}
                for col in FILLABLE_COLS:
                    if col not in result:
                        continue
                    old_val = str(record.get(col, "")).strip()
                    new_val = str(result[col]).strip()
                    if old_val == new_val:
                        continue  # no change
                    old_empty = old_val.lower() in EMPTY_VALS
                    new_empty = new_val.lower() in EMPTY_VALS

                    # Skip if AI output is noise
                    if new_empty:
                        continue
                    # Skip Location that is all Unknown
                    if col == "Location":
                        parts = [p.strip() for p in new_val.split(",")]
                        if all("unknown" in p.lower() for p in parts):
                            continue
                    # Accept: old was empty → fill; old was dirty → clean
                    if old_empty or old_val != new_val:
                        applied[col] = new_val
                        self._filled += 1
                results.append((row_idx, applied))
            except json.JSONDecodeError as e:
                self.error.emit(
                    f"Row {row_idx+1}: AI returned invalid JSON: {str(e)[:200]}")
            except Exception as e:
                msg = str(e)[:300]
                self.error.emit(f"Row {row_idx+1}: {msg}")

            self.progress.emit(i + 1, total)

        # Commit results back to the worker's output
        self._results = results
        self.finished.emit(self._filled)

    @property
    def results(self) -> list:
        return getattr(self, '_results', [])

    def _build_kwargs(self, model, sys_prompt, user_content):
        kwargs = {
            "model": model,
            "messages": [
                {"role": "system", "content": sys_prompt},
                {"role": "user", "content": user_content},
            ],
        }
        if "v4" in model.lower():
            kwargs["stream"] = False
            if "pro" in model.lower():
                kwargs["reasoning_effort"] = "high"
                kwargs["extra_body"] = {"thinking": {"type": "enabled"}}
            else:
                kwargs["temperature"] = 0.1
        elif "reasoner" in model.lower() or "r1" in model.lower():
            pass
        else:
            kwargs["temperature"] = 0.1
        return kwargs

    @staticmethod
    def _clean_json(text: str) -> str:
        text = re.sub(r'<think>.*?</think>', '', text, flags=re.DOTALL).strip()
        return re.sub(r'^```(?:json)?\s*|\s*```$', '', text,
                      flags=re.IGNORECASE)
