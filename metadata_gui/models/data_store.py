"""Central data store wrapping metadata TSV files with pandas."""

import os
from typing import Optional, List, Set
import numpy as np
import pandas as pd

CORE13_COLS = [
    "Run", "ReleaseDate", "CollectionDate", "Location", "Source",
    "Tissue", "Age_GrowthStage", "ScientificName", "TaxID",
    "LibrarySource", "CenterName", "BioProject", "PMID"
]

NA_PLACEHOLDERS = {"NA", "N/A", "Not_Provided", "not collected",
                   "missing", "none", "unknown", "", " "}


class MetadataStore:
    """In-memory store for metadata records with change tracking."""

    def __init__(self):
        self._df: Optional[pd.DataFrame] = None
        self._full_df: Optional[pd.DataFrame] = None
        self._filepath: Optional[str] = None
        self._full_filepath: Optional[str] = None
        self._modified = False
        self._deleted_rows: List[int] = []
        self._ai_filled: set = set()  # (row, col_name) tuples
        self._imported: set = set()    # row indices from online search import
        # Cache
        self._unique_cache: dict = {}
        self._stats_cache: Optional[dict] = None

    def _invalidate_cache(self):
        self._unique_cache.clear()
        self._stats_cache = None

    @property
    def dataframe(self) -> pd.DataFrame:
        if self._df is None:
            self._df = pd.DataFrame(columns=CORE13_COLS)
        return self._df

    @property
    def full_dataframe(self) -> pd.DataFrame:
        return self._full_df if self._full_df is not None else self.dataframe

    @property
    def is_loaded(self) -> bool:
        return self._df is not None and len(self._df) > 0

    @property
    def is_modified(self) -> bool:
        return self._modified

    @property
    def row_count(self) -> int:
        return len(self._df) if self._df is not None else 0

    @property
    def column_count(self) -> int:
        return len(self._df.columns) if self._df is not None else 0

    @property
    def columns(self) -> List[str]:
        return list(self._df.columns) if self._df is not None else CORE13_COLS

    def load(self, filepath: str) -> bool:
        if not os.path.exists(filepath):
            return False
        try:
            df = pd.read_csv(filepath, sep="\t", dtype=str, keep_default_na=False)
            df.fillna("", inplace=True)
            self._df = df.astype(object)
            self._filepath = filepath
            self._modified = False
            self._deleted_rows = []
            self._invalidate_cache()
            return True
        except Exception:
            return False

    def load_full(self, filepath: str) -> bool:
        if not os.path.exists(filepath):
            return False
        try:
            df = pd.read_csv(filepath, sep="\t", dtype=str, keep_default_na=False)
            df.fillna("", inplace=True)
            self._full_df = df.astype(object)
            self._full_filepath = filepath
            return True
        except Exception:
            return False

    def get_cell(self, row: int, col: int) -> str:
        if self._df is None:
            return ""
        return str(self._df.iat[row, col])

    def set_cell(self, row: int, col: int, value: str):
        if self._df is None:
            return
        self._df.iat[row, col] = value
        self._modified = True
        self._invalidate_cache()

    def get_row(self, row: int) -> dict:
        if self._df is None:
            return {}
        return self._df.iloc[row].to_dict()

    def set_row(self, row: int, data: dict):
        if self._df is None:
            return
        for col, val in data.items():
            if col in self._df.columns:
                self._df.iat[row, self._df.columns.get_loc(col)] = val
        self._modified = True
        self._invalidate_cache()

    def add_row(self, data: dict = None):
        if self._df is None:
            return
        new = {c: data.get(c, "") if data else "" for c in self._df.columns}
        self._df = pd.concat(
            [self._df, pd.DataFrame([new])], ignore_index=True)
        self._modified = True
        self._invalidate_cache()

    def append_dataframe(self, df: pd.DataFrame) -> int:
        """Append rows from a DataFrame, aligning columns. Returns new count."""
        if self._df is None or df.empty:
            return self.row_count
        old_count = self.row_count
        cols = list(self._df.columns)
        records = []
        for _, row in df.iterrows():
            rec = {}
            for c in cols:
                rec[c] = str(row[c]) if c in df.columns and pd.notna(row[c]) else ""
            records.append(rec)
        self._df = pd.concat(
            [self._df, pd.DataFrame(records, columns=cols)], ignore_index=True)
        # Mark all new rows as imported
        for i in range(old_count, self.row_count):
            self._imported.add(i)
        self._modified = True
        self._invalidate_cache()
        return self.row_count

    def is_imported(self, row: int) -> bool:
        return row in self._imported

    def delete_rows(self, indices: List[int]):
        if self._df is None or not indices:
            return
        self._df.drop(self._df.index[list(indices)], inplace=True)
        self._df.reset_index(drop=True, inplace=True)
        self._modified = True
        self._invalidate_cache()

    # ── AI fill tracking ──
    def mark_ai_filled(self, row: int, col_name: str):
        self._ai_filled.add((row, col_name))

    def is_ai_filled(self, row: int, col_name: str) -> bool:
        return (row, col_name) in self._ai_filled

    def clear_ai_marks(self):
        self._ai_filled.clear()

    # ── Fast native-Python search ──
    def search(self, keyword: str, columns: Optional[List[str]] = None) -> pd.DataFrame:
        if self._df is None:
            return pd.DataFrame()
        kw = keyword.lower()
        cols = columns or list(self._df.columns)
        n = len(self._df)
        mask = np.zeros(n, dtype=bool)
        for c in cols:
            if c in self._df.columns:
                arr = self._df[c].values
                mask |= np.array([kw in str(v).lower() for v in arr])
        return self._df.iloc[mask].copy() if mask.any() else pd.DataFrame(columns=self._df.columns)

    def filter_by(self, column: str, value: str) -> pd.DataFrame:
        if self._df is None or column not in self._df.columns:
            return pd.DataFrame()
        vl = value.lower()
        arr = self._df[column].values
        mask = np.array([vl in str(v).lower() for v in arr])
        return self._df.iloc[mask].copy() if mask.any() else pd.DataFrame(columns=self._df.columns)

    # ── Cached stats ──
    def get_stats(self) -> dict:
        if self._stats_cache is not None:
            return self._stats_cache
        if self._df is None:
            return {}
        stats = {"total_rows": len(self._df)}
        for col in self._df.columns:
            arr = self._df[col].values
            filled = sum(1 for v in arr if str(v).strip() != "")
            stats[f"{col}_filled"] = int(filled)
            stats[f"{col}_missing"] = int(len(self._df) - filled)
        self._stats_cache = stats
        return stats

    # ── Cached unique values ──
    def get_unique_values(self, column: str, limit: int = 200) -> List[str]:
        if column in self._unique_cache:
            return self._unique_cache[column]
        if self._df is None or column not in self._df.columns:
            return []
        na_lower = {s.lower() for s in NA_PLACEHOLDERS}
        seen: Set[str] = set()
        for v in self._df[column].values:
            s = str(v).strip()
            if s.lower() not in na_lower and s not in seen:
                seen.add(s)
        result = sorted(seen)[:limit]
        self._unique_cache[column] = result
        return result

    # ── Persistence ──
    def save(self, filepath: Optional[str] = None) -> bool:
        path = filepath or self._filepath
        if not path or self._df is None:
            return False
        try:
            self._df.to_csv(path, sep="\t", index=False)
            self._modified = False
            self._filepath = path
            return True
        except Exception:
            return False

    def export_excel(self, filepath: str) -> bool:
        if self._df is None:
            return False
        try:
            self._df.to_excel(filepath, index=False, engine="openpyxl")
            return True
        except Exception:
            return False

    def import_file(self, filepath: str) -> bool:
        try:
            if filepath.endswith((".xlsx", ".xls")):
                df = pd.read_excel(filepath, dtype=str, keep_default_na=False)
            elif filepath.endswith(".csv"):
                df = pd.read_csv(filepath, dtype=str, keep_default_na=False)
            else:
                df = pd.read_csv(filepath, sep="\t", dtype=str, keep_default_na=False)
            df.fillna("", inplace=True)
            self._df = df.astype(object)
            self._filepath = filepath
            self._modified = False
            self._invalidate_cache()
            return True
        except Exception:
            return False

    @staticmethod
    def is_missing(val) -> bool:
        if val is None:
            return True
        s = str(val).strip().lower()
        return s in NA_PLACEHOLDERS or s == ""
