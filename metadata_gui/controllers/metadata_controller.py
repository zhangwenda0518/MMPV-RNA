"""Metadata controller — bridges data store and views with business logic."""

import os
from typing import Optional, List, Dict

import pandas as pd

from models.data_store import MetadataStore


class MetadataController:
    """Controller handling metadata operations and view coordination."""

    def __init__(self, store: MetadataStore):
        self._store = store

    # ── Load / Save ─────────────────────────────────

    def load_file(self, path: str) -> bool:
        return self._store.load(path)

    def load_default_data(self) -> bool:
        """Try to auto-discover and load data."""
        possible_dirs = []
        # From controller's location
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        project = os.path.dirname(base)
        possible_dirs.append(os.path.join(
            project, "public_metadata_pipeline",
            "public_data_pipeline_output", "info"))
        # Direct
        possible_dirs.append(os.path.join(
            project, "public_metadata_pipeline_output", "info"))

        for d in possible_dirs:
            core = os.path.join(d, "Global_Unified_Metadata_Core13.tsv")
            if os.path.isfile(core):
                ok = self._store.load(core)
                full = os.path.join(d, "Global_Unified_Metadata_Full.tsv")
                if os.path.isfile(full):
                    self._store.load_full(full)
                return ok
        return False

    def save(self, path: Optional[str] = None) -> bool:
        return self._store.save(path)

    def export_excel(self, path: str) -> bool:
        return self._store.export_excel(path)

    def import_file(self, path: str) -> bool:
        return self._store.import_file(path)

    # ── Query ──────────────────────────────────────

    def search(self, keyword: str) -> pd.DataFrame:
        return self._store.search(keyword)

    def filter_by(self, column: str, value: str) -> pd.DataFrame:
        return self._store.filter_by(column, value)

    def get_stats(self) -> dict:
        return self._store.get_stats()

    def get_unique_values(self, column: str) -> List[str]:
        return self._store.get_unique_values(column)

    # ── Edit ───────────────────────────────────────

    def update_cell(self, row: int, col: int, value: str):
        self._store.set_cell(row, col, value)

    def update_row(self, row: int, data: Dict[str, str]):
        self._store.set_row(row, data)

    def delete_rows(self, indices: List[int]):
        self._store.delete_rows(indices)

    def add_row(self, data: Dict[str, str] = None) -> int:
        old_count = self._store.row_count
        self._store.add_row(data)
        return old_count  # new row index

    # ── Batch completion ───────────────────────────

    def fill_missing_from_cache(self, cache_dir: str) -> int:
        """Try to fill missing values from JSON cache files."""
        if not os.path.isdir(cache_dir):
            return 0

        import json
        filled = 0
        mapping = {
            "SubmissionDate": "ReleaseDate",
            "Organization": "CenterName",
            "PRJ": "BioProject",
            "TaxID": "TaxID",
            "ScientificName": "ScientificName",
        }

        for idx in range(self._store.row_count):
            run = self._store.get_cell(
                idx, self._store.dataframe.columns.get_loc("Run")
                if "Run" in self._store.dataframe.columns else 0)
            if not run:
                continue

            cache_path = os.path.join(cache_dir, f"{run}.json")
            if not os.path.isfile(cache_path):
                continue

            try:
                with open(cache_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
            except Exception:
                continue

            for json_key, col_name in mapping.items():
                if col_name not in self._store.dataframe.columns:
                    continue
                val = data.get(json_key, "")
                if val and self._store.is_missing(
                    self._store.get_cell(
                        idx,
                        self._store.dataframe.columns.get_loc(col_name)
                    )
                ):
                    self._store.set_cell(
                        idx,
                        self._store.dataframe.columns.get_loc(col_name),
                        str(val)
                    )
                    filled += 1

        return filled
