"""Info panel — Run ID list input for gsa_sra.info.py deep extraction."""

import os
import re
import pandas as pd

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QPushButton,
    QLabel, QGroupBox, QProgressBar, QTableView, QHeaderView,
    QAbstractItemView, QMessageBox,
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QAction

from models.data_store import MetadataStore
from views.search_view import DeepExtractWorker, SearchResultModel


class InfoPanel(QWidget):
    """Info tab: paste Run IDs or pull from Search, run gsa_sra.info.py."""

    import_requested = Signal(pd.DataFrame)

    def __init__(self, info_store: MetadataStore, search_store: MetadataStore):
        super().__init__()
        self._info_store = info_store
        self._search_store = search_store
        self._worker = None
        self._result_df = pd.DataFrame()
        self._setup_ui()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        # ── Section 1: Input ───────────────────────
        inp_group = QGroupBox("Run ID Input")
        ilay = QVBoxLayout(inp_group)

        ilay.addWidget(QLabel(
            "Paste Run IDs (comma/space/newline separated), "
            "or pull from Search results:"))

        self._text_edit = QTextEdit()
        self._text_edit.setPlaceholderText(
            "SRR31651831\nCRR1126132\nERR13512478\n...")
        self._text_edit.setMaximumHeight(120)
        ilay.addWidget(self._text_edit)

        btn_row = QHBoxLayout()
        btn_row.addWidget(QLabel("Source:"))

        self._from_search_btn = QPushButton("Pull from Search Results")
        self._from_search_btn.clicked.connect(self._pull_from_search)
        self._from_search_btn.setStyleSheet(
            "QPushButton { color: #0072B2; }")
        btn_row.addWidget(self._from_search_btn)

        btn_row.addStretch()

        self._extract_btn = QPushButton("Deep Extract (gsa_sra.info.py)")
        self._extract_btn.setDefault(True)
        self._extract_btn.clicked.connect(self._start_extract)
        self._extract_btn.setStyleSheet(
            "QPushButton { font-weight: bold; background: #0072B2; "
            "color: white; padding: 6px 16px; border-radius: 4px; }")
        btn_row.addWidget(self._extract_btn)
        ilay.addLayout(btn_row)

        layout.addWidget(inp_group)

        # ── Progress ───────────────────────────────
        self._progress = QProgressBar()
        self._progress.setVisible(False)
        self._progress.setTextVisible(False)
        layout.addWidget(self._progress)

        # ── Status + import ────────────────────────
        sr = QHBoxLayout()
        self._status_label = QLabel("")
        sr.addWidget(self._status_label)
        sr.addStretch()
        self._import_btn = QPushButton("Import to I-Browse")
        self._import_btn.setVisible(False)
        self._import_btn.clicked.connect(self._do_import)
        self._import_btn.setStyleSheet("QPushButton { font-weight: bold; }")
        sr.addWidget(self._import_btn)
        layout.addLayout(sr)

        # ── Results preview ───────────────────────
        self._result_model = SearchResultModel()
        self._result_table = QTableView()
        self._result_table.setModel(self._result_model)
        self._result_table.setSortingEnabled(True)
        self._result_table.setAlternatingRowColors(True)
        self._result_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._result_table.horizontalHeader().setStretchLastSection(True)
        self._result_table.setMaximumHeight(200)
        self._result_table.setEditTriggers(QAbstractItemView.NoEditTriggers)

        copy_act = QAction("Copy Cell", self._result_table)
        copy_act.setShortcut("Ctrl+C")
        copy_act.triggered.connect(self._copy_cell)
        self._result_table.addAction(copy_act)

        layout.addWidget(self._result_table, 1)

    def _pull_from_search(self):
        """Pull Run IDs from Search store."""
        if not self._search_store.is_loaded:
            self._status_label.setText("No Search data available.")
            return
        df = self._search_store.dataframe
        if "Run" not in df.columns:
            return
        ids = [str(r).strip() for r in df["Run"].tolist()
               if str(r).strip() and str(r).strip().lower() != "not_provided"]
        if not ids:
            self._status_label.setText("No Run IDs found in Search data.")
            return
        self._text_edit.setPlainText("\n".join(ids))
        self._status_label.setText(f"Pulled {len(ids)} Run IDs from Search.")

    def _start_extract(self):
        text = self._text_edit.toPlainText().strip()
        if not text:
            self._status_label.setText("Please paste Run IDs first.")
            return
        ids = re.split(r'[,\s\n;]+', text)
        ids = [i.strip() for i in ids if i.strip()]
        if not ids:
            return

        self._extract_btn.setEnabled(False)
        self._import_btn.setVisible(False)
        self._progress.setVisible(True)
        self._progress.setRange(0, 0)
        self._result_df = pd.DataFrame()
        self._result_model.set_dataframe(pd.DataFrame())

        api_key = os.environ.get("DEEPSEEK_API_KEY", "")
        ncbi_api = os.environ.get("NCBI_API_KEY", "")

        self._worker = DeepExtractWorker(ids, api_key, ncbi_api)
        self._worker.progress.connect(lambda m: self._status_label.setText(m))
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, tsv_path: str):
        self._extract_btn.setEnabled(True)
        self._progress.setVisible(False)
        try:
            df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
            df.fillna("", inplace=True)
            self._result_df = df
            self._result_model.set_dataframe(df)
            self._status_label.setText(
                f"Extracted {len(df)} records. Click Import to load into I-Browse.")
            self._import_btn.setVisible(True)
        except Exception as e:
            self._status_label.setText(f"Failed to load results: {e}")

    def _on_error(self, msg: str):
        self._extract_btn.setEnabled(True)
        self._progress.setVisible(False)
        self._status_label.setText(f"Error: {msg}")

    def _do_import(self):
        if self._result_df.empty:
            return
        self.import_requested.emit(self._result_df)
        self._import_btn.setVisible(False)
        self._status_label.setText(f"Imported {len(self._result_df)} records to I-Browse.")

    def _copy_cell(self):
        idx = self._result_table.currentIndex()
        if idx.isValid():
            from PySide6.QtWidgets import QApplication
            QApplication.clipboard().setText(
                str(self._result_model.data(idx, Qt.DisplayRole)))
