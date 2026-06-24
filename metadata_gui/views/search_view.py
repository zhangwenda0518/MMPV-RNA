"""Search panel — online SRA/GSA query + local metadata filtering + preview."""

import pandas as pd
import numpy as np

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton,
    QTableView, QHeaderView, QAbstractItemView, QLabel,
    QComboBox, QGroupBox, QGridLayout, QProgressBar,
    QMessageBox, QSplitter, QTextEdit,
)
from PySide6.QtCore import Qt, QThread, Signal, QAbstractTableModel, QModelIndex, QUrl
from PySide6.QtGui import QAction, QColor, QBrush, QFont, QDesktopServices

from models.data_store import MetadataStore
from views.metadata_table import LINK_COLUMNS, _build_url


# ── Background search worker ──────────────────────
class SearchWorker(QThread):
    finished = Signal(dict)
    progress = Signal(str)

    def __init__(self, db: str, query: str, source: str):
        super().__init__()
        self._db = db
        self._query = query
        self._source = source

    def run(self):
        from controllers.search_bridge import search_sra, search_gsa, search_both
        import os
        outdir = os.path.join(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))),
            "..", "search_results")
        os.makedirs(outdir, exist_ok=True)
        self.progress.emit(f"Searching {self._db.upper()} for '{self._query}'...")
        try:
            if self._db in ("sra", "both"):
                self.progress.emit("Querying NCBI SRA...")
                r_sra = search_sra(self._query, self._source, outdir)
                self.progress.emit(f"SRA found: {len(r_sra.get('df', pd.DataFrame()))} runs")
            if self._db in ("gsa", "both"):
                self.progress.emit("Querying CNCB GSA...")
                r_gsa = search_gsa(self._query, self._source, outdir)
                self.progress.emit(f"GSA found: {len(r_gsa.get('df', pd.DataFrame()))} runs")
            if self._db == "sra":
                r_sra["db"] = "sra"; self.finished.emit(r_sra)
            elif self._db == "gsa":
                r_gsa["db"] = "gsa"; self.finished.emit(r_gsa)
            else:
                self.finished.emit({"sra": r_sra, "gsa": r_gsa, "ok": True, "db": "both"})
        except Exception as e:
            self.finished.emit({"ok": False, "error": str(e), "db": self._db})


# ── Simple read-only model for search results ─────
class SearchResultModel(QAbstractTableModel):
    """Lightweight read-only model for displaying search results."""

    def __init__(self):
        super().__init__()
        self._df = pd.DataFrame()
        self._cols = []

    def set_dataframe(self, df: pd.DataFrame):
        self.beginResetModel()
        self._df = df.fillna("")
        self._cols = list(self._df.columns)
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()) -> int:
        return len(self._df)

    def columnCount(self, parent=QModelIndex()) -> int:
        return len(self._cols)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        r, c = index.row(), index.column()
        if r >= len(self._df) or c >= len(self._cols):
            return None
        val = str(self._df.iat[r, c])
        col_name = self._cols[c] if c < len(self._cols) else ""

        if role in (Qt.DisplayRole, Qt.EditRole):
            return val
        elif role == Qt.ForegroundRole:
            if col_name in LINK_COLUMNS and val.strip():
                return QBrush(QColor(0, 100, 200))
            return None
        elif role == Qt.FontRole:
            if col_name in LINK_COLUMNS and val.strip():
                font = QFont()
                font.setUnderline(True)
                return font
            return None
        elif role == Qt.ToolTipRole:
            url = _build_url(col_name, val)
            if url:
                return f"[{col_name}] {val}\n\nClick to open: {url}"
            return f"[{col_name}] {val}"
        return None

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self._cols[section] if section < len(self._cols) else ""
            else:
                return str(section + 1)
        return None

    @property
    def dataframe(self) -> pd.DataFrame:
        return self._df


# ── Search panel ──────────────────────────────────
# ── Deep extract worker ──────────────────────────
class DeepExtractWorker(QThread):
    finished = Signal(str)    # tsv_path
    error = Signal(str)
    progress = Signal(str)

    def __init__(self, run_ids: list, api_key: str = "", ncbi_api: str = ""):
        super().__init__()
        self._ids = run_ids
        self._api_key = api_key
        self._ncbi_api = ncbi_api

    def run(self):
        from controllers.search_bridge import deep_extract
        import os
        outdir = os.path.join(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))),
            "..", "info_results")
        os.makedirs(outdir, exist_ok=True)
        self.progress.emit(f"Extracting {len(self._ids)} runs via gsa_sra.info.py...")
        r = deep_extract(self._ids, outdir,
                         deepseek_api=self._api_key or None,
                         ncbi_api=self._ncbi_api or None)
        if r["ok"]:
            self.finished.emit(r["tsv_path"])
        else:
            self.error.emit(r["error"])


class SearchPanel(QWidget):
    """Search tab: online DB search + preview + local filters."""

    import_requested = Signal(pd.DataFrame)
    send_to_info = Signal(list)   # emits Run ID list for Info tab

    def __init__(self, store: MetadataStore):
        super().__init__()
        self._store = store
        self._worker = None
        self._search_result_df = pd.DataFrame()
        self._setup_ui()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        # ═══ Section 1: Online Search ═══════════════
        online_group = QGroupBox("Online Database Search (SRA / GSA)")
        olay = QVBoxLayout(online_group)

        r1 = QHBoxLayout()
        r1.addWidget(QLabel("Species:"))
        self._species_input = QLineEdit()
        self._species_input.setPlaceholderText("e.g. Lycium barbarum")
        self._species_input.returnPressed.connect(self._start_online_search)
        r1.addWidget(self._species_input, 2)

        r1.addWidget(QLabel("Source:"))
        self._source_combo = QComboBox()
        self._source_combo.setEditable(True)
        self._source_combo.addItems(["TRANSCRIPTOMIC", "All", "GENOMIC", "METAGENOMIC", "OTHER"])
        self._source_combo.setCurrentText("TRANSCRIPTOMIC")
        r1.addWidget(self._source_combo, 1)

        r1.addWidget(QLabel("DB:"))
        self._db_combo = QComboBox()
        self._db_combo.addItems(["both", "sra", "gsa"])
        r1.addWidget(self._db_combo)

        self._search_btn = QPushButton("Search Online")
        self._search_btn.setDefault(True)
        self._search_btn.clicked.connect(self._start_online_search)
        r1.addWidget(self._search_btn)

        self._stop_btn = QPushButton("Stop")
        self._stop_btn.setVisible(False)
        self._stop_btn.clicked.connect(self._stop_search)
        self._stop_btn.setStyleSheet(
            "QPushButton { color: #C44E52; font-weight: bold; }")
        r1.addWidget(self._stop_btn)
        olay.addLayout(r1)

        self._progress = QProgressBar()
        self._progress.setVisible(False)
        self._progress.setTextVisible(False)
        olay.addWidget(self._progress)

        # Search log output area
        self._search_log = QTextEdit()
        self._search_log.setReadOnly(True)
        self._search_log.setMaximumHeight(80)
        self._search_log.setPlaceholderText("Search progress will appear here...")
        self._search_log.setStyleSheet("QTextEdit { font-family: Consolas; font-size: 11px; background: #f8f8f8; }")
        olay.addWidget(self._search_log)

        # Import row
        ir = QHBoxLayout()
        self._result_count = QLabel("")
        ir.addWidget(self._result_count)
        ir.addStretch()
        self._import_btn = QPushButton("Import to Browse Table")
        self._import_btn.setVisible(False)
        self._import_btn.clicked.connect(self._do_import)
        self._import_btn.setStyleSheet("QPushButton { font-weight: bold; }")
        ir.addWidget(self._import_btn)
        # Send to Info button
        self._send_info_btn = QPushButton("Send to Info")
        self._send_info_btn.setVisible(False)
        self._send_info_btn.clicked.connect(self._send_to_info)
        self._send_info_btn.setStyleSheet(
            "QPushButton { color: #0072B2; }")
        self._send_info_btn.setToolTip(
            "Send Run IDs to Info tab for deep extraction via gsa_sra.info.py")
        ir.addWidget(self._send_info_btn)
        olay.addLayout(ir)

        layout.addWidget(online_group)

        # ═══ Online results table ═══════════════════
        self._online_model = SearchResultModel()
        self._online_table = QTableView()
        self._online_table.setModel(self._online_model)
        self._online_table.setSortingEnabled(True)
        self._online_table.setAlternatingRowColors(True)
        self._online_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._online_table.horizontalHeader().setStretchLastSection(True)
        self._online_table.setMaximumHeight(200)
        self._online_table.setEditTriggers(QAbstractItemView.NoEditTriggers)

        online_copy = QAction("Copy Cell", self._online_table)
        online_copy.setShortcut("Ctrl+C")
        online_copy.triggered.connect(self._copy_online_cell)
        self._online_table.addAction(online_copy)
        self._online_table.clicked.connect(self._on_online_click)

        layout.addWidget(self._online_table)

        # ═══ Section 2: Local Filters ═══════════════
        local_group = QGroupBox("Local Metadata Filter")
        llay = QVBoxLayout(local_group)

        qs = QHBoxLayout()
        self._quick_input = QLineEdit()
        self._quick_input.setPlaceholderText("Search loaded records... (e.g. 'leaf', 'Lycium')")
        self._quick_input.setClearButtonEnabled(True)
        self._quick_input.returnPressed.connect(self._do_local_filter)
        qs.addWidget(self._quick_input, 1)
        filter_btn = QPushButton("Filter")
        filter_btn.clicked.connect(self._do_local_filter)
        qs.addWidget(filter_btn)
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self._clear_local)
        qs.addWidget(clear_btn)
        llay.addLayout(qs)

        grid = QGridLayout()
        grid.setSpacing(4)
        self._filters = {}
        filter_fields = [
            "ScientificName", "Tissue", "Source", "Location",
            "BioProject", "LibrarySource", "CenterName",
        ]
        for i, field in enumerate(filter_fields):
            row, col = i % 3, (i // 3) * 2
            lbl = QLabel(field + ":")
            combo = QComboBox()
            combo.setEditable(True)
            combo.addItem("(any)")
            vals = self._store.get_unique_values(field, 80)
            for v in vals:
                combo.addItem(v)
            combo.setMinimumWidth(150)
            combo.currentTextChanged.connect(self._do_local_filter)
            grid.addWidget(lbl, row, col)
            grid.addWidget(combo, row, col + 1)
            self._filters[field] = combo
        llay.addLayout(grid)

        self._local_count = QLabel("")
        llay.addWidget(self._local_count)

        layout.addWidget(local_group)

        # ═══ Local results table ════════════════════
        from views.metadata_table import MetadataTableModel
        self._local_model = MetadataTableModel(self._store)
        self._local_table = QTableView()
        self._local_table.setModel(self._local_model)
        self._local_table.setSortingEnabled(True)
        self._local_table.setAlternatingRowColors(True)
        self._local_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._local_table.horizontalHeader().setStretchLastSection(True)
        self._local_table.setEditTriggers(QAbstractItemView.NoEditTriggers)

        local_copy = QAction("Copy Cell", self._local_table)
        local_copy.setShortcut("Ctrl+C")
        local_copy.triggered.connect(self._copy_local_cell)
        self._local_table.addAction(local_copy)
        self._local_table.clicked.connect(self._on_local_click)

        layout.addWidget(self._local_table, 1)

    # ── Online search ──────────────────────────────
    def _start_online_search(self):
        species = self._species_input.text().strip()
        if not species:
            self._search_log.append("Please enter a species name.")
            return
        self._search_btn.setEnabled(False)
        self._stop_btn.setVisible(True)
        self._import_btn.setVisible(False)
        self._send_info_btn.setVisible(False)
        self._progress.setVisible(True)
        self._progress.setRange(0, 0)
        self._search_log.clear()
        self._search_log.append(f"> Searching for: {species}")
        self._result_count.setText("")
        self._search_result_df = pd.DataFrame()
        self._online_model.set_dataframe(pd.DataFrame())

        db = self._db_combo.currentText()
        src = self._source_combo.currentText().strip()
        self._search_log.append(f"  DB: {db.upper()}, Source: {src or 'All'}")
        self._worker = SearchWorker(db, species, src)
        self._worker.progress.connect(self._on_search_progress)
        self._worker.finished.connect(self._on_search_done)
        self._worker.start()

    def _on_search_progress(self, msg: str):
        self._search_log.append(f"  {msg}")

    def _stop_search(self):
        if self._worker and self._worker.isRunning():
            self._worker.terminate()
            self._worker.wait(1000)
        self._search_btn.setEnabled(True)
        self._stop_btn.setVisible(False)
        self._progress.setVisible(False)
        self._search_log.append("  [STOPPED by user]")

    def _on_search_done(self, result: dict):
        self._search_btn.setEnabled(True)
        self._stop_btn.setVisible(False)
        self._progress.setVisible(False)

        if not result.get("ok"):
            err = result.get("error", "")
            self._search_log.append(f"  FAILED: {err}")
            self._result_count.setText(f"Search failed")
            return

        db = result.get("db", "")
        if db == "both":
            sra = result.get("sra", {})
            gsa = result.get("gsa", {})
            sra_n = len(sra.get("df", pd.DataFrame()))
            gsa_n = len(gsa.get("df", pd.DataFrame()))
            self._search_log.append(f"  SRA: {sra_n} runs, GSA: {gsa_n} runs")
            dfs = [d for d in [sra.get("df", pd.DataFrame()), gsa.get("df", pd.DataFrame())] if not d.empty]
            self._search_result_df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
        else:
            self._search_result_df = result.get("df", pd.DataFrame())
            self._search_log.append(f"  {db.upper()}: {len(self._search_result_df)} runs")
            if not self._search_result_df.empty:
                cols = list(self._search_result_df.columns)
                self._search_log.append(f"  Columns: {', '.join(cols[:6])}...")

        total = len(self._search_result_df)
        if total > 0:
            self._online_model.set_dataframe(self._search_result_df)
            self._search_log.append(f"  Done: {total} total results")
            self._result_count.setText(f"{total} results")
            self._import_btn.setVisible(True)
            self._send_info_btn.setVisible(True)
        else:
            self._online_model.set_dataframe(pd.DataFrame())
            self._search_log.append(f"  No results found")
            self._result_count.setText("No results")
            self._import_btn.setVisible(False)
            self._send_info_btn.setVisible(False)

    def _do_import(self):
        if self._search_result_df.empty:
            return
        self.import_requested.emit(self._search_result_df)
        self._import_btn.setVisible(False)
        self._send_info_btn.setVisible(False)
        self._result_count.setText(
            f"Imported {len(self._search_result_df)} records to S-Browse.")

    def _send_to_info(self):
        """Send Run IDs from search results to Info tab."""
        if self._search_result_df.empty or "Run" not in self._search_result_df.columns:
            return
        ids = [str(r).strip() for r in self._search_result_df["Run"].tolist()
               if str(r).strip()]
        if ids:
            self.send_to_info.emit(ids)

    # ── Local filter ───────────────────────────────
    def _do_local_filter(self):
        df = self._store.dataframe
        if df is None or len(df) == 0:
            return
        kw = self._quick_input.text().strip()
        mask = pd.Series(True, index=df.index)
        if kw:
            kwl = kw.lower()
            kmask = pd.Series(False, index=df.index)
            for c in df.columns:
                kmask |= df[c].apply(lambda x: kwl in str(x).lower() if pd.notna(x) else False)
            mask &= kmask
        for field, combo in self._filters.items():
            val = combo.currentText().strip()
            if val and val != "(any)":
                vl = val.lower()
                m = df[field].apply(lambda x: vl in str(x).lower() if pd.notna(x) else False)
                mask &= m
        n = mask.sum()
        self._local_count.setText(f"Showing: {n} / {len(df)} records")

    def _clear_local(self):
        self._quick_input.clear()
        for combo in self._filters.values():
            combo.setCurrentIndex(0)
        self._local_count.clear()

    def _on_online_click(self, index):
        if not index.isValid():
            return
        col_name = self._online_model._cols[index.column()] if index.column() < len(self._online_model._cols) else ""
        if col_name not in LINK_COLUMNS:
            return
        val = str(self._online_model._df.iat[index.row(), index.column()]).strip()
        url = _build_url(col_name, val)
        if url:
            QDesktopServices.openUrl(QUrl(url))

    def _on_local_click(self, index):
        if not index.isValid():
            return
        # local table uses MetadataTableModel, check column name
        cols = self._store.columns
        col_name = cols[index.column()] if index.column() < len(cols) else ""
        if col_name not in LINK_COLUMNS:
            return
        val = self._store.get_cell(index.row(), index.column()).strip()
        url = _build_url(col_name, val)
        if url:
            QDesktopServices.openUrl(QUrl(url))

    def _copy_online_cell(self):
        idx = self._online_table.currentIndex()
        if idx.isValid():
            from PySide6.QtWidgets import QApplication
            QApplication.clipboard().setText(str(self._online_model.data(idx, Qt.DisplayRole)))

    def _copy_local_cell(self):
        idx = self._local_table.currentIndex()
        if idx.isValid():
            from PySide6.QtWidgets import QApplication
            QApplication.clipboard().setText(str(self._local_model.data(idx, Qt.DisplayRole)))
