"""Editable metadata table with search/filter and missing-value highlighting."""

from typing import Optional, List

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton,
    QTableView, QHeaderView, QAbstractItemView, QLabel,
    QComboBox, QMessageBox,
)
from PySide6.QtCore import (
    Qt, QAbstractTableModel, QModelIndex, Signal,
    QSortFilterProxyModel, QUrl,
)
from PySide6.QtGui import QColor, QBrush, QFont, QAction, QDesktopServices

from models.data_store import MetadataStore


LINK_COLUMNS = {"Run", "BioProject", "PMID", "BioSample", "TaxID"}


def _build_url(col_name: str, val: str) -> str:
    """Build external URL for a given column + value."""
    v = val.strip()
    if not v:
        return ""
    col_lower = col_name.lower()

    if col_lower == "run":
        if v.startswith(("SRR", "ERR", "DRR", "SRX", "ERX", "DRX")):
            return f"https://www.ncbi.nlm.nih.gov/sra/{v}"
        elif v.startswith("CRR"):
            return f"https://ngdc.cncb.ac.cn/gsa/search?searchTerm={v}"
        return f"https://www.ncbi.nlm.nih.gov/sra/?term={v}"

    if col_lower == "bioproject":
        if v.startswith("PRJNA"):
            return f"https://www.ncbi.nlm.nih.gov/bioproject/{v}"
        elif v.startswith("PRJCA") or v.startswith("PRJEB"):
            return f"https://ngdc.cncb.ac.cn/bioproject/browse/{v}"
        return f"https://www.ncbi.nlm.nih.gov/bioproject/?term={v}"

    if col_lower == "pmid":
        if v.isdigit():
            return f"https://pubmed.ncbi.nlm.nih.gov/{v}/"
        return f"https://pubmed.ncbi.nlm.nih.gov/?term={v}"

    if col_lower == "biosample":
        if v.startswith("SAMN") or v.startswith("SAME"):
            return f"https://www.ncbi.nlm.nih.gov/biosample/{v}"
        return f"https://ngdc.cncb.ac.cn/biosample/{v}"

    if col_lower == "taxid":
        if v.isdigit():
            return f"https://www.ncbi.nlm.nih.gov/taxonomy/{v}"
        return ""

    return ""


class MetadataTableModel(QAbstractTableModel):
    """Qt model adapter for MetadataStore DataFrame, with link support."""

    data_changed = Signal()

    def __init__(self, store: MetadataStore):
        super().__init__()
        self._store = store

    def rowCount(self, parent=QModelIndex()) -> int:
        return self._store.row_count

    def columnCount(self, parent=QModelIndex()) -> int:
        return self._store.column_count

    def _get_col_name(self, col: int) -> str:
        cols = self._store.columns
        return cols[col] if col < len(cols) else ""

    def _is_link_col(self, col: int) -> bool:
        return self._get_col_name(col) in LINK_COLUMNS

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        row, col = index.row(), index.column()
        if row >= self._store.row_count or col >= self._store.column_count:
            return None
        val = self._store.get_cell(row, col)
        col_name = self._get_col_name(col)

        if role == Qt.DisplayRole:
            return val if val != "" else ""
        elif role == Qt.EditRole:
            return val
        elif role == Qt.BackgroundRole:
            if self._store.is_ai_filled(row, col_name):
                return QBrush(QColor(220, 255, 220))  # light green: AI filled
            if self._store.is_imported(row):
                return QBrush(QColor(230, 240, 255))  # light blue: imported
            if self._store.is_missing(val):
                return QBrush(QColor(255, 255, 200))  # light yellow: missing
            return None
        elif role == Qt.ForegroundRole:
            if self._store.is_missing(val):
                return QBrush(QColor(180, 50, 50))
            if self._is_link_col(col) and val.strip():
                return QBrush(QColor(0, 100, 200))  # blue link
            return None
        elif role == Qt.FontRole:
            if self._is_link_col(col) and val.strip():
                font = QFont()
                font.setUnderline(True)
                return font
            return None
        elif role == Qt.ToolTipRole:
            url = _build_url(col_name, val)
            if url:
                return f"[{col_name}] {val}\n\nClick to open in browser: {url}"
            return f"[{col_name}] {val}" if val else f"[{col_name}] <empty>"
        return None

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                cols = self._store.columns
                return cols[section] if section < len(cols) else str(section)
            else:
                return str(section + 1)  # 1-based row numbers
        return None

    def flags(self, index):
        default = Qt.ItemIsEnabled | Qt.ItemIsSelectable
        if index.isValid():
            default |= Qt.ItemIsEditable
        return default

    def setData(self, index, value, role=Qt.EditRole):
        if role == Qt.EditRole and index.isValid():
            self._store.set_cell(index.row(), index.column(), str(value))
            self.dataChanged.emit(index, index, [Qt.DisplayRole, Qt.BackgroundRole])
            self.data_changed.emit()
            return True
        return False

    def refresh(self):
        """Full model reset."""
        self.beginResetModel()
        self.endResetModel()


class StableSortProxy(QSortFilterProxyModel):
    """Proxy that preserves original row numbers in vertical header."""

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Vertical:
            src_idx = self.mapToSource(self.index(section, 0))
            return str(src_idx.row() + 1) if src_idx.isValid() else str(section + 1)
        return super().headerData(section, orientation, role)


class MetadataTableView(QWidget):
    """Browse tab — sortable, filterable, editable table."""

    row_selected = Signal(int)   # emits row index

    def __init__(self, store: MetadataStore):
        super().__init__()
        self._store = store
        self._filtered_df = None
        self._setup_ui()
        self.refresh()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        # ── Top bar: filter ─────────────────────────
        top = QHBoxLayout()

        self._filter_input = QLineEdit()
        self._filter_input.setPlaceholderText("Search all fields... (Enter to filter, Esc to clear)")
        self._filter_input.setClearButtonEnabled(True)
        self._filter_input.returnPressed.connect(self._apply_filter)
        top.addWidget(self._filter_input, 1)

        self._col_combo = QComboBox()
        self._col_combo.addItem("All Columns")
        for c in self._store.columns:
            self._col_combo.addItem(c)
        self._col_combo.currentIndexChanged.connect(self._apply_filter)
        top.addWidget(self._col_combo)

        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self._clear_filter)
        top.addWidget(clear_btn)

        top.addStretch()
        self._fill_btn = QPushButton("Fill from Cache")
        self._fill_btn.clicked.connect(self._fill_missing)
        self._fill_btn.setToolTip("Auto-fill missing values from local JSON web cache")
        top.addWidget(self._fill_btn)

        self._ai_btn = QPushButton("AI Fill Missing")
        self._ai_btn.clicked.connect(self._ai_fill)
        self._ai_btn.setToolTip("Use DeepSeek AI to clean and complete missing metadata")
        self._ai_btn.setStyleSheet(
            "QPushButton { font-weight: bold; color: #0072B2; }")
        top.addWidget(self._ai_btn)

        self._export_btn = QPushButton("Export Table")
        self._export_btn.clicked.connect(self._export_table)
        self._export_btn.setToolTip("Export current table data to TSV/CSV/Excel")
        top.addWidget(self._export_btn)

        self._summary_btn = QPushButton("AI Summary")
        self._summary_btn.clicked.connect(self._ai_summary)
        self._summary_btn.setToolTip("Generate SCI writing summary from data statistics")
        self._summary_btn.setStyleSheet(
            "QPushButton { font-weight: bold; color: #0072B2; }")
        top.addWidget(self._summary_btn)

        layout.addLayout(top)

        # ── Stats bar ───────────────────────────────
        self._stats_label = QLabel()
        layout.addWidget(self._stats_label)

        # ── Table ───────────────────────────────────
        self._model = MetadataTableModel(self._store)
        self._proxy = StableSortProxy()
        self._proxy.setSourceModel(self._model)

        self._table = QTableView()
        self._table.setModel(self._proxy)
        self._table.setSortingEnabled(True)
        self._table.setAlternatingRowColors(True)
        self._table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self._table.horizontalHeader().setMinimumSectionSize(60)
        self._table.verticalHeader().setVisible(True)
        self._table.setEditTriggers(
            QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed)
        self._table.setContextMenuPolicy(Qt.ActionsContextMenu)

        # Context menu actions
        copy_act = QAction("Copy Cell", self._table)
        copy_act.setShortcut("Ctrl+C")
        copy_act.triggered.connect(self._copy_cell)
        self._table.addAction(copy_act)

        delete_act = QAction("Delete Selected Rows", self._table)
        delete_act.setShortcut("Delete")
        delete_act.triggered.connect(self._delete_selected)
        self._table.addAction(delete_act)

        self._table.selectionModel().selectionChanged.connect(self._on_selection)
        self._table.clicked.connect(self._on_cell_clicked)

        self._model.data_changed.connect(self._update_stats)
        layout.addWidget(self._table, 1)

    # ── Filter ─────────────────────────────────────
    def _apply_filter(self):
        kw = self._filter_input.text().strip()
        if not kw:
            self._filtered_df = None
            self._model.refresh()
            self._update_stats()
            return

        col_idx = self._col_combo.currentIndex()
        if col_idx == 0:
            self._filtered_df = self._store.search(kw)
        else:
            col_name = self._col_combo.currentText()
            self._filtered_df = self._store.filter_by(col_name, kw)

        self._model.refresh()
        self._update_stats()

    def _clear_filter(self):
        self._filter_input.clear()
        self._filtered_df = None
        self._model.refresh()
        self._update_stats()

    def _update_stats(self):
        total = self._store.row_count
        if total == 0:
            self._stats_label.setText("No data loaded. Use File > Open to load a metadata TSV.")
            return
        missing_count = 0
        for col in self._store.columns:
            if col in self._store.dataframe.columns:
                try:
                    cnt = int(self._store.dataframe[col].apply(
                        lambda x: self._store.is_missing(x)).sum())
                except (ValueError, TypeError):
                    cnt = 0
                missing_count += cnt
        pct = (missing_count / (total * self._store.column_count) * 100) if total > 0 else 0
        self._stats_label.setText(
            f"Total: {total} rows  |  {self._store.column_count} columns  |  "
            f"Missing cells: {missing_count} ({pct:.1f}%)  |  "
            f"Double-click to edit"
        )

    # ── Selection ──────────────────────────────────
    def _on_selection(self):
        indexes = self._table.selectionModel().selectedRows()
        if indexes:
            source_idx = self._proxy.mapToSource(indexes[0])
            self.row_selected.emit(source_idx.row())

    def _on_cell_clicked(self, index):
        """Open browser for link columns (Run, BioProject, PMID, etc.)."""
        source_idx = self._proxy.mapToSource(index)
        col = source_idx.column()
        cols = self._store.columns
        col_name = cols[col] if col < len(cols) else ""
        if col_name not in LINK_COLUMNS:
            return
        val = self._store.get_cell(source_idx.row(), col).strip()
        if not val:
            return
        url = _build_url(col_name, val)
        if url:
            QDesktopServices.openUrl(QUrl(url))

    @property
    def selected_rows(self) -> List[int]:
        rows = set()
        for idx in self._table.selectionModel().selectedRows():
            rows.add(self._proxy.mapToSource(idx).row())
        return sorted(rows)

    def get_selected_row(self) -> int:
        rows = self.selected_rows
        return rows[0] if rows else -1

    # ── Actions ────────────────────────────────────
    def _copy_cell(self):
        idx = self._table.currentIndex()
        if idx.isValid():
            from PySide6.QtWidgets import QApplication
            val = str(self._proxy.data(idx, Qt.DisplayRole))
            QApplication.clipboard().setText(val)

    def _delete_selected(self):
        rows = self.selected_rows
        if not rows:
            return
        r = QMessageBox.question(
            self, "Delete Rows",
            f"Delete {len(rows)} selected row(s)? This can be undone by reloading.",
            QMessageBox.Yes | QMessageBox.No)
        if r == QMessageBox.Yes:
            self._store.delete_rows(rows)
            self._model.refresh()

    def _fill_missing(self):
        """Auto-fill missing metadata from JSON web cache."""
        from controllers.metadata_controller import MetadataController
        import os

        ctrl = MetadataController(self._store)
        # Find cache directory
        gui_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        project_dir = os.path.dirname(gui_dir)
        cache_dir = os.path.join(
            project_dir, "public_metadata_pipeline",
            "public_data_pipeline_output", "info",
            "GSA_Results", "0_web_cache")

        if not os.path.isdir(cache_dir):
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.information(self, "Fill Missing",
                                    "JSON cache directory not found.\n"
                                    "Run gsa_sra.info.py first to populate cache.")
            return

        filled = ctrl.fill_missing_from_cache(cache_dir)
        self._model.refresh()
        self._update_stats()

        from PySide6.QtWidgets import QMessageBox
        if filled > 0:
            QMessageBox.information(self, "Fill Missing",
                                    f"Filled {filled} missing values from JSON cache.")
        else:
            QMessageBox.information(self, "Fill Missing",
                                    "No missing values could be filled from cache.")

    def _ai_fill(self):
        """Use DeepSeek AI to complete missing metadata fields."""
        from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout,
                                        QLineEdit, QComboBox, QDialogButtonBox,
                                        QProgressBar, QLabel, QMessageBox)
        import os

        # ── API key dialog ────────────────
        dlg = QDialog(self)
        dlg.setWindowTitle("AI Fill — DeepSeek Configuration")
        dlg.setMinimumWidth(420)
        dl = QVBoxLayout(dlg)

        dl.addWidget(QLabel("DeepSeek API Key:"))
        key_edit = QLineEdit()
        key_edit.setEchoMode(QLineEdit.Password)
        key_edit.setPlaceholderText("sk-... or set DEEPSEEK_API_KEY env var")
        # Try env var first
        env_key = os.environ.get("DEEPSEEK_API_KEY", "")
        if env_key:
            key_edit.setText(env_key)
            key_edit.setEnabled(False)
            dl.addWidget(QLabel("(using DEEPSEEK_API_KEY from environment)"))
        dl.addWidget(key_edit)

        dl.addWidget(QLabel("API Base URL:"))
        base_edit = QLineEdit("https://api.deepseek.com")
        dl.addWidget(base_edit)

        dl.addWidget(QLabel("Model:"))
        model_combo = QComboBox()
        model_combo.addItems([
            "deepseek-chat", "deepseek-v4-flash",
            "deepseek-v4-pro", "deepseek-reasoner"])
        dl.addWidget(model_combo)

        dl.addWidget(QLabel(
            "AI will clean and fill missing: Location, Tissue, Source, "
            "Age_GrowthStage, ScientificName, LibrarySource."))

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept)
        btns.rejected.connect(dlg.reject)
        dl.addWidget(btns)

        if dlg.exec() != QDialog.Accepted:
            return

        api_key = key_edit.text().strip()
        api_base = base_edit.text().strip()
        model = model_combo.currentText()

        if not api_key:
            QMessageBox.warning(self, "AI Fill",
                                "Please provide a DeepSeek API key.\n"
                                "Or set DEEPSEEK_API_KEY environment variable.")
            return

        # ── Collect all records with any fillable-column data ──
        from controllers.ai_completer import AICompleteWorker, FILLABLE_COLS, EMPTY_VALS

        records = []
        store_cols = set(self._store.columns)
        for idx in range(self._store.row_count):
            row = self._store.get_row(idx)
            # Send to AI if any fillable column has some data (to clean or fill)
            has_data = any(
                c in store_cols and str(row.get(c, "")).strip()
                for c in FILLABLE_COLS
            )
            if has_data:
                records.append((idx, row))

        if not records:
            QMessageBox.information(self, "AI Fill",
                                    "No records to process.")
            return

        # ── Progress dialog ────────────────
        pdlg = QDialog(self)
        pdlg.setWindowTitle("AI Fill — Processing...")
        pdlg.setMinimumWidth(380)
        pl = QVBoxLayout(pdlg)
        pl.addWidget(QLabel(
            f"Sending {len(records)} records to AI for cleaning & completion..."))
        progress_bar = QProgressBar()
        progress_bar.setRange(0, len(records))
        progress_bar.setTextVisible(False)
        pl.addWidget(progress_bar)
        status_label = QLabel("")
        pl.addWidget(status_label)
        pdlg.show()

        self._ai_worker = AICompleteWorker(
            records, api_key, api_base, model)
        self._ai_worker.progress.connect(
            lambda cur, tot: (
                progress_bar.setValue(cur),
                status_label.setText(f"Processing {cur}/{tot}...")
            ))
        self._ai_worker.error.connect(
            lambda e: status_label.setText(f"Error: {e}"))

        def on_done(filled_count):
            for row_idx, applied in self._ai_worker.results:
                self._store.set_row(row_idx, applied)
                for col_name in applied:
                    self._store.mark_ai_filled(row_idx, col_name)
            self._model.refresh()
            self._update_stats()
            pdlg.accept()
            QMessageBox.information(
                self, "AI Complete",
                f"AI cleaned/filled {filled_count} fields across "
                f"{len(self._ai_worker.results)} records.\n"
                f"Green cells = AI modified.")

        self._ai_worker.finished.connect(on_done)
        self._ai_worker.start()

    def _export_table(self):
        from PySide6.QtWidgets import QFileDialog
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Table", "metadata_export.tsv",
            "TSV (*.tsv);;CSV (*.csv);;Excel (*.xlsx)")
        if not path:
            return
        try:
            if path.endswith(".xlsx"):
                self._store.export_excel(path)
            else:
                self._store.save(path)
            QMessageBox.information(self, "Export", f"Exported to {path}")
        except Exception as e:
            QMessageBox.warning(self, "Export Failed", str(e))

    def _ai_summary(self):
        """Generate SCI writing summary using AI based on data statistics."""
        from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout,
                                        QLineEdit, QTextEdit, QPushButton,
                                        QLabel, QDialogButtonBox)
        import os, json
        from collections import Counter

        df = self._store.dataframe
        if df is None or len(df) == 0:
            QMessageBox.information(self, "AI Summary", "No data available.")
            return

        # ── Collect statistics ──
        stats = {"total_runs": len(df)}
        cols = {c.lower(): c for c in df.columns}

        def col(name):
            key = name.lower()
            return cols.get(key)

        def safe_col(name):
            c = col(name)
            if c and c in df.columns:
                series = df[c].apply(lambda x: str(x).strip())
                series = series[~series.apply(self._store.is_missing)]
                return Counter(series)
            return Counter()

        db_counts = safe_col("Database") if "Database" in cols else Counter()
        species_counts = safe_col("ScientificName")
        tissue_counts = safe_col("Tissue")
        source_counts = safe_col("Source")
        location_counts = safe_col("Location")
        center_counts = safe_col("CenterName")
        age_counts = safe_col("Age_GrowthStage")

        # Collection years
        years = set()
        date_col = col("ReleaseDate") or col("CollectionDate")
        if date_col and date_col in df.columns:
            for v in df[date_col]:
                s = str(v).strip()
                parts = s.replace("/","-").split("-")
                if parts[0].isdigit() and len(parts[0])==4:
                    years.add(parts[0])

        # Data volume
        total_gb = 0; gb_count = 0
        size_col = col("FileSize_GB") or col("FileSize_MB")
        if size_col and size_col in df.columns:
            for v in df[size_col]:
                try:
                    total_gb += float(str(v).strip())
                    gb_count += 1
                except (ValueError, TypeError):
                    pass
        vol_info = ""
        if gb_count > 0:
            vol_info = f"Data volume: {total_gb:.1f} GB total ({gb_count} runs with size data, avg {total_gb/gb_count:.1f} GB/run)"

        stats_text = f"""Total records: {stats['total_runs']}
{vol_info}
Database sources: {', '.join(f'{k}({v})' for k,v in db_counts.most_common())}
Species: {', '.join(f'{k}({v})' for k,v in species_counts.most_common())}
Tissues: {', '.join(f'{k}({v})' for k,v in tissue_counts.most_common())}
Source types: {', '.join(f'{k}({v})' for k,v in source_counts.most_common())}
Locations: {', '.join(f'{k}({v})' for k,v in location_counts.most_common(10))}
Institutions: {', '.join(f'{k}({v})' for k,v in center_counts.most_common(10))}
Growth stages: {', '.join(f'{k}({v})' for k,v in age_counts.most_common(5))}
Collection years: {', '.join(sorted(years)) if years else 'N/A'}"""

        # ── API key dialog ──
        dlg = QDialog(self)
        dlg.setWindowTitle("AI Summary — DeepSeek API")
        dlg.setMinimumWidth(500)
        dl = QVBoxLayout(dlg)

        dl.addWidget(QLabel("DeepSeek API Key:"))
        key_edit = QLineEdit()
        key_edit.setEchoMode(QLineEdit.Password)
        key_edit.setPlaceholderText("sk-... or set DEEPSEEK_API_KEY env var")
        env_key = os.environ.get("DEEPSEEK_API_KEY", "")
        if env_key:
            key_edit.setText(env_key); key_edit.setEnabled(False)
            dl.addWidget(QLabel("(using DEEPSEEK_API_KEY from environment)"))
        dl.addWidget(key_edit)

        dl.addWidget(QLabel("Model:"))
        from PySide6.QtWidgets import QComboBox
        model_combo = QComboBox()
        model_combo.addItems(["deepseek-chat", "deepseek-v4-flash"])
        dl.addWidget(model_combo)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept); btns.rejected.connect(dlg.reject)
        dl.addWidget(btns)

        if dlg.exec() != QDialog.Accepted:
            return

        api_key = key_edit.text().strip()
        if not api_key:
            QMessageBox.warning(self, "AI Summary", "API key required.")
            return

        # ── Call AI ──
        prompt = f"""You are a scientific writer preparing a manuscript for a virome/metagenomics study.
Based on the following metadata statistics of public sequencing runs collected for analysis, write a concise paragraph (150-250 words) suitable for the "Data Collection" or "Sample Information" section of a scientific paper.

{stats_text}

Requirements:
- Write in formal scientific English, past tense
- Include total number of runs, database split (SRA/GSA), and total data volume (in GB)
- Include species covered, tissue types, geographic locations, and collection time span
- Mention key institutions that contributed data
- Note the sequencing type (transcriptomic/genomic) and average data volume per run
- End with a note on data availability

Output ONLY the paragraph, no markdown, no headings."""

        try:
            from openai import OpenAI
            client = OpenAI(api_key=api_key, base_url="https://api.deepseek.com")
            kwargs = {"model": model_combo.currentText(), "messages": [
                {"role": "system", "content": "You are a scientific writer. Output only the requested paragraph."},
                {"role": "user", "content": prompt}
            ]}
            if "v4" in model_combo.currentText().lower():
                kwargs["temperature"] = 0.3
            else:
                kwargs["temperature"] = 0.3

            response = client.chat.completions.create(**kwargs)
            summary = response.choices[0].message.content or ""

            # ── Show result ──
            result_dlg = QDialog(self)
            result_dlg.setWindowTitle("AI Summary")
            result_dlg.setMinimumSize(600, 300)
            rl = QVBoxLayout(result_dlg)
            text_edit = QTextEdit()
            text_edit.setPlainText(summary)
            text_edit.setReadOnly(False)
            rl.addWidget(text_edit)
            rh = QHBoxLayout()
            copy_btn = QPushButton("Copy")
            def do_copy():
                from PySide6.QtWidgets import QApplication as QA
                QA.clipboard().setText(text_edit.toPlainText())
            copy_btn.clicked.connect(do_copy)
            rh.addStretch(); rh.addWidget(copy_btn)
            close_btn = QPushButton("Close"); close_btn.clicked.connect(result_dlg.accept)
            rh.addWidget(close_btn); rl.addLayout(rh)
            result_dlg.exec()
        except Exception as e:
            QMessageBox.warning(self, "AI Summary Error", str(e))

    def refresh(self):
        self._model.refresh()
        self._update_stats()
