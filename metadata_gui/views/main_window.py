"""Main application window — dual pipeline: Search + Info."""

import os

from PySide6.QtWidgets import (
    QMainWindow, QTabWidget, QMenuBar, QMenu, QToolBar,
    QStatusBar, QLabel, QFileDialog, QMessageBox, QWidget,
    QVBoxLayout, QHBoxLayout,
)
from PySide6.QtCore import Qt, QSettings, QSize
from PySide6.QtGui import QAction, QKeySequence, QIcon

from models.data_store import MetadataStore


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Metadata Manager - MMPV-RNA")

        self.search_store = MetadataStore()
        self.info_store = MetadataStore()
        self._active_store = self.search_store

        self._loading = True
        self._switching = False
        self._settings = QSettings("MMPV-RNA", "MetadataManager")
        self._last_dir = ""

        self._setup_menu()
        self._setup_statusbar()
        self._setup_tabs()
        self._restore_state()

    # ── Menu ──────────────────────────────────────
    def _setup_menu(self):
        mb = self.menuBar()
        fm = mb.addMenu("&File")
        fm.addAction("&Open TSV...\tCtrl+O", self._on_open)
        fm.addAction("&Save\tCtrl+S", self._on_save)
        fm.addAction("Save &As...\tCtrl+Shift+S", self._on_save_as)
        fm.addSeparator()
        fm.addAction("&Import (CSV/Excel)...\tCtrl+I", self._on_import)
        fm.addAction("&Export Excel...\tCtrl+E", self._on_export_excel)
        fm.addSeparator()
        fm.addAction("&Quit\tCtrl+Q", self.close)
        vm = mb.addMenu("&View")
        for i, name in enumerate([
            "Search", "S-Browse", "S-Edit", "S-Viz",
            "Info", "I-Browse", "I-Edit", "I-Viz",
        ]):
            vm.addAction(f"&{name}", lambda idx=i: self._tabs.setCurrentIndex(idx))
        mb.addMenu("&Help").addAction("&About", self._show_about)

    # ── Toolbar ───────────────────────────────────
    # ── Statusbar ─────────────────────────────────
    def _setup_statusbar(self):
        self._status = QStatusBar()
        self.setStatusBar(self._status)
        self._row_label = QLabel("")
        self._status.addPermanentWidget(self._row_label)

    def set_status(self, text: str, timeout: int = 0):
        self._status.showMessage(text, timeout)

    def _update_row_label(self):
        s = self.search_store.row_count
        i = self.info_store.row_count
        act = "Search" if self._active_store is self.search_store else "Info"
        self._row_label.setText(
            f"Search: {s} rows  |  Info: {i} rows  |  Active: {act}")

    # ── Tabs ──────────────────────────────────────
    def _setup_tabs(self):
        self._tabs = QTabWidget()
        self._tabs.setDocumentMode(True)
        self.setCentralWidget(self._tabs)

        # ── Search pipeline: 0-3 ──
        self._tabs.addTab(self._make_placeholder("Search", "Enter species name to search SRA/GSA."),
                          "Search")
        self._tabs.addTab(self._make_placeholder("S-Browse", "Search results will appear here."),
                          "S-Browse")
        self._tabs.addTab(self._make_placeholder("S-Edit", "Select a row in S-Browse to edit."),
                          "S-Edit")
        self._tabs.addTab(self._make_placeholder("S-Viz", "Visualize search results."),
                          "S-Viz")

        # ── Info pipeline: 4-7 ──
        self._tabs.addTab(self._make_placeholder("Info", "Paste Run IDs or import from Search."),
                          "Info")
        self._tabs.addTab(self._make_placeholder("I-Browse", "Info extraction results appear here."),
                          "I-Browse")
        self._tabs.addTab(self._make_placeholder("I-Edit", "Select a row in I-Browse to edit."),
                          "I-Edit")
        self._tabs.addTab(self._make_placeholder("I-Viz", "Visualize info extraction results."),
                          "I-Viz")

        self._tabs.currentChanged.connect(self._on_tab_changed)

    def _make_placeholder(self, title, hint):
        w = QWidget()
        l = QVBoxLayout(w)
        lb = QLabel(f"<h2>{title}</h2><p>{hint}</p>")
        lb.setAlignment(Qt.AlignCenter)
        l.addWidget(lb)
        return w

    # ── Tab switch ────────────────────────────────
    def _on_tab_changed(self, index):
        if self._loading or self._switching:
            return

        # Search pipeline
        if index == 0:   # Search
            self._active_store = self.search_store
            self._activate_search()
        elif index == 1:  # S-Browse
            self._active_store = self.search_store
            self._activate_browse(1, self.search_store, "S-Browse")
        elif index == 2:  # S-Edit
            self._active_store = self.search_store
            self._activate_edit(2, self.search_store, "S-Edit")
        elif index == 3:  # S-Viz
            self._active_store = self.search_store
            self._activate_viz(3, self.search_store, "S-Viz")

        # Info pipeline
        elif index == 4:  # Info
            self._active_store = self.info_store
            self._activate_info()
        elif index == 5:  # I-Browse
            self._active_store = self.info_store
            self._activate_browse(5, self.info_store, "I-Browse")
        elif index == 6:  # I-Edit
            self._active_store = self.info_store
            self._activate_edit(6, self.info_store, "I-Edit")
        elif index == 7:  # I-Viz
            self._active_store = self.info_store
            self._activate_viz(7, self.info_store, "I-Viz")

        self._update_row_label()

    # ── Tab activators ────────────────────────────
    def _activate_search(self):
        from views.search_view import SearchPanel
        # Don't recreate if already active
        w = self._tabs.widget(0)
        if isinstance(w, SearchPanel):
            return
        if self.search_store.is_loaded or True:  # always allow search
            self._switching = True
            self._tabs.blockSignals(True)
            try:
                sp = SearchPanel(self.search_store)
                sp.import_requested.connect(self._on_search_import)
                sp.send_to_info.connect(self._on_send_to_info)
                self._tabs.removeTab(0)
                self._tabs.insertTab(0, sp, "Search")
                self._tabs.setCurrentIndex(0)
            finally:
                self._tabs.blockSignals(False)
                self._switching = False

    def _activate_info(self):
        from views.info_view import InfoPanel
        w = self._tabs.widget(4)
        if isinstance(w, InfoPanel):
            return
        self._switching = True
        self._tabs.blockSignals(True)
        try:
            ip = InfoPanel(self.info_store, self.search_store)
            ip.import_requested.connect(self._on_info_import)
            self._tabs.removeTab(4)
            self._tabs.insertTab(4, ip, "Info")
            self._tabs.setCurrentIndex(4)
        finally:
            self._tabs.blockSignals(False)
            self._switching = False

    def _activate_browse(self, tab_idx, store, label):
        from views.metadata_table import MetadataTableView
        w = self._tabs.widget(tab_idx)
        if isinstance(w, MetadataTableView):
            return
        if store.is_loaded:
            self._switching = True
            self._tabs.blockSignals(True)
            try:
                bt = MetadataTableView(store)
                self._tabs.removeTab(tab_idx)
                self._tabs.insertTab(tab_idx, bt, label)
                self._tabs.setCurrentIndex(tab_idx)
            finally:
                self._tabs.blockSignals(False)
                self._switching = False

    def _activate_edit(self, tab_idx, store, label):
        from views.detail_panel import DetailEditPanel
        w = self._tabs.widget(tab_idx)
        if isinstance(w, DetailEditPanel):
            return
        if store.is_loaded:
            self._switching = True
            self._tabs.blockSignals(True)
            try:
                dp = DetailEditPanel(store)
                self._tabs.removeTab(tab_idx)
                self._tabs.insertTab(tab_idx, dp, label)
                self._tabs.setCurrentIndex(tab_idx)
            finally:
                self._tabs.blockSignals(False)
                self._switching = False

    def _activate_viz(self, tab_idx, store, label):
        from views.viz_panel import VisualizationPanel
        w = self._tabs.widget(tab_idx)
        if isinstance(w, VisualizationPanel):
            return
        if store.is_loaded:
            self._switching = True
            self._tabs.blockSignals(True)
            try:
                vp = VisualizationPanel(store)
                self._tabs.removeTab(tab_idx)
                self._tabs.insertTab(tab_idx, vp, label)
                self._tabs.setCurrentIndex(tab_idx)
            finally:
                self._tabs.blockSignals(False)
                self._switching = False

    # ── Import handlers ───────────────────────────
    def _on_search_import(self, df):
        self._replace_store(self.search_store, df)
        self._refresh_pipeline_tabs(0, self.search_store)
        self._update_row_label()

    def _on_info_import(self, df):
        self._replace_store(self.info_store, df)
        self._refresh_pipeline_tabs(4, self.info_store)
        self._update_row_label()

    def _on_send_to_info(self, run_ids):
        """Receive Run IDs from Search tab, send to Info tab input."""
        # Info tab handles this internally via its store reference
        self._tabs.setCurrentIndex(4)  # Switch to Info tab

    @staticmethod
    def _replace_store(store, df):
        import pandas as pd
        cols = list(store.dataframe.columns)
        records = [{c: str(row[c]) if c in df.columns and pd.notna(row[c]) else ""
                    for c in cols} for _, row in df.iterrows()]
        store._df = pd.DataFrame(records, columns=cols)
        store._modified = True
        store._invalidate_cache()
        store._imported.clear()
        store._ai_filled.clear()

    def _refresh_pipeline_tabs(self, base_idx, store):
        """Reset Browse/Edit/Viz tabs for a pipeline."""
        self._loading = True
        self._tabs.blockSignals(True)
        try:
            for i in range(base_idx + 1, base_idx + 4):
                self._tabs.removeTab(base_idx + 1)
            self._tabs.insertTab(base_idx + 1,
                self._make_placeholder("Browse", "Click to load data."),
                f"S-Browse" if base_idx == 0 else "I-Browse")
            self._tabs.insertTab(base_idx + 2,
                self._make_placeholder("Edit", "Select a row to edit."),
                f"S-Edit" if base_idx == 0 else "I-Edit")
            self._tabs.insertTab(base_idx + 3,
                self._make_placeholder("Viz", "Load data for charts."),
                f"S-Viz" if base_idx == 0 else "I-Viz")
        finally:
            self._tabs.blockSignals(False)
            self._loading = False

    # ── File actions ──────────────────────────────
    def _on_open(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open TSV", "", "TSV Files (*.tsv *.txt);;All Files (*)")
        if path and self._active_store.load(path):
            self._active_store._imported.clear()
            self._active_store._ai_filled.clear()
            self.set_status(f"Loaded {self._active_store.row_count} rows")
            self._update_row_label()
            base = 0 if self._active_store is self.search_store else 4
            self._refresh_pipeline_tabs(base, self._active_store)

    def _on_save(self):
        s = self._active_store
        if s.is_modified:
            s.save()
            self.set_status("Saved", 3000)

    def _on_save_as(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save As", "", "TSV (*.tsv)")
        if path and self._active_store.save(path):
            self.set_status(f"Saved to {path}", 3000)

    def _on_import(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Import", "", "Supported (*.tsv *.csv *.xlsx *.xls);;All (*)")
        if path and self._active_store.import_file(path):
            self._active_store._imported.clear()
            self._active_store._ai_filled.clear()
            self.set_status(f"Imported {self._active_store.row_count} records")
            self._update_row_label()
            base = 0 if self._active_store is self.search_store else 4
            self._refresh_pipeline_tabs(base, self._active_store)

    def _on_export_excel(self):
        if not self._active_store.is_loaded:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Excel", "export.xlsx", "Excel (*.xlsx)")
        if path:
            self._active_store.export_excel(path)

    def _show_about(self):
        QMessageBox.about(self, "About",
            "<h3>Metadata Manager v2.0</h3>"
            "<p>Dual pipeline: Search (gsa_sra.search.py) + Info (gsa_sra.info.py).</p>")

    # ── State ─────────────────────────────────────
    def _restore_state(self):
        geo = self._settings.value("geometry")
        if geo: self.restoreGeometry(geo)
        state = self._settings.value("windowState")
        if state: self.restoreState(state)

        # Load demo data for both pipelines (same 13-column format)
        self._load_search_demo()
        self._load_info_demo()

        self._active_store = self.search_store
        self._update_row_label()
        # Auto-activate Search panel BEFORE releasing loading lock
        self._activate_search()
        self._loading = False

    def _load_search_demo(self):
        """Search demo: mimics REAL gsa_sra.search.py detailed output (13 cols)."""
        import pandas as pd
        demos = pd.DataFrame([
            {"Database":"SRA","Run":"SRR31651831","BioProject":"PRJNA1218117",
             "BioSample":"SAMN45678901","ScientificName":"Lycium barbarum",
             "Tissue":"root","Age_GrowthStage":"2 years",
             "Location":"China: Ningxia, Yinchuan","LibraryStrategy":"RNA-Seq",
             "LibrarySource":"TRANSCRIPTOMIC","Platform":"ILLUMINA",
             "CenterName":"Ningxia University","ReleaseDate":"2025-06-13"},
            {"Database":"SRA","Run":"SRR33301106","BioProject":"PRJNA1219886",
             "BioSample":"SAMN45678902","ScientificName":"Lycium ruthenicum",
             "Tissue":"fruit","Age_GrowthStage":"3 years | mature fruit stage",
             "Location":"China: Inner Mongolia, Alxa","LibraryStrategy":"RNA-Seq",
             "LibrarySource":"TRANSCRIPTOMIC","Platform":"ILLUMINA",
             "CenterName":"Inner Mongolia University","ReleaseDate":"2025-11-20"},
            {"Database":"SRA","Run":"SRR33389501","BioProject":"PRJNA1236100",
             "BioSample":"SAMN45678903","ScientificName":"Lycium chinense",
             "Tissue":"leaf","Age_GrowthStage":"",
             "Location":"China: Qinghai, Xining","LibraryStrategy":"RNA-Seq",
             "LibrarySource":"TRANSCRIPTOMIC","Platform":"ILLUMINA",
             "CenterName":"Qinghai University","ReleaseDate":"2025-12-15"},
            {"Database":"SRA","Run":"SRR29563412","BioProject":"PRJNA1146557",
             "BioSample":"SAMN45678904","ScientificName":"Lycium barbarum",
             "Tissue":"leaf","Age_GrowthStage":"",
             "Location":"China: Gansu, Lanzhou","LibraryStrategy":"RNA-Seq",
             "LibrarySource":"TRANSCRIPTOMIC","Platform":"ILLUMINA",
             "CenterName":"Gansu Agricultural University","ReleaseDate":"2025-04-12"},
            {"Database":"SRA","Run":"ERR13512478","BioProject":"PRJEB78201",
             "BioSample":"SAME7890123","ScientificName":"Lycium barbarum",
             "Tissue":"leaf","Age_GrowthStage":"1 year",
             "Location":"United Kingdom: England, London","LibraryStrategy":"WGS",
             "LibrarySource":"GENOMIC","Platform":"ILLUMINA",
             "CenterName":"Royal Botanic Gardens Kew","ReleaseDate":"2025-03-20"},
            {"Database":"SRA","Run":"DRR512340","BioProject":"PRJDB16890",
             "BioSample":"SAMD01234567","ScientificName":"Lycium chinense",
             "Tissue":"flower","Age_GrowthStage":"2 years | flowering stage",
             "Location":"Japan: Tokyo, Hachioji","LibraryStrategy":"RNA-Seq",
             "LibrarySource":"TRANSCRIPTOMIC","Platform":"ILLUMINA",
             "CenterName":"University of Tokyo","ReleaseDate":"2025-01-10"},
            {"Database":"GSA","Run":"CRR1126132","BioProject":"PRJCA025572",
             "BioSample":"SAMC3551095","ScientificName":"Lycium barbarum",
             "Tissue":"leaf","Age_GrowthStage":"3 year",
             "Location":"","LibraryStrategy":"RNA-Seq",
             "LibrarySource":"TRANSCRIPTOMIC","Platform":"ILLUMINA",
             "CenterName":"Beijing Forestry University","ReleaseDate":"2024-04-25"},
            {"Database":"GSA","Run":"CRR1126144","BioProject":"PRJCA025572",
             "BioSample":"SAMC3551107","ScientificName":"Lycium barbarum",
             "Tissue":"pistil","Age_GrowthStage":"3 year",
             "Location":"","LibraryStrategy":"RNA-Seq",
             "LibrarySource":"TRANSCRIPTOMIC","Platform":"ILLUMINA",
             "CenterName":"Beijing Forestry University","ReleaseDate":"2024-04-25"},
            {"Database":"GSA","Run":"CRR1128966","BioProject":"PRJCA025589",
             "BioSample":"SAMC3552001","ScientificName":"Lycium barbarum",
             "Tissue":"anther","Age_GrowthStage":"Archeocyte stage",
             "Location":"","LibraryStrategy":"RNA-Seq",
             "LibrarySource":"TRANSCRIPTOMIC","Platform":"ILLUMINA",
             "CenterName":"North Minzu University","ReleaseDate":"2026-04-28"},
            {"Database":"SRA","Run":"SRR32890123","BioProject":"PRJNA1256789",
             "BioSample":"SAMN45678905","ScientificName":"Lycium barbarum",
             "Tissue":"fruit","Age_GrowthStage":"",
             "Location":"USA: Maryland, Beltsville","LibraryStrategy":"RNA-Seq",
             "LibrarySource":"TRANSCRIPTOMIC","Platform":"ILLUMINA",
             "CenterName":"USDA-ARS","ReleaseDate":"2025-07-25"},
        ])
        self.search_store._df = demos.astype(object)
        # Add data volume columns
        vol_data = {"FileSize_MB": ["156", "210", "89", "134", "298", "67", "45", "45", "52", "178"],
                     "Bases": ["4.2G", "5.8G", "2.1G", "3.5G", "7.9G", "1.8G", "1.2G", "1.2G", "1.4G", "4.5G"],
                     "Spots": ["14.2M", "19.3M", "7.1M", "11.8M", "26.3M", "6.0M", "4.0M", "4.0M", "4.7M", "15.0M"]}
        for col, vals in vol_data.items():
            self.search_store._df[col] = vals
        self.search_store._modified = False

    def _load_info_demo(self):
        """Info demo: mimics gsa_sra.info.py output — AI-filled, rich metadata."""
        import pandas as pd
        demos = pd.DataFrame([
            {"Run":"SRR31651831","ReleaseDate":"2025-06-13",
             "CollectionDate":"2024-07-03","Location":"China, Ningxia, Yinchuan_AI",
             "Source":"Ningqi No.5","Tissue":"root","Age_GrowthStage":"2 years",
             "ScientificName":"Lycium barbarum","TaxID":"112863",
             "LibrarySource":"TRANSCRIPTOMIC","CenterName":"Ningxia University",
             "BioProject":"PRJNA1218117","BioSample":"SAMN45678901","PMID":"40401760"},
            {"Run":"SRR33301106","ReleaseDate":"2025-11-20",
             "CollectionDate":"2024-09-01","Location":"China, Inner Mongolia, Hohhot_AI",
             "Source":"Ningqi No.7","Tissue":"fruit",
             "Age_GrowthStage":"3 years | mature fruit stage",
             "ScientificName":"Lycium ruthenicum","TaxID":"112864",
             "LibrarySource":"TRANSCRIPTOMIC","CenterName":"Inner Mongolia University",
             "BioProject":"PRJNA1219886","BioSample":"SAMN45678902","PMID":"40439086"},
            {"Run":"CRR1126132","ReleaseDate":"2024-04-25",
             "CollectionDate":"2024-04-24","Location":"China, Beijing, Beijing_AI",
             "Source":"Ningqi No.1","Tissue":"leaf","Age_GrowthStage":"3 years",
             "ScientificName":"Lycium barbarum","TaxID":"112863",
             "LibrarySource":"TRANSCRIPTOMIC","CenterName":"Beijing Forestry University",
             "BioProject":"PRJCA025572","BioSample":"SAMC3551095","PMID":""},
            {"Run":"CRR1128966","ReleaseDate":"2026-04-28",
             "CollectionDate":"2024-04-28","Location":"China, Ningxia, Yinchuan_AI",
             "Source":"Ningqi No.1","Tissue":"anther",
             "Age_GrowthStage":"archeocyte stage",
             "ScientificName":"Lycium barbarum","TaxID":"112863",
             "LibrarySource":"TRANSCRIPTOMIC","CenterName":"North Minzu University",
             "BioProject":"PRJCA025589","BioSample":"SAMC3552001","PMID":""},
            {"Run":"SRR29563412","ReleaseDate":"2025-04-12",
             "CollectionDate":"2024-05-20","Location":"China, Gansu, Lanzhou_AI",
             "Source":"wild","Tissue":"leaf","Age_GrowthStage":"",
             "ScientificName":"Lycium barbarum","TaxID":"112863",
             "LibrarySource":"TRANSCRIPTOMIC","CenterName":"Gansu Agricultural University",
             "BioProject":"PRJNA1146557","BioSample":"SAMN45678904","PMID":"39867230"},
            {"Run":"ERR13512478","ReleaseDate":"2025-03-20",
             "CollectionDate":"2023-11-01","Location":"United Kingdom, England, London_AI",
             "Source":"cultivated","Tissue":"leaf","Age_GrowthStage":"1 year",
             "ScientificName":"Lycium barbarum","TaxID":"112863",
             "LibrarySource":"GENOMIC","CenterName":"Royal Botanic Gardens Kew",
             "BioProject":"PRJEB78201","BioSample":"SAME7890123","PMID":"39501234"},
            {"Run":"DRR512340","ReleaseDate":"2025-01-10",
             "CollectionDate":"2023-07-15","Location":"Japan, Tokyo, Hachioji_AI",
             "Source":"cultivated","Tissue":"flower",
             "Age_GrowthStage":"2 years | flowering stage",
             "ScientificName":"Lycium chinense","TaxID":"112883",
             "LibrarySource":"TRANSCRIPTOMIC","CenterName":"University of Tokyo",
             "BioProject":"PRJDB16890","BioSample":"SAMD01234567","PMID":"40123456"},
            {"Run":"SRR32890123","ReleaseDate":"2025-07-25",
             "CollectionDate":"2024-04-10","Location":"USA, Maryland, Beltsville_AI",
             "Source":"wild","Tissue":"fruit","Age_GrowthStage":"",
             "ScientificName":"Lycium barbarum","TaxID":"112863",
             "LibrarySource":"TRANSCRIPTOMIC","CenterName":"USDA-ARS",
             "BioProject":"PRJNA1256789","BioSample":"SAMN45678905","PMID":""},
        ])
        self.info_store._df = demos.astype(object)
        self.info_store._modified = False

    def closeEvent(self, event):
        if self.search_store.is_modified or self.info_store.is_modified:
            r = QMessageBox.question(self, "Unsaved Changes",
                "Save changes before closing?",
                QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
            if r == QMessageBox.Save:
                for s in [self.search_store, self.info_store]:
                    if s.is_modified and not s.save():
                        QMessageBox.warning(self, "Error", "Save failed.")
            elif r == QMessageBox.Cancel:
                event.ignore()
                return
        self._settings.setValue("geometry", self.saveGeometry())
        self._settings.setValue("windowState", self.saveState())
        super().closeEvent(event)
