#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
submission_gui.py — Virome Submission Pipeline 独立桌面 GUI
============================================================

功能:
  - 加载 unified_metadata.csv, 查看/编辑/批量填充/验证
  - 双击编辑, 缺失/占位符高亮
  - 批量填充常用字段 (bioproject, biosample, authors, ...)
  - 必填字段验证报告
  - 导出 CSV / TSV / Excel
  - 预览 source.src, miuvig.tsv, assembly.tsv 等产出文件

用法:
  cd virome_submission_pipeline
  python submission_gui.py [unified_metadata.csv]
"""

import sys
import os
import re

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QTableView, QHeaderView, QAbstractItemView,
    QMenuBar, QToolBar, QStatusBar, QLabel, QLineEdit, QPushButton,
    QComboBox, QGroupBox, QFormLayout, QFileDialog, QMessageBox,
    QSplitter, QSizePolicy, QTextEdit, QProgressBar,
)
from PySide6.QtCore import Qt, QAbstractTableModel, QModelIndex, Signal, QSettings, QSize
from PySide6.QtGui import QAction, QColor, QBrush, QFont


# ══════════════════════════════════════════════════════════════
# 字段定义 (与 unified_metadata.py 一致)
# ══════════════════════════════════════════════════════════════

UNIFIED_COLUMNS = [
    ("organism",               True,  "NCBI Taxonomy 物种名"),
    ("sequence_name",          True,  "FASTA 中的序列 ID"),
    ("authors",                True,  "引用作者: Last, First; ..."),
    ("collection_date",        True,  "采集日期 (YYYY-MM-DD)"),
    ("bioproject",             True,  "BioProject ID (PRJNA...)"),
    ("src-Isolate",            True,  "唯一分离株标识"),
    ("src-geo_loc_name",       True,  "采集地点 (Country:Region)"),
    ("src-Lat_Lon",            True,  "经纬度 (XX.XX N XXX.XX E)"),
    ("src-Host",               False, "宿主物种名"),
    ("src-Segment",            False, "分段病毒的片段编号"),
    ("src-Isolation-source",   False, "分离来源描述"),
    ("src-Note",               False, "额外备注"),
    ("src-Tissue_type",        False, "组织类型"),
    ("src-Collected_by",       False, "采集人"),
    ("src-Cultivar",           False, "栽培品种"),
    ("src-Dev_stage",          False, "发育阶段"),
    ("gb-sample_name",         True,  "GenBank 记录名 (≤50字符)"),
    ("gb-title",               False, "提交标题"),
    ("sra",                    False, "SRA 登录号 (SRR...)"),
    ("biosample",              True,  "BioSample ID (SAMN...)"),
    ("cmt-Assembly_Method",    True,  "组装方法"),
    ("cmt-Sequencing_Technology", True, "测序平台"),
    ("cmt-Genome_Coverage",    False, "基因组覆盖度"),
    ("cmt-Annotation_Pipeline", False, "注释流程"),
    ("bs-isolate",             False, "BioSample: 分离株标识"),
    ("bs-geo_loc_name",        False, "BioSample: 采集地点"),
    ("bs-host",                False, "BioSample: 宿主"),
    ("bs-isolation_source",    False, "BioSample: 分离来源"),
]

REQUIRED_COLS = [name for name, req, _ in UNIFIED_COLUMNS if req]
COL_DESC = {name: desc for name, _, desc in UNIFIED_COLUMNS}

NA_PLACEHOLDERS = {"NA", "N/A", "Not_Provided", "not collected",
                   "missing", "none", "unknown", "", " "}
PLACEHOLDER_RE = re.compile(
    r'XXXX|YYYY|PRJNAXXXX|Country:Region|SAMNXXXXXXXX|XX\.\d+', re.IGNORECASE)


def is_placeholder(val) -> bool:
    if val is None:
        return True
    s = str(val).strip()
    if s.lower() in NA_PLACEHOLDERS or s == "":
        return True
    return bool(PLACEHOLDER_RE.search(s))


# ══════════════════════════════════════════════════════════════
# Data Store
# ══════════════════════════════════════════════════════════════

class SubmissionStore:
    """In-memory store for unified_metadata.csv with change tracking."""

    def __init__(self):
        self._df: pd.DataFrame = pd.DataFrame(columns=[c for c, _, _ in UNIFIED_COLUMNS])
        self._filepath: str = ""
        self._modified = False

    @property
    def dataframe(self) -> pd.DataFrame:
        return self._df

    @property
    def is_loaded(self) -> bool:
        return len(self._df) > 0

    @property
    def is_modified(self) -> bool:
        return self._modified

    @property
    def row_count(self) -> int:
        return len(self._df)

    @property
    def columns(self) -> list:
        return list(self._df.columns)

    def load(self, filepath: str) -> bool:
        try:
            if filepath.endswith(".csv"):
                df = pd.read_csv(filepath, dtype=str, keep_default_na=False)
            elif filepath.endswith((".xlsx", ".xls")):
                df = pd.read_excel(filepath, dtype=str, keep_default_na=False)
            else:
                df = pd.read_csv(filepath, sep="\t", dtype=str, keep_default_na=False)
            df.fillna("", inplace=True)
            self._df = df.astype(object)
            self._filepath = filepath
            self._modified = False
            return True
        except Exception as e:
            print(f"Load error: {e}")
            return False

    def get_cell(self, row: int, col: int) -> str:
        return str(self._df.iat[row, col])

    def set_cell(self, row: int, col: int, value: str):
        self._df.iat[row, col] = value
        self._modified = True

    def get_row(self, row: int) -> dict:
        return self._df.iloc[row].to_dict()

    def set_row(self, row: int, data: dict):
        for col, val in data.items():
            if col in self._df.columns:
                self._df.iat[row, self._df.columns.get_loc(col)] = val
        self._modified = True

    def delete_rows(self, indices: list):
        self._df.drop(self._df.index[list(indices)], inplace=True)
        self._df.reset_index(drop=True, inplace=True)
        self._modified = True

    # ── Batch operations ──
    def batch_replace(self, column: str, old_value: str, new_value: str) -> int:
        if column not in self._df.columns:
            return 0
        mask = self._df[column].astype(str).str.strip() == old_value.strip()
        count = int(mask.sum())
        if count:
            self._df.loc[mask, column] = new_value
            self._modified = True
        return count

    def batch_fill_placeholder(self, column: str, new_value: str) -> int:
        if column not in self._df.columns:
            return 0
        mask = self._df[column].apply(is_placeholder)
        count = int(mask.sum())
        if count:
            self._df.loc[mask, column] = new_value
            self._modified = True
        return count

    def validate(self) -> list:
        issues = []
        for col in REQUIRED_COLS:
            if col not in self._df.columns:
                issues.append({"column": col, "missing_count": -1, "examples": []})
                continue
            mask = self._df[col].apply(is_placeholder)
            n = int(mask.sum())
            if n > 0:
                examples = self._df.loc[mask, col].unique()[:3].tolist()
                issues.append({"column": col, "missing_count": n, "examples": examples})
        return issues

    # ── Persistence ──
    def save(self, filepath: str = None) -> bool:
        path = filepath or self._filepath
        if not path:
            return False
        try:
            sep = "," if path.lower().endswith(".csv") else "\t"
            self._df.to_csv(path, sep=sep, index=False, encoding="utf-8-sig")
            self._modified = False
            self._filepath = path
            return True
        except Exception:
            return False

    def export_excel(self, filepath: str) -> bool:
        try:
            self._df.to_excel(filepath, index=False, engine="openpyxl")
            return True
        except Exception:
            return False

    def stats(self) -> dict:
        if not self.is_loaded:
            return {}
        total = len(self._df)
        filled = 0
        total_cells = total * len(self._df.columns)
        for col in self._df.columns:
            filled += (~self._df[col].apply(is_placeholder)).sum()
        return {"total": total, "filled": int(filled), "cells": total_cells,
                "pct": filled / total_cells * 100 if total_cells else 0}


# ══════════════════════════════════════════════════════════════
# Table Model
# ══════════════════════════════════════════════════════════════

class SubmissionTableModel(QAbstractTableModel):
    data_changed = Signal()

    def __init__(self, store: SubmissionStore):
        super().__init__()
        self._store = store

    def rowCount(self, parent=QModelIndex()):
        return self._store.row_count

    def columnCount(self, parent=QModelIndex()):
        return len(self._store.columns)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        r, c = index.row(), index.column()
        if r >= self._store.row_count or c >= len(self._store.columns):
            return None
        val = self._store.get_cell(r, c)

        if role in (Qt.DisplayRole, Qt.EditRole):
            return val
        elif role == Qt.BackgroundRole:
            if is_placeholder(val):
                return QBrush(QColor(255, 243, 224))  # orange tint
            return None
        elif role == Qt.ForegroundRole:
            if is_placeholder(val):
                return QBrush(QColor(230, 81, 0))
            return None
        elif role == Qt.ToolTipRole:
            col_name = self._store.columns[c]
            desc = COL_DESC.get(col_name, "")
            req = "REQUIRED" if col_name in REQUIRED_COLS else "optional"
            return f"[{col_name}] ({req}) {desc}\nValue: {val}"
        return None

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                cols = self._store.columns
                return cols[section] if section < len(cols) else ""
            return str(section + 1)
        elif role == Qt.ForegroundRole and orientation == Qt.Horizontal:
            cols = self._store.columns
            if section < len(cols) and cols[section] in REQUIRED_COLS:
                return QBrush(QColor(183, 28, 28))  # dark red for required
        elif role == Qt.FontRole and orientation == Qt.Horizontal:
            cols = self._store.columns
            if section < len(cols) and cols[section] in REQUIRED_COLS:
                f = QFont()
                f.setBold(True)
                return f
        return None

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

    def setData(self, index, value, role=Qt.EditRole):
        if role == Qt.EditRole and index.isValid():
            self._store.set_cell(index.row(), index.column(), str(value))
            self.dataChanged.emit(index, index, [Qt.DisplayRole, Qt.BackgroundRole])
            self.data_changed.emit()
            return True
        return False

    def refresh(self):
        self.beginResetModel()
        self.endResetModel()


# ══════════════════════════════════════════════════════════════
# Main Window
# ══════════════════════════════════════════════════════════════

class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Virome Submission Manager")
        self.resize(1440, 900)

        self.store = SubmissionStore()
        self._settings = QSettings("MMPV-RNA", "ViromeSubmission")

        self._setup_menu()
        self._setup_toolbar()
        self._setup_statusbar()
        self._setup_central()

        self._restore_state()

    # ── Menu ──────────────────────────────────────────
    def _setup_menu(self):
        mb = self.menuBar()

        fm = mb.addMenu("&File")
        fm.addAction("&Open...\tCtrl+O", self._on_open)
        fm.addAction("&Save\tCtrl+S", self._on_save)
        fm.addAction("Save &As...\tCtrl+Shift+S", self._on_save_as)
        fm.addSeparator()
        fm.addAction("Export &Excel...\tCtrl+E", self._on_export_excel)
        fm.addSeparator()
        fm.addAction("&Quit\tCtrl+Q", self.close)

        em = mb.addMenu("&Edit")
        em.addAction("&Find...\tCtrl+F", self._on_find)

    # ── Toolbar ───────────────────────────────────────
    def _setup_toolbar(self):
        tb = QToolBar("Main")
        tb.setIconSize(QSize(20, 20))
        tb.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.addToolBar(tb)
        tb.addAction("Open", self._on_open)
        tb.addAction("Save", self._on_save)
        tb.addSeparator()
        tb.addAction("Export Excel", self._on_export_excel)

    # ── Statusbar ─────────────────────────────────────
    def _setup_statusbar(self):
        sb = QStatusBar()
        self.setStatusBar(sb)
        self._row_label = QLabel("No data")
        sb.addPermanentWidget(self._row_label)

    def _update_status(self):
        s = self.store.stats()
        if not s:
            self._row_label.setText("No data loaded")
            return
        mod = " *" if self.store.is_modified else ""
        self._row_label.setText(
            f"{s['total']} records | {s['filled']}/{s['cells']} filled ({s['pct']:.1f}%){mod}")

    # ── Central ───────────────────────────────────────
    def _setup_central(self):
        splitter = QSplitter(Qt.Horizontal)

        # Left: table
        self._model = SubmissionTableModel(self.store)
        self._table = QTableView()
        self._table.setModel(self._model)
        self._table.setSortingEnabled(True)
        self._table.setAlternatingRowColors(True)
        self._table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self._table.setEditTriggers(
            QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed)
        self._table.verticalHeader().setDefaultSectionSize(24)

        # Context menu
        delete_act = QAction("Delete Selected Rows", self._table)
        delete_act.setShortcut("Delete")
        delete_act.triggered.connect(self._delete_selected)
        self._table.addAction(delete_act)

        # Filter bar above table
        filter_widget = QWidget()
        filter_layout = QVBoxLayout(filter_widget)
        filter_layout.setContentsMargins(0, 0, 0, 0)
        filter_layout.setSpacing(2)

        filter_row = QHBoxLayout()
        self._filter_input = QLineEdit()
        self._filter_input.setPlaceholderText("Search all fields... (Enter)")
        self._filter_input.setClearButtonEnabled(True)
        self._filter_input.returnPressed.connect(self._apply_filter)
        filter_row.addWidget(self._filter_input, 1)
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self._clear_filter)
        filter_row.addWidget(clear_btn)
        filter_layout.addLayout(filter_row)

        filter_layout.addWidget(self._table, 1)
        splitter.addWidget(filter_widget)

        # Right: tools
        tools = QTabWidget()
        tools.setMaximumWidth(400)
        tools.setMinimumWidth(260)

        tools.addTab(self._build_batch_fill_tab(), "Batch Fill")
        tools.addTab(self._build_validate_tab(), "Validate")
        tools.addTab(self._build_preview_tab(), "Preview")
        tools.addTab(self._build_info_tab(), "Column Info")

        splitter.addWidget(tools)
        splitter.setSizes([1000, 380])

        self.setCentralWidget(splitter)

    # ── Batch Fill tab ────────────────────────────────
    def _build_batch_fill_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)

        # Manual batch
        g = QGroupBox("Batch Fill")
        fl = QFormLayout(g)

        self._fill_col_combo = QComboBox()
        self._fill_col_combo.setEditable(True)
        for c, _, _ in UNIFIED_COLUMNS:
            self._fill_col_combo.addItem(c)
        fl.addRow("Column:", self._fill_col_combo)

        self._fill_old = QLineEdit()
        self._fill_old.setPlaceholderText("(empty = fill all placeholders)")
        fl.addRow("Old value:", self._fill_old)

        self._fill_new = QLineEdit()
        self._fill_new.setPlaceholderText("New value")
        fl.addRow("New value:", self._fill_new)

        btn = QPushButton("Apply")
        btn.clicked.connect(self._do_batch_fill)
        fl.addRow(btn)
        lay.addWidget(g)

        # Quick fill
        g2 = QGroupBox("Quick Fill")
        ql = QFormLayout(g2)
        for field, placeholder in [
            ("bioproject", "PRJNA..."),
            ("biosample", "SAMN..."),
            ("authors", "Last, First; ..."),
            ("collection_date", "YYYY-MM-DD"),
            ("src-geo_loc_name", "Country:Region"),
            ("src-Lat_Lon", "XX.XX N XXX.XX E"),
            ("src-Host", "Host species"),
            ("cmt-Assembly_Method", "MEGAHIT;1.2.9;..."),
            ("cmt-Sequencing_Technology", "Illumina NovaSeq 6000"),
        ]:
            row = QHBoxLayout()
            inp = QLineEdit()
            inp.setPlaceholderText(placeholder)
            b = QPushButton("Fill")
            b.setFixedWidth(50)
            b.clicked.connect(
                lambda checked, col=field, i=inp: self._quick_fill(col, i.text()))
            row.addWidget(inp, 1)
            row.addWidget(b)
            ql.addRow(field + ":", row)
        lay.addWidget(g2)
        lay.addStretch()
        return w

    # ── Validate tab ──────────────────────────────────
    def _build_validate_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)

        btn = QPushButton("Run Validation")
        btn.setStyleSheet("QPushButton { font-weight: bold; padding: 8px; }")
        btn.clicked.connect(self._do_validate)
        lay.addWidget(btn)

        self._validate_output = QTextEdit()
        self._validate_output.setReadOnly(True)
        self._validate_output.setFont(QFont("Consolas", 9))
        lay.addWidget(self._validate_output, 1)
        return w

    # ── Preview tab (source.src / miuvig / assembly) ──
    def _build_preview_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)

        lay.addWidget(QLabel("Preview generated files from submission directory:"))

        row = QHBoxLayout()
        self._preview_dir_input = QLineEdit()
        self._preview_dir_input.setPlaceholderText("Submission output directory...")
        row.addWidget(self._preview_dir_input, 1)
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self._browse_preview_dir)
        row.addWidget(browse_btn)
        load_btn = QPushButton("Load")
        load_btn.clicked.connect(self._load_preview)
        row.addWidget(load_btn)
        lay.addLayout(row)

        self._preview_combo = QComboBox()
        self._preview_combo.addItems(["source.src", "miuvig.tsv", "assembly.tsv",
                                      "biosample_template.tsv", "validation_report.txt"])
        self._preview_combo.currentTextChanged.connect(self._show_preview_file)
        lay.addWidget(self._preview_combo)

        self._preview_text = QTextEdit()
        self._preview_text.setReadOnly(True)
        self._preview_text.setFont(QFont("Consolas", 9))
        lay.addWidget(self._preview_text, 1)
        return w

    # ── Column Info tab ───────────────────────────────
    def _build_info_tab(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        self._info_label = QLabel("Load data first")
        self._info_label.setWordWrap(True)
        self._info_label.setAlignment(Qt.AlignTop)
        lay.addWidget(self._info_label, 1)
        return w

    # ── File actions ──────────────────────────────────
    def _on_open(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open unified_metadata.csv", "",
            "CSV (*.csv);;TSV (*.tsv *.txt);;Excel (*.xlsx);;All (*)")
        if path:
            if self.store.load(path):
                self._model.refresh()
                self._update_status()
                self._update_info()
                # Auto-detect submission dir
                parent = os.path.dirname(path)
                if os.path.isfile(os.path.join(parent, "source.src")):
                    self._preview_dir_input.setText(parent)
                elif os.path.isfile(os.path.join(parent, "miuvig.tsv")):
                    self._preview_dir_input.setText(parent)
            else:
                QMessageBox.warning(self, "Error", f"Failed to load: {path}")

    def _on_save(self):
        if self.store.is_modified:
            if self.store.save():
                self._update_status()
            else:
                QMessageBox.warning(self, "Error", "Save failed")
        else:
            self.statusBar().showMessage("No changes", 2000)

    def _on_save_as(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save As", "unified_metadata.csv",
            "CSV (*.csv);;TSV (*.tsv);;All (*)")
        if path:
            if self.store.save(path):
                self._update_status()
            else:
                QMessageBox.warning(self, "Error", "Save failed")

    def _on_export_excel(self):
        if not self.store.is_loaded:
            QMessageBox.information(self, "Info", "No data")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Excel", "unified_metadata_export.xlsx",
            "Excel (*.xlsx)")
        if path:
            if self.store.export_excel(path):
                self.statusBar().showMessage(f"Exported to {path}", 3000)
            else:
                QMessageBox.warning(self, "Error", "Export failed")

    # ── Edit actions ──────────────────────────────────
    def _on_find(self):
        self._filter_input.setFocus()
        self._filter_input.selectAll()

    def _apply_filter(self):
        kw = self._filter_input.text().strip().lower()
        if not kw:
            self._table.setRowHidden(0, False)
            for r in range(self.store.row_count):
                self._table.setRowHidden(r, False)
            return
        for r in range(self.store.row_count):
            row_data = self.store.get_row(r)
            match = any(kw in str(v).lower() for v in row_data.values())
            self._table.setRowHidden(r, not match)

    def _clear_filter(self):
        self._filter_input.clear()
        for r in range(self.store.row_count):
            self._table.setRowHidden(r, False)

    def _delete_selected(self):
        rows = sorted(set(idx.row() for idx in self._table.selectionModel().selectedRows()),
                      reverse=True)
        if not rows:
            return
        r = QMessageBox.question(
            self, "Delete", f"Delete {len(rows)} row(s)?",
            QMessageBox.Yes | QMessageBox.No)
        if r == QMessageBox.Yes:
            self.store.delete_rows(rows)
            self._model.refresh()
            self._update_status()

    # ── Batch fill ────────────────────────────────────
    def _do_batch_fill(self):
        col = self._fill_col_combo.currentText().strip()
        new_val = self._fill_new.text().strip()
        old_val = self._fill_old.text().strip()
        if not col or not new_val:
            QMessageBox.information(self, "Info", "Select column and enter new value")
            return
        if old_val:
            count = self.store.batch_replace(col, old_val, new_val)
        else:
            count = self.store.batch_fill_placeholder(col, new_val)
        self._model.refresh()
        self._update_status()
        QMessageBox.information(self, "Batch Fill", f"Updated {count} cells in '{col}'")

    def _quick_fill(self, column: str, value: str):
        if not value:
            return
        count = self.store.batch_fill_placeholder(column, value)
        self._model.refresh()
        self._update_status()
        self.statusBar().showMessage(f"Filled {count} cells in '{column}'", 3000)

    # ── Validate ──────────────────────────────────────
    def _do_validate(self):
        issues = self.store.validate()
        if not issues:
            self._validate_output.setPlainText("All required fields OK! Ready to submit.")
            return
        lines = [f"Validation: {len(issues)} issue(s)\n"]
        for item in issues:
            col = item["column"]
            n = item["missing_count"]
            desc = COL_DESC.get(col, "")
            if n == -1:
                lines.append(f"  [MISSING COLUMN] {col} — {desc}")
            else:
                ex = ", ".join(repr(e) for e in item["examples"][:3])
                lines.append(f"  [{n} missing] {col} — {desc}\n    e.g. {ex}")
        lines.append("\nFix all placeholders before submitting to GenBank.")
        self._validate_output.setPlainText("\n".join(lines))

    # ── Preview ───────────────────────────────────────
    def _browse_preview_dir(self):
        d = QFileDialog.getExistingDirectory(self, "Select Submission Directory")
        if d:
            self._preview_dir_input.setText(d)
            self._load_preview()

    def _load_preview(self):
        d = self._preview_dir_input.text().strip()
        if not d or not os.path.isdir(d):
            return
        # Refresh combo with actually existing files
        self._preview_combo.clear()
        for fname in ["source.src", "miuvig.tsv", "assembly.tsv",
                      "biosample_template.tsv", "validation_report.txt",
                      "unified_metadata.csv"]:
            if os.path.isfile(os.path.join(d, fname)):
                self._preview_combo.addItem(fname)
        self._show_preview_file()

    def _show_preview_file(self):
        d = self._preview_dir_input.text().strip()
        fname = self._preview_combo.currentText()
        if not d or not fname:
            self._preview_text.clear()
            return
        fpath = os.path.join(d, fname)
        if os.path.isfile(fpath):
            try:
                with open(fpath, "r", encoding="utf-8", errors="replace") as f:
                    content = f.read(50000)  # cap at 50KB
                self._preview_text.setPlainText(content)
            except Exception as e:
                self._preview_text.setPlainText(f"Error: {e}")
        else:
            self._preview_text.setPlainText(f"Not found: {fpath}")

    # ── Column info ───────────────────────────────────
    def _update_info(self):
        if not self.store.is_loaded:
            self._info_label.setText("Load data first")
            return
        lines = []
        for col in self.store.columns:
            filled = (~self.store.dataframe[col].apply(is_placeholder)).sum()
            total = len(self.store.dataframe)
            req = "REQUIRED" if col in REQUIRED_COLS else "optional"
            desc = COL_DESC.get(col, "")
            bar_len = 20
            filled_len = int(filled / total * bar_len) if total else 0
            bar = "█" * filled_len + "░" * (bar_len - filled_len)
            lines.append(
                f"<b>{col}</b> <span style='color:#999'>[{req}]</span><br>"
                f"&nbsp;&nbsp;{desc}<br>"
                f"&nbsp;&nbsp;{bar} {filled}/{total}")
        self._info_label.setText("<br>".join(lines))

    # ── State ─────────────────────────────────────────
    def _restore_state(self):
        geo = self._settings.value("geometry")
        if geo:
            self.restoreGeometry(geo)
        # Auto-load from command line arg
        if len(sys.argv) > 1:
            path = sys.argv[1]
            if os.path.isfile(path) and self.store.load(path):
                self._model.refresh()
                self._update_status()
                self._update_info()
                parent = os.path.dirname(path)
                if os.path.isfile(os.path.join(parent, "source.src")):
                    self._preview_dir_input.setText(parent)
        else:
            # Auto-detect unified_metadata.csv in common locations
            script_dir = os.path.dirname(os.path.abspath(__file__))
            for candidate in [
                os.path.join(script_dir, "submission", "unified_metadata.csv"),
                os.path.join(script_dir, "submission_metadata", "unified_metadata.csv"),
            ]:
                if os.path.isfile(candidate):
                    if self.store.load(candidate):
                        self._model.refresh()
                        self._update_status()
                        self._update_info()
                        self._preview_dir_input.setText(os.path.dirname(candidate))
                    break
            else:
                self._load_demo_data()
        self._update_status()

    def _load_demo_data(self):
        """Load demo submission data for testing."""
        import pandas as pd
        demo_rows = [
            {
                "organism": "Betacytorhabdovirus lycii",
                "sequence_name": "CRR123456_Betacytorhabdovirus_lycii_contig1",
                "authors": "Zhang, Wenda; Li, Ming; Wang, Fang",
                "collection_date": "YYYY-MM-DD",
                "bioproject": "PRJNAXXXXXX",
                "src-Isolate": "Betacytorhabdovirus_lycii_CRR123456",
                "src-geo_loc_name": "China:Ningxia",
                "src-Lat_Lon": "38.47 N 106.27 E",
                "src-Host": "Lycium barbarum",
                "src-Segment": "",
                "src-Isolation-source": "plant virome",
                "src-Note": "",
                "src-Tissue_type": "root",
                "src-Collected_by": "Ningxia University",
                "src-Cultivar": "Ningqi No.5",
                "src-Dev_stage": "2 years",
                "gb-sample_name": "Betacytorhabdovirus_lycii_CRR123456",
                "gb-title": "",
                "sra": "CRR123456",
                "biosample": "SAMNXXXXXXXX",
                "cmt-Assembly_Method": "MEGAHIT;1.2.9;default parameters",
                "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
                "cmt-Genome_Coverage": "42.5x",
                "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
                "bs-isolate": "Betacytorhabdovirus_lycii_CRR123456",
                "bs-geo_loc_name": "China:Ningxia",
                "bs-host": "Lycium barbarum",
                "bs-isolation_source": "plant virome",
            },
            {
                "organism": "Betacytorhabdovirus lycii",
                "sequence_name": "CRR123456_Betacytorhabdovirus_lycii_contig2",
                "authors": "Zhang, Wenda; Li, Ming; Wang, Fang",
                "collection_date": "YYYY-MM-DD",
                "bioproject": "PRJNAXXXXXX",
                "src-Isolate": "Betacytorhabdovirus_lycii_CRR123456",
                "src-geo_loc_name": "China:Ningxia",
                "src-Lat_Lon": "38.47 N 106.27 E",
                "src-Host": "Lycium barbarum",
                "src-Segment": "2",
                "src-Isolation-source": "plant virome",
                "src-Note": "segment 2 of multipartite virus",
                "src-Tissue_type": "root",
                "src-Collected_by": "Ningxia University",
                "src-Cultivar": "Ningqi No.5",
                "src-Dev_stage": "2 years",
                "gb-sample_name": "Betacytorhabdovirus_lycii_CRR123456_seg2",
                "gb-title": "",
                "sra": "CRR123456",
                "biosample": "SAMNXXXXXXXX",
                "cmt-Assembly_Method": "MEGAHIT;1.2.9;default parameters",
                "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
                "cmt-Genome_Coverage": "38.1x",
                "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
                "bs-isolate": "Betacytorhabdovirus_lycii_CRR123456",
                "bs-geo_loc_name": "China:Ningxia",
                "bs-host": "Lycium barbarum",
                "bs-isolation_source": "plant virome",
            },
            {
                "organism": "Torradovirus Ningxiaense",
                "sequence_name": "SRR31651831_Torradovirus_contig1",
                "authors": "Zhang, Wenda; Li, Ming; Wang, Fang",
                "collection_date": "2024-07-03",
                "bioproject": "PRJNA1218117",
                "src-Isolate": "Torradovirus_Ningxiaense_SRR31651831",
                "src-geo_loc_name": "China:Ningxia:Yinchuan",
                "src-Lat_Lon": "38.47 N 106.27 E",
                "src-Host": "Lycium barbarum",
                "src-Segment": "",
                "src-Isolation-source": "plant virome",
                "src-Note": "",
                "src-Tissue_type": "leaf",
                "src-Collected_by": "Ningxia University",
                "src-Cultivar": "Ningqi No.7",
                "src-Dev_stage": "3 years",
                "gb-sample_name": "Torradovirus_Ningxiaense_SRR31651831",
                "gb-title": "Torradovirus Ningxiaense genome sequencing",
                "sra": "SRR31651831",
                "biosample": "SAMN56789012",
                "cmt-Assembly_Method": "MEGAHIT;1.2.9;default parameters",
                "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
                "cmt-Genome_Coverage": "56.3x",
                "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
                "bs-isolate": "Torradovirus_Ningxiaense_SRR31651831",
                "bs-geo_loc_name": "China:Ningxia:Yinchuan",
                "bs-host": "Lycium barbarum",
                "bs-isolation_source": "plant virome",
            },
            {
                "organism": "Mint virus X",
                "sequence_name": "SRR33301106_MintVirusX_contig1",
                "authors": "Zhang, Wenda; Li, Ming; Wang, Fang",
                "collection_date": "2024-09-01",
                "bioproject": "PRJNA1219886",
                "src-Isolate": "MintVirusX_SRR33301106",
                "src-geo_loc_name": "China:Inner Mongolia:Alxa",
                "src-Lat_Lon": "39.08 N 105.73 E",
                "src-Host": "Lycium ruthenicum",
                "src-Segment": "",
                "src-Isolation-source": "plant virome",
                "src-Note": "",
                "src-Tissue_type": "fruit",
                "src-Collected_by": "Inner Mongolia University",
                "src-Cultivar": "wild",
                "src-Dev_stage": "3 years, mature fruit",
                "gb-sample_name": "MintVirusX_SRR33301106",
                "gb-title": "Mint virus X from Lycium ruthenicum",
                "sra": "SRR33301106",
                "biosample": "SAMN67890123",
                "cmt-Assembly_Method": "MEGAHIT;1.2.9;default parameters",
                "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
                "cmt-Genome_Coverage": "23.7x",
                "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
                "bs-isolate": "MintVirusX_SRR33301106",
                "bs-geo_loc_name": "China:Inner Mongolia:Alxa",
                "bs-host": "Lycium ruthenicum",
                "bs-isolation_source": "plant virome",
            },
            {
                "organism": "unclassified virus",
                "sequence_name": "SRR33389501_novel_virus_contig1",
                "authors": "Author, First",
                "collection_date": "2023-08-10",
                "bioproject": "PRJNAXXXXXX",
                "src-Isolate": "novel_virus_SRR33389501",
                "src-geo_loc_name": "Country:Region",
                "src-Lat_Lon": "XX.XX N XXX.XX E",
                "src-Host": "",
                "src-Segment": "",
                "src-Isolation-source": "plant virome",
                "src-Note": "novel virus, no close reference",
                "src-Tissue_type": "leaf",
                "src-Collected_by": "",
                "src-Cultivar": "",
                "src-Dev_stage": "",
                "gb-sample_name": "novel_virus_SRR33389501",
                "gb-title": "",
                "sra": "SRR33389501",
                "biosample": "SAMNXXXXXXXX",
                "cmt-Assembly_Method": "MEGAHIT;1.2.9;default parameters",
                "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
                "cmt-Genome_Coverage": "",
                "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
                "bs-isolate": "novel_virus_SRR33389501",
                "bs-geo_loc_name": "Country:Region",
                "bs-host": "",
                "bs-isolation_source": "plant virome",
            },
            {
                "organism": "Potexvirus lycii",
                "sequence_name": "SRR30124789_Potexvirus_lycii_contig1",
                "authors": "Zhang, Wenda; Li, Ming; Wang, Fang",
                "collection_date": "2024-06-15",
                "bioproject": "PRJNA1189201",
                "src-Isolate": "Potexvirus_lycii_SRR30124789",
                "src-geo_loc_name": "China:Xinjiang:Urumqi",
                "src-Lat_Lon": "43.79 N 87.58 E",
                "src-Host": "Lycium barbarum",
                "src-Segment": "",
                "src-Isolation-source": "plant virome",
                "src-Note": "",
                "src-Tissue_type": "fruit",
                "src-Collected_by": "Xinjiang University",
                "src-Cultivar": "cultivated",
                "src-Dev_stage": "5 years, ripening",
                "gb-sample_name": "Potexvirus_lycii_SRR30124789",
                "gb-title": "Potexvirus lycii from goji berry",
                "sra": "SRR30124789",
                "biosample": "SAMN78901234",
                "cmt-Assembly_Method": "MEGAHIT;1.2.9;default parameters",
                "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
                "cmt-Genome_Coverage": "61.2x",
                "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
                "bs-isolate": "Potexvirus_lycii_SRR30124789",
                "bs-geo_loc_name": "China:Xinjiang:Urumqi",
                "bs-host": "Lycium barbarum",
                "bs-isolation_source": "plant virome",
            },
        ]
        cols = [c for c, _, _ in UNIFIED_COLUMNS]
        records = []
        for row in demo_rows:
            records.append({c: row.get(c, "") for c in cols})
        self.store._df = pd.DataFrame(records, columns=cols).astype(object)
        self.store._filepath = ""
        self.store._modified = False
        self._model.refresh()
        self._update_status()
        self._update_info()
        self.statusBar().showMessage("Demo data loaded (6 virus sequences)", 5000)

    def closeEvent(self, event):
        if self.store.is_modified:
            r = QMessageBox.question(
                self, "Unsaved Changes", "Save before closing?",
                QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
            if r == QMessageBox.Save:
                if not self.store.save():
                    r2 = QMessageBox.question(
                        self, "Save Failed", "Close anyway?",
                        QMessageBox.Yes | QMessageBox.No)
                    if r2 != QMessageBox.Yes:
                        event.ignore()
                        return
            elif r == QMessageBox.Cancel:
                event.ignore()
                return
        self._settings.setValue("geometry", self.saveGeometry())
        super().closeEvent(event)


# ══════════════════════════════════════════════════════════════
# Entry
# ══════════════════════════════════════════════════════════════

def main():
    app = QApplication(sys.argv)
    app.setApplicationName("Virome Submission Manager")
    app.setOrganizationName("MMPV-RNA")
    app.setFont(QFont("Segoe UI", 10))
    app.setStyle("Fusion")

    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
