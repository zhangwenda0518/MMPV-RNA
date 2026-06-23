"""Detail / Edit panel — full record view with inline editing + links."""

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLineEdit,
    QTextEdit, QPushButton, QLabel, QScrollArea, QGroupBox,
    QMessageBox, QSplitter, QSizePolicy,
)
from PySide6.QtCore import Qt, Signal, QUrl
from PySide6.QtGui import QFont, QDesktopServices

from models.data_store import MetadataStore
from views.metadata_table import LINK_COLUMNS, _build_url


class DetailEditPanel(QWidget):
    """Detail/Edit tab — row-level metadata viewer and editor."""

    def __init__(self, store: MetadataStore):
        super().__init__()
        self._store = store
        self._current_row = -1
        self._editors = {}
        self._setup_ui()

    def _setup_ui(self):
        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        splitter = QSplitter(Qt.Horizontal)

        # ── Left: Record selector ──────────────────
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)

        nav_label = QLabel("<b>Record Navigator</b>")
        left_layout.addWidget(nav_label)

        nav_layout = QHBoxLayout()
        self._prev_btn = QPushButton("< Prev")
        self._prev_btn.clicked.connect(self._go_prev)
        nav_layout.addWidget(self._prev_btn)

        self._row_label = QLabel("Row: -")
        self._row_label.setAlignment(Qt.AlignCenter)
        nav_layout.addWidget(self._row_label)

        self._next_btn = QPushButton("Next >")
        self._next_btn.clicked.connect(self._go_next)
        nav_layout.addWidget(self._next_btn)

        left_layout.addLayout(nav_layout)

        self._run_list_label = QLabel()
        self._run_list_label.setWordWrap(True)
        left_layout.addWidget(self._run_list_label)

        left_layout.addStretch()

        splitter.addWidget(left)

        # ── Right: Form editor ─────────────────────
        right_scroll = QScrollArea()
        right_scroll.setWidgetResizable(True)
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)

        form_group = QGroupBox("Record Details")
        self._form_layout = QFormLayout(form_group)
        self._form_layout.setSpacing(8)
        right_layout.addWidget(form_group)

        btn_layout = QHBoxLayout()
        self._save_btn = QPushButton("Save Changes")
        self._save_btn.clicked.connect(self._save)
        self._save_btn.setEnabled(False)
        btn_layout.addWidget(self._save_btn)

        self._reset_btn = QPushButton("Reset")
        self._reset_btn.clicked.connect(self._reset)
        self._reset_btn.setEnabled(False)
        btn_layout.addWidget(self._reset_btn)

        right_layout.addLayout(btn_layout)

        # JSON cache section
        json_group = QGroupBox("Raw Cache Data (if available)")
        json_layout = QVBoxLayout(json_group)
        self._json_view = QTextEdit()
        self._json_view.setReadOnly(True)
        self._json_view.setFont(QFont("Consolas", 9))
        self._json_view.setMaximumHeight(200)
        json_layout.addWidget(self._json_view)
        right_layout.addWidget(json_group)

        right_layout.addStretch()
        right_scroll.setWidget(right_widget)
        splitter.addWidget(right_scroll)

        splitter.setSizes([250, 750])
        layout.addWidget(splitter)

    def load_record(self, row: int):
        """Load a row into the form."""
        if row < 0 or row >= self._store.row_count:
            return
        self._current_row = row
        data = self._store.get_row(row)

        # Clear old editors
        self._clear_form()

        for col, val in data.items():
            val_str = str(val) if not self._store.is_missing(val) else ""

            if col in LINK_COLUMNS and val_str:
                # Link field: editor + open button
                row_widget = QWidget()
                row_layout = QHBoxLayout(row_widget)
                row_layout.setContentsMargins(0, 0, 0, 0)
                row_layout.setSpacing(4)

                editor = QLineEdit()
                editor.setText(val_str)
                editor.textChanged.connect(self._on_edit)
                row_layout.addWidget(editor, 1)

                open_btn = QPushButton("Open")
                open_btn.setFixedWidth(50)
                url = _build_url(col, val_str)
                if url:
                    open_btn.clicked.connect(
                        lambda checked, u=url: QDesktopServices.openUrl(QUrl(u)))
                    open_btn.setStyleSheet(
                        "QPushButton { color: #0072B2; font-weight: bold; }")
                else:
                    open_btn.setEnabled(False)
                row_layout.addWidget(open_btn)

                self._form_layout.addRow(col, row_widget)
                self._editors[col] = editor
            else:
                editor = QLineEdit()
                editor.setText(val_str)
                editor.textChanged.connect(self._on_edit)
                if self._store.is_ai_filled(row, col):
                    editor.setStyleSheet("background-color: #dcffdc;")
                elif self._store.is_missing(val):
                    editor.setStyleSheet("background-color: #ffffc8;")
                    editor.setPlaceholderText(f"(missing) {col}")
                self._form_layout.addRow(col, editor)
                self._editors[col] = editor

        self._row_label.setText(f"Row: {row + 1} / {self._store.row_count}")
        run_val = data.get("Run", data.get("Run", ""))
        sci_name = data.get("ScientificName", "")
        self._run_list_label.setText(
            f"<b>Run:</b> {run_val}<br><b>Species:</b> {sci_name}")

        self._save_btn.setEnabled(False)
        self._reset_btn.setEnabled(False)

        # Try to load JSON cache
        self._load_json_cache(run_val)

    def _load_json_cache(self, run_id: str):
        """Try to load cached JSON metadata for this run."""
        import os
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cache_dir = os.path.join(
            base, "public_metadata_pipeline",
            "public_data_pipeline_output", "info",
            "GSA_Results", "0_web_cache")
        if not os.path.isdir(cache_dir):
            cache_dir = os.path.join(base, "..", "public_metadata_pipeline",
                                     "public_data_pipeline_output", "info",
                                     "GSA_Results", "0_web_cache")

        json_path = os.path.join(cache_dir, f"{run_id}.json")
        if os.path.isfile(json_path):
            try:
                import json
                with open(json_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                self._json_view.setPlainText(
                    json.dumps(data, indent=2, ensure_ascii=False))
            except Exception:
                self._json_view.clear()
        else:
            self._json_view.clear()

    def _clear_form(self):
        while self._form_layout.rowCount() > 0:
            self._form_layout.removeRow(0)
        self._editors.clear()

    def _on_edit(self):
        self._save_btn.setEnabled(True)
        self._reset_btn.setEnabled(True)

    def _save(self):
        if self._current_row < 0:
            return
        data = {col: e.text() for col, e in self._editors.items()}
        self._store.set_row(self._current_row, data)
        self._save_btn.setEnabled(False)
        self._reset_btn.setEnabled(False)

        # Re-highlight
        for col, editor in self._editors.items():
            if self._store.is_missing(editor.text()):
                editor.setStyleSheet("background-color: #ffffc8;")
            else:
                editor.setStyleSheet("")

    def _reset(self):
        if self._current_row >= 0:
            self.load_record(self._current_row)

    def _go_prev(self):
        if self._current_row > 0:
            self.load_record(self._current_row - 1)

    def _go_next(self):
        if self._current_row < self._store.row_count - 1:
            self.load_record(self._current_row + 1)
