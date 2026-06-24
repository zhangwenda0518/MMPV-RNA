"""Visualization panel — 2-col scrollable charts + export + SCI figure."""

import os
import tempfile
from collections import Counter

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QComboBox, QGroupBox, QScrollArea, QSizePolicy,
    QProgressBar, QFileDialog, QFrame, QGridLayout,
)
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QPixmap

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np

from models.data_store import MetadataStore

FIG_W = 7
FIG_DPI = 96
BAR_ITEM_H = 0.25
BAR_MIN_H = 2.5
FULL_W = 12
HEAT_ITEM_H = 0.35
TIME_H = 3.0


class SCIPlotWorker(QThread):
    finished = Signal(str, str)
    error = Signal(str)
    progress = Signal(str)

    def __init__(self, csv_path: str, outdir: str):
        super().__init__()
        self._csv, self._outdir = csv_path, outdir

    def run(self):
        import subprocess, sys
        gui_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        project_dir = os.path.dirname(gui_dir)
        script = os.path.join(project_dir, "public_metadata_pipeline", "gsa_sra.plot.py")
        self.progress.emit("Generating SCI figure...")
        try:
            r = subprocess.run([sys.executable, script, "-i", self._csv, "-o", self._outdir],
                               capture_output=True, text=True, timeout=120,
                               encoding="utf-8", errors="replace",
                               cwd=os.path.dirname(script))
            if r.returncode != 0:
                self.error.emit(r.stderr[-500:] or "Unknown error"); return
            png = os.path.join(self._outdir, "Combined_Landscape_Full.png")
            pdf = os.path.join(self._outdir, "Combined_Landscape_Full.pdf")
            self.finished.emit(png, pdf) if os.path.isfile(png) else self.error.emit(f"Missing: {png}")
        except Exception as e:
            self.error.emit(str(e))


class ZoomablePreview(QScrollArea):
    def __init__(self):
        super().__init__()
        self._label = QLabel(); self._label.setAlignment(Qt.AlignCenter)
        self.setWidget(self._label); self.setWidgetResizable(False)
        self._scale = 0.5; self._pixmap = None

    def set_image(self, path: str):
        self._pixmap = QPixmap(path); self._apply_scale()

    def zoom_in(self): self._scale = min(2.0, self._scale+0.12); self._apply_scale()
    def zoom_out(self): self._scale = max(0.12, self._scale-0.12); self._apply_scale()
    def fit_width(self):
        if self._pixmap and not self._pixmap.isNull():
            self._scale = (self.viewport().width()-30)/self._pixmap.width(); self._apply_scale()

    def _apply_scale(self):
        if self._pixmap and not self._pixmap.isNull():
            w, h = max(1,int(self._pixmap.width()*self._scale)), max(1,int(self._pixmap.height()*self._scale))
            pm = self._pixmap.scaled(w, h, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            self._label.setPixmap(pm); self._label.resize(pm.size())


class VisualizationPanel(QWidget):
    def __init__(self, store: MetadataStore):
        super().__init__()
        self._store = store; self._worker = None; self._temp_files = []
        self._setup_ui(); self._plot_all()

    def _setup_ui(self):
        main = QVBoxLayout(self); main.setContentsMargins(0, 0, 0, 0)

        # ── Control bar ──
        ctrl = QHBoxLayout(); ctrl.setContentsMargins(8, 8, 8, 4)
        ctrl.addWidget(QLabel("<b>Visualization</b>"))
        self._field_combo = QComboBox()
        for c in self._store.columns: self._field_combo.addItem(c)
        self._field_combo.setCurrentText("Tissue")
        ctrl.addWidget(QLabel("Field:")); ctrl.addWidget(self._field_combo)
        refresh_btn = QPushButton("Refresh Charts"); refresh_btn.clicked.connect(self._plot_all)
        ctrl.addWidget(refresh_btn)
        ctrl.addSpacing(16)
        self._sci_btn = QPushButton("Generate SCI Figure")
        self._sci_btn.setStyleSheet("QPushButton{font-weight:bold;background:#0072B2;color:white;padding:6px 16px;border-radius:4px}")
        self._sci_btn.clicked.connect(self._generate_sci_figure); ctrl.addWidget(self._sci_btn)
        ctrl.addStretch(); main.addLayout(ctrl)

        # ── SCI progress ──
        self._sci_progress = QProgressBar(); self._sci_progress.setVisible(False)
        self._sci_progress.setRange(0,0); self._sci_progress.setMaximumHeight(3)
        self._sci_progress.setTextVisible(False); main.addWidget(self._sci_progress)

        # ── SCI preview ──
        self._sci_frame = QFrame(); self._sci_frame.setVisible(False); self._sci_frame.setFrameShape(QFrame.StyledPanel)
        sf = QVBoxLayout(self._sci_frame); sf.setContentsMargins(4,4,4,4)
        zc = QHBoxLayout(); self._sci_status = QLabel(""); zc.addWidget(self._sci_status); zc.addStretch()
        zc.addWidget(QLabel("Zoom:"))
        for txt, fn in [("-",lambda:self._preview.zoom_out()),("+",lambda:self._preview.zoom_in()),("Fit",lambda:self._preview.fit_width())]:
            b=QPushButton(txt); b.setFixedWidth(36); b.clicked.connect(fn); zc.addWidget(b)
        for txt, fn in [("Save PNG",self._save_png),("Save PDF",self._save_pdf)]:
            b=QPushButton(txt); b.clicked.connect(fn); zc.addWidget(b)
        sf.addLayout(zc)
        self._preview = ZoomablePreview(); self._preview.setMinimumHeight(250)
        sf.addWidget(self._preview, 1); main.addWidget(self._sci_frame)

        # ── Separator ──
        sep = QFrame(); sep.setFrameShape(QFrame.HLine); sep.setStyleSheet("QFrame{color:#ccc}"); main.addWidget(sep)

        # ── Charts scroll area ──
        scroll = QScrollArea(); scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self._chart_widget = QWidget()
        self._chart_layout = QVBoxLayout(self._chart_widget)
        self._chart_layout.setContentsMargins(4, 8, 4, 24); self._chart_layout.setSpacing(12)
        scroll.setWidget(self._chart_widget); main.addWidget(scroll, 1)

    # ── SCI figure ─────────────────────────────
    def _generate_sci_figure(self):
        df = self._store.dataframe
        if df is None or len(df)==0: return
        fd, csv = tempfile.mkstemp(suffix=".csv"); os.close(fd); self._temp_files.append(csv)
        df.to_csv(csv, index=False)
        outdir = tempfile.mkdtemp(prefix="sci_"); self._temp_files.append(outdir)
        self._sci_btn.setEnabled(False); self._sci_progress.setVisible(True)
        self._sci_status.setText("Generating..."); self._sci_frame.setVisible(True)
        self._worker = SCIPlotWorker(csv, outdir)
        self._worker.progress.connect(lambda m: self._sci_status.setText(m))
        self._worker.finished.connect(self._on_sci_done)
        self._worker.error.connect(self._on_sci_error); self._worker.start()

    def _on_sci_done(self, png, pdf):
        self._sci_btn.setEnabled(True); self._sci_progress.setVisible(False)
        self._sci_status.setText("SCI figure ready")
        self._png_path = png; self._pdf_path = pdf; self._preview.set_image(png)

    def _on_sci_error(self, msg):
        self._sci_btn.setEnabled(True); self._sci_progress.setVisible(False); self._sci_status.setText(f"Error: {msg}")

    def _save_png(self):
        if hasattr(self,'_png_path') and os.path.isfile(self._png_path):
            dest,_ = QFileDialog.getSaveFileName(self,"Save PNG","SCI_Figure.png","PNG(*.png)")
            if dest: import shutil; shutil.copy2(self._png_path, dest)

    def _save_pdf(self):
        if hasattr(self,'_pdf_path') and os.path.isfile(self._pdf_path):
            dest,_ = QFileDialog.getSaveFileName(self,"Save PDF","SCI_Figure.pdf","PDF(*.pdf)")
            if dest: import shutil; shutil.copy2(self._pdf_path, dest)

    # ── Charts ─────────────────────────────────
    def _plot_all(self):
        # Clear
        while self._chart_layout.count():
            item = self._chart_layout.takeAt(0)
            if item.widget(): item.widget().deleteLater()
            elif item.layout(): self._cl(item.layout())

        df = self._store.dataframe
        if df is None or len(df) == 0:
            self._chart_layout.addWidget(QLabel("No data to visualize.")); return
        field = self._field_combo.currentText()

        # ── 2-column grid for regular charts ──
        grid = QGridLayout(); grid.setSpacing(12)
        charts = []
        if field in df.columns:
            charts.append((field, field+" Distribution"))
        for c in ["ScientificName","Tissue","Source","LibrarySource","Location"]:
            if c in df.columns and c != field:
                charts.append((c, c+" Distribution"))
        for i, (col, title) in enumerate(charts):
            grid.addWidget(self._mk_bar(df, col, title), i//2, i%2)
        self._chart_layout.addLayout(grid)

        # ── Full-width charts ──
        self._chart_layout.addWidget(self._mk_heatmap(df))
        if "CollectionDate" in df.columns:
            self._chart_layout.addWidget(self._mk_timeline(df, "CollectionDate", "Collection Timeline"))
        self._chart_layout.addWidget(self._mk_stats(df))
        self._chart_layout.addStretch()

    def _cl(self, layout):
        while layout.count():
            item = layout.takeAt(0)
            if item.widget(): item.widget().deleteLater()
            elif item.layout(): self._cl(item.layout())

    def _mk_bar(self, df, col, title):
        series = df[col].apply(lambda x: str(x).strip())
        series = series[~series.apply(self._store.is_missing)]
        g = QGroupBox(title); g.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.MinimumExpanding)
        lay = QVBoxLayout(g)
        if len(series) == 0: lay.addWidget(QLabel("No data")); return g
        counts = Counter(series); most = counts.most_common(20)
        labels = [k if len(k) <= 30 else k[:29]+"…" for k,_ in most]
        values = [v for _, v in most]
        h = max(BAR_MIN_H, len(labels)*BAR_ITEM_H)
        fig = Figure(figsize=(FIG_W, h), dpi=FIG_DPI); ax = fig.add_subplot(111)
        ax.barh(range(len(labels)), values, color="#4C72B0", height=0.7)
        ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("Count", fontsize=10); ax.invert_yaxis()
        for i, v in enumerate(values):
            ax.text(v+max(values)*0.01, i, str(v), va="center", fontsize=8)
        fig.tight_layout(pad=0.3)
        canvas = FigureCanvas(fig); canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        fig_h = fig.get_size_inches()[1]
        canvas.setFixedHeight(int(fig_h*fig.dpi*1.1))
        canvas.wheelEvent = lambda e: e.ignore()
        lay.addWidget(canvas)
        # Save button
        btn = QPushButton("Save PNG"); btn.setMaximumWidth(100)
        btn.clicked.connect(lambda checked, f=fig, t=title: self._save_chart(f, t))
        hb = QHBoxLayout(); hb.addStretch(); hb.addWidget(btn); lay.addLayout(hb)
        return g

    def _mk_heatmap(self, df):
        g = QGroupBox("Missing Data Overview"); g.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.MinimumExpanding)
        lay = QVBoxLayout(g); n = len(df)
        if n == 0: lay.addWidget(QLabel("No data")); return g
        cols = list(df.columns); data = np.zeros((len(cols),1)); labels = []
        for i, col in enumerate(cols):
            s = df[col].apply(lambda x: str(x).strip())
            data[i,0] = s.apply(self._store.is_missing).sum()/n*100; labels.append(col)
        h = max(3.0, len(labels)*HEAT_ITEM_H)
        fig = Figure(figsize=(FULL_W, h), dpi=FIG_DPI); ax = fig.add_subplot(111)
        cmap = plt.cm.RdYlGn_r
        ax.imshow(data.T, aspect="auto", cmap=cmap, vmin=0, vmax=100)
        ax.set_yticks([0]); ax.set_yticklabels(["Missing %"], fontsize=10)
        ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
        for i in range(len(labels)):
            v = data[i,0]; ax.text(i, 0, f"{v:.0f}%", ha="center", va="center",
                                    fontsize=9, fontweight="bold", color="white" if v>50 else "black")
        fig.tight_layout(pad=0.3)
        canvas = FigureCanvas(fig); canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        canvas.setFixedHeight(int(fig.get_size_inches()[1]*fig.dpi*1.1))
        canvas.wheelEvent = lambda e: e.ignore()
        lay.addWidget(canvas)
        btn = QPushButton("Save PNG"); btn.setMaximumWidth(100)
        btn.clicked.connect(lambda checked, f=fig: self._save_chart(f, "Missing_Data"))
        hb = QHBoxLayout(); hb.addStretch(); hb.addWidget(btn); lay.addLayout(hb)
        return g

    def _mk_timeline(self, df, col, title):
        g = QGroupBox(title); lay = QVBoxLayout(g)
        series = df[col].apply(lambda x: str(x).strip())
        series = series[~series.apply(self._store.is_missing)]
        if len(series) == 0: lay.addWidget(QLabel("No date data")); return g
        year_counts = Counter()
        for v in series:
            parts = v.replace("/","-").split("-")
            if parts[0].isdigit() and len(parts[0])==4: year_counts[parts[0]] += 1
        if not year_counts: lay.addWidget(QLabel("No parseable dates")); return g
        years = sorted(year_counts.keys()); counts = [year_counts[y] for y in years]
        fig = Figure(figsize=(FULL_W, TIME_H), dpi=FIG_DPI); ax = fig.add_subplot(111)
        ax.bar(years, counts, color="#55A868", width=0.6)
        ax.set_xlabel("Year", fontsize=10); ax.set_ylabel("Records", fontsize=10)
        for i,(y,c) in enumerate(zip(years,counts)): ax.text(i, c+max(counts)*0.02, str(c), ha="center", fontsize=9)
        fig.tight_layout(pad=0.3)
        canvas = FigureCanvas(fig); canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        canvas.setFixedHeight(int(fig.get_size_inches()[1]*fig.dpi*1.1))
        canvas.wheelEvent = lambda e: e.ignore()
        lay.addWidget(canvas)
        btn = QPushButton("Save PNG"); btn.setMaximumWidth(100)
        btn.clicked.connect(lambda checked, f=fig, t=title: self._save_chart(f, t))
        hb = QHBoxLayout(); hb.addStretch(); hb.addWidget(btn); lay.addLayout(hb)
        return g

    def _mk_stats(self, df):
        g = QGroupBox("Data Summary"); lay = QVBoxLayout(g)
        total, ncols = len(df), len(df.columns); cells = total*ncols
        filled = sum((~df[c].apply(lambda x: str(x).strip()).apply(self._store.is_missing)).sum() for c in df.columns)
        missing = cells-filled; complete = sum(1 for _,row in df.iterrows() if all(not self._store.is_missing(v) for v in row))
        pct = filled/cells*100
        lbl = QLabel(f"Total: {total} records  |  Fields: {ncols}  |  Filled: {filled} ({pct:.1f}%)  |  Missing: {missing}\nComplete records: {complete} ({complete/total*100:.1f}%)")
        lbl.setStyleSheet("font-size:12px;padding:12px;line-height:1.8"); lay.addWidget(lbl)
        return g

    def _save_chart(self, fig, name):
        safe = "".join(c if c.isalnum() or c in "._- " else "_" for c in name).replace(" ","_")
        dest, _ = QFileDialog.getSaveFileName(self, "Save Chart", f"{safe}.png", "PNG (*.png)")
        if dest: fig.savefig(dest, dpi=150, bbox_inches="tight")
