#!/usr/bin/env python3
"""Metadata Manager - Desktop GUI for public metadata pipeline data."""

import sys, os, traceback, datetime

LOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "error.log")
def _log(msg):
    try:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(f"[{datetime.datetime.now()}] {msg}\n")
    except: pass

try:
    _log("main.py starting")
    import pandas as pd
    import matplotlib
    matplotlib.use("QtAgg")
    _log("pandas+matplotlib OK")

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    from PySide6.QtWidgets import QApplication
    from PySide6.QtGui import QFont
    _log("PySide6 OK")

    app = QApplication(sys.argv)
    app.setApplicationName("Metadata Manager")
    app.setOrganizationName("MMPV-RNA")
    app.setFont(QFont("Segoe UI", 10))
    app.setStyle("Fusion")
    _log("QApplication OK")

    from views.main_window import MainWindow
    _log("MainWindow module loaded")
    window = MainWindow()
    _log(f"MainWindow created, search_rows={window.search_store.row_count}")
    window.resize(1400, 900)
    window.show()
    _log("Window shown, entering event loop")
    sys.exit(app.exec())
except Exception:
    _log(f"FATAL:\n{traceback.format_exc()}")

