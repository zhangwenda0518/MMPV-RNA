#!/usr/bin/env python3
"""Metadata Manager - Desktop GUI for public metadata pipeline data."""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    from PySide6.QtWidgets import QApplication
    from PySide6.QtGui import QFont

    app = QApplication(sys.argv)
    app.setApplicationName("Metadata Manager")
    app.setOrganizationName("MMPV-RNA")
    app.setFont(QFont("Segoe UI", 10))
    app.setStyle("Fusion")

    from views.main_window import MainWindow
    window = MainWindow()
    window.resize(1400, 900)
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
