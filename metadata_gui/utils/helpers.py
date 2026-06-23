"""Utility helpers for the metadata GUI application."""

import os
import sys
from pathlib import Path
from typing import Optional


def find_project_root() -> Optional[str]:
    """Find the MMPV-RNA project root directory."""
    current = Path(__file__).resolve().parent.parent
    # Check if we're in metadata_gui/
    if current.name == "metadata_gui":
        parent = current.parent
        if (parent / "public_metadata_pipeline").is_dir():
            return str(parent)
        if (parent / "pipeline_release").is_dir():
            return str(parent)
    return None


def find_data_dir() -> Optional[str]:
    """Auto-discover the metadata output directory."""
    root = find_project_root()
    if not root:
        return None

    candidates = [
        os.path.join(root, "public_metadata_pipeline",
                     "public_data_pipeline_output", "info"),
        os.path.join(root, "public_data_pipeline_output", "info"),
    ]
    for d in candidates:
        if os.path.isdir(d):
            return d
    return None


def format_file_size(size_bytes: int) -> str:
    """Format bytes to human-readable string."""
    for unit in ("B", "KB", "MB", "GB"):
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} TB"
