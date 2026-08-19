#!/usr/bin/env python3
"""Thin command-line wrapper for the module-local OpenSearch report implementation."""
from pathlib import Path
import runpy

runpy.run_path(
    str(Path(__file__).resolve().parents[1] / "modules/local/opensearch_summary/bin/opensearch_summary.py"),
    run_name="__main__",
)
