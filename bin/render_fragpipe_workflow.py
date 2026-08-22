#!/usr/bin/env python3
"""Render a FragPipe workflow for one input sample.

The supplied workflow remains the source of truth for search parameters.  Only
settings that depend on the execution environment or the input file type are
replaced.  All other lines, including comments and ordering, are preserved.
"""
from __future__ import annotations

import argparse
from pathlib import Path

OVERRIDE_KEYS = {
    "database.db-path",
    "fragpipe-config.tools-folder",
    "fragpipe-config.bin-diann",
    "fragpipe-config.bin-python",
    "crystalc.run-crystalc",
    "workflow.input.data-type.im-ms",
    "workflow.input.data-type.regular-ms",
    "msfragger.write_calibrated_mzml",
}
REQUIRED_KEYS = {"database.db-path", "msfragger.run-msfragger"}


def parse_bool(value: str) -> str:
    value = value.strip().lower()
    if value not in {"true", "false"}:
        raise ValueError(f"expected true/false, got {value!r}")
    return value


def read_values(lines: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in lines:
        text = line.strip()
        if not text or text.startswith("#") or "=" not in text:
            continue
        key, value = text.split("=", 1)
        values[key] = value
    return values


def render(args: argparse.Namespace) -> None:
    source = Path(args.source)
    destination = Path(args.destination)
    if not source.is_file():
        raise FileNotFoundError(source)

    lines = source.read_text(encoding="utf-8-sig").replace("\r\n", "\n").splitlines()
    values = read_values(lines)

    database = args.database or values.get("database.db-path", "")
    if not database:
        raise ValueError("FragPipe workflow has no database.db-path and --database was not supplied")

    if args.crystalc == "auto":
        crystalc = "false" if args.im_ms == "true" else values.get("crystalc.run-crystalc", "false")
    else:
        crystalc = parse_bool(args.crystalc)

    im_ms = parse_bool(args.im_ms)
    regular_ms = "false" if im_ms == "true" else "true"

    overrides = {
        "database.db-path": database,
        "workflow.input.data-type.im-ms": im_ms,
        "workflow.input.data-type.regular-ms": regular_ms,
        # Downstream Casanovo/AA_stat/QC need an mzML when FragPipe can create one.
        "msfragger.write_calibrated_mzml": "true" if args.write_calibrated else "false",
        "crystalc.run-crystalc": crystalc,
    }
    if args.tools_folder:
        overrides["fragpipe-config.tools-folder"] = args.tools_folder
    if args.diann:
        overrides["fragpipe-config.bin-diann"] = args.diann
    if args.python:
        overrides["fragpipe-config.bin-python"] = args.python

    kept = []
    for line in lines:
        text = line.strip()
        if text and not text.startswith("#") and "=" in text:
            key = text.split("=", 1)[0]
            if key in OVERRIDE_KEYS:
                continue
        kept.append(line)

    prefix = [f"{key}={value}" for key, value in overrides.items()]
    rendered = "\n".join(prefix + kept).rstrip() + "\n"
    final_values = read_values(rendered.splitlines())
    missing = sorted(k for k in REQUIRED_KEYS if not final_values.get(k))
    if missing:
        raise ValueError("rendered workflow is missing required settings: " + ", ".join(missing))

    destination.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--destination", required=True)
    ap.add_argument("--database", default="")
    ap.add_argument("--tools-folder", default="")
    ap.add_argument("--diann", default="")
    ap.add_argument("--python", default="")
    ap.add_argument("--im-ms", required=True, choices=["true", "false"])
    ap.add_argument("--crystalc", default="auto", choices=["auto", "true", "false"])
    ap.add_argument("--write-calibrated", action="store_true")
    args = ap.parse_args()
    render(args)
