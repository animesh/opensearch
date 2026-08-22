#!/usr/bin/env python3
"""Count spectra and MS2 spectra in mzML without loading the document into memory."""
from __future__ import annotations
import gzip
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

MS_LEVEL_ACCESSION = "MS:1000511"


def open_xml(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def count(path: Path) -> tuple[int, int]:
    total = 0
    ms2 = 0
    with open_xml(path) as fh:
        for _, elem in ET.iterparse(fh, events=("end",)):
            if elem.tag.rsplit("}", 1)[-1] != "spectrum":
                continue
            total += 1
            level = None
            for child in elem.iter():
                if child.tag.rsplit("}", 1)[-1] == "cvParam" and child.attrib.get("accession") == MS_LEVEL_ACCESSION:
                    level = child.attrib.get("value")
                    break
            if level == "2":
                ms2 += 1
            elem.clear()
    return total, ms2


def main() -> int:
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} INPUT.mzML OUTPUT.tsv", file=sys.stderr)
        return 2
    mzml = Path(sys.argv[1])
    out = Path(sys.argv[2])
    if not mzml.is_file() or mzml.stat().st_size == 0:
        print(f"ERROR: mzML does not exist or is empty: {mzml}", file=sys.stderr)
        return 2
    try:
        total, ms2 = count(mzml)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8") as fh:
            fh.write("sample\ttotal_spectra\tms2_spectra\n")
            fh.write(f"{mzml.stem.removesuffix('_calibrated')}\t{total}\t{ms2}\n")
        print(f"Counted {total} spectra; {ms2} MS2 spectra: {mzml}")
        return 0
    except Exception as exc:
        print(f"ERROR: failed to count mzML spectra: {exc}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
