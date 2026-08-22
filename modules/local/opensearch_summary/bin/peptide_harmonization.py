#!/usr/bin/env python3
"""Harmonize FragPipe and Casanovo peptide/peptidoform representations.

The comparison deliberately keeps two concepts separate:
1. sequence identity after removing modifications, and
2. peptidoform identity after position-aware modification normalization.

Unknown modifications are retained as unknown and are never guessed.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

# Monoisotopic masses (Da) for common Unimod entries observed in Casanovo output.
UNIMOD_MASS = {
    "1": 42.010565,       # Acetyl
    "4": 57.021464,       # Carbamidomethyl
    "5": 43.005814,       # Carbamyl
    "7": 0.984016,        # Deamidated
    "28": -17.026549,     # Gln -> pyro-Glu
    "35": 15.994915,      # Oxidation
    "40": 79.956815,      # Sulfo
    "141": 41.026549,     # Amidination
    "351": 3.995796,      # Trp -> Kynurenin
    "374": -1.002625,     # Dehydro
    "385": -17.026549,    # Ammonia-loss (N-term)
}

NAME_MASS = {
    "acetyl": 42.010565,
    "carbamidomethyl": 57.021464,
    "carbamyl": 43.005814,
    "deamidated": 0.984016,
    "oxidation": 15.994915,
    "sulfo": 79.956815,
    "amidination": 41.026549,
    "ammonia-loss": -17.026549,
    "ammonia loss": -17.026549,
}


def _mass_from_name(name: str):
    key = re.sub(r"\s+", " ", name.strip().lower())
    for k, mass in NAME_MASS.items():
        if k in key:
            return mass
    return None


def normalize_sequence(seq: str, il_equivalent: bool = False) -> str:
    seq = re.sub(r"[^A-Za-z]", "", seq or "").upper()
    return seq.replace("I", "L") if il_equivalent else seq


def parse_fragpipe_mods(peptide: str, assigned: str):
    """Return sorted (position, mass) modifications from Assigned Modifications."""
    mods = []
    for item in (assigned or "").split(","):
        item = item.strip()
        if not item:
            continue
        m = re.fullmatch(r"N-term\((-?\d+(?:\.\d+)?)\)", item)
        if m:
            mods.append((0, float(m.group(1))))
            continue
        m = re.fullmatch(r"(\d+)[A-Z]\((-?\d+(?:\.\d+)?)\)", item)
        if m:
            mods.append((int(m.group(1)), float(m.group(2))))
            continue
    return tuple(sorted(mods))


def parse_casanovo_mods(modifications: str, proforma: str, sequence: str):
    """Return sorted (position, mass) modifications from Casanovo mzTab.

    Position 0 is the peptide N-terminus. Residue positions are 1-based.
    """
    mods = []
    text = (modifications or "").strip()
    if text and text.lower() != "null":
        for item in text.split(";"):
            m = re.match(r"\s*(\d+)-(.+?)(?::UNIMOD:(\d+))?\s*$", item)
            if not m:
                continue
            pos = int(m.group(1))
            name = m.group(2)
            uid = m.group(3)
            mass = UNIMOD_MASS.get(uid) if uid else _mass_from_name(name)
            if mass is not None:
                mods.append((pos, mass))
            else:
                mods.append((pos, f"UNKNOWN:{name}"))
    return tuple(sorted(mods, key=lambda x: (x[0], str(x[1]))))


def mods_equal(a, b, tolerance=0.05):
    if len(a) != len(b):
        return False
    for (pa, ma), (pb, mb) in zip(a, b):
        if pa != pb:
            return False
        if isinstance(ma, str) or isinstance(mb, str):
            if ma != mb:
                return False
        elif abs(ma - mb) > tolerance:
            return False
    return True


def parse_fp_key(spectrum: str):
    m = re.search(r"\.(\d+)\.\d+\.(\d+)$", spectrum or "")
    return (int(m.group(1)), int(m.group(2))) if m else None


def parse_cas_key(ref: str, charge: str, scan_number: str = ""):
    if scan_number and re.fullmatch(r"\d+", scan_number.strip()):
        return int(scan_number), int(float(charge))
    m = re.search(r"\bscan=(\d+)", ref or "")
    if m:
        return int(m.group(1)), int(float(charge))
    m = re.search(r"(?:index=)(\d+)", ref or "")
    if m:
        # index is zero-based and is not guaranteed to equal an instrument scan;
        # retain it only as a fallback because the same convention must be used
        # by both tools to make a valid spectrum key.
        return int(m.group(1)) + 1, int(float(charge))
    return None


def _better_fp(a, b):
    """Choose one FragPipe PSM per spectrum deterministically."""
    def score(r):
        q = float(r.get("Qvalue") or 1.0)
        prob = float(r.get("Probability") or 0.0)
        hyper = float(r.get("Hyperscore") or 0.0)
        return (-q, prob, hyper)
    return a if score(a) >= score(b) else b


def load_fragpipe(psm_path: Path):
    by_key = {}
    with psm_path.open(encoding="utf-8", errors="replace", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if str(r.get("Is Decoy", "")).lower() == "true" or str(r.get("Is Contaminant", "")).lower() == "true":
                continue
            key = parse_fp_key(r.get("Spectrum", ""))
            if not key:
                continue
            if key not in by_key or _better_fp(r, by_key[key]) is r:
                by_key[key] = r
    return by_key


def load_casanovo(mztab_path: Path, score_threshold=0.5):
    by_key = {}
    with mztab_path.open(encoding="utf-8", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("PSH\t"):
                header = line.rstrip("\n").split("\t")
                continue
            if not header or not line.startswith("PSM\t"):
                continue
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            try:
                score = float(r.get("search_engine_score[1]", "nan"))
            except ValueError:
                continue
            if score < score_threshold:
                continue
            key = parse_cas_key(r.get("spectra_ref", ""), r.get("charge", "0"), r.get("opt_global_cv_MS:1003057_scan_number", ""))
            # v5.1.x uses opt_ms_run[1]_proforma; v5.2 uses opt_global_cv_MS:1003169_proforma_peptidoform_sequence.
            proforma = r.get("opt_global_cv_MS:1003169_proforma_peptidoform_sequence") or r.get("opt_ms_run[1]_proforma") or r.get("sequence", "")
            r["_proforma"] = proforma
            r["_score"] = score
            if key and key not in by_key:
                by_key[key] = r
    return by_key


def compare(frag_by_key, cas_by_key, tolerance=0.05):
    rows = []
    counts = defaultdict(int)
    for key in sorted(set(frag_by_key) & set(cas_by_key)):
        f = frag_by_key[key]
        c = cas_by_key[key]
        fs = normalize_sequence(f.get("Peptide", ""))
        cs = normalize_sequence(c.get("sequence", ""))
        fs_il = normalize_sequence(fs, True)
        cs_il = normalize_sequence(cs, True)
        fmods = parse_fragpipe_mods(fs, f.get("Assigned Modifications", ""))
        cmods = parse_casanovo_mods(c.get("modifications", ""), c.get("_proforma", ""), cs)
        same_seq = fs == cs
        same_il = fs_il == cs_il
        same_mod = same_seq and mods_equal(fmods, cmods, tolerance)
        if same_mod:
            category = "same_sequence_exact_peptidoform"
        elif same_seq:
            category = "same_sequence_modification_disagreement"
        else:
            category = "different_sequence"
        counts[category] += 1
        rows.append({
            "scan": key[0], "charge": key[1], "category": category,
            "fragpipe_sequence": fs, "casanovo_sequence": cs,
            "same_sequence": same_seq, "same_sequence_IL_equivalent": same_il,
            "fragpipe_mods": ";".join(f"{p}:{m:.4f}" if not isinstance(m, str) else f"{p}:{m}" for p, m in fmods),
            "casanovo_mods": ";".join(f"{p}:{m:.4f}" if not isinstance(m, str) else f"{p}:{m}" for p, m in cmods),
            "fragpipe_modified_peptide": f.get("Modified Peptide", ""),
            "casanovo_proforma": c.get("_proforma", ""),
            "casanovo_score": c.get("_score", ""),
        })
    return rows, dict(counts)
