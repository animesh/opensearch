#!/usr/bin/env python3
"""
OpenSearch integrated MultiQC summary v3.

Key changes from v2:
  * Uses FragPipe psm.tsv + peptide.tsv + protein.tsv directly.
  * Uses Casanovo mzTab sequence AND ProForma columns.
  * Separates spectrum-level and peptide-level comparisons.
  * Reports:
      - Casanovo attempted / MS2
      - Casanovo >=0.5 / MS2
      - Casanovo >=0.5 / attempted
      - Casanovo confirmation by FragPipe
      - FragPipe confirmation by Casanovo
  * Peptide comparison has three levels:
      1. stripped amino-acid sequence
      2. stripped sequence after optional I/L equivalence
      3. modification-aware peptidoform comparison
  * Modification-aware comparison harmonizes common FragPipe integer-mass
    notation with Casanovo ProForma modification names/masses.
  * Keeps an explicit "same sequence, different modification" category.
  * Uses only high-confidence Casanovo PSMs (default >=0.50) for the
    primary de-novo peptide comparison.
  * Adds protein recovery of the smaller run alongside Jaccard.
  * Produces per-metric QC flags rather than one dominant flag.
  * Keeps machine-readable summary.tsv and provenance.tsv.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

MOD_TOL_DA = 0.05
DEFAULT_CASANOVO_THRESHOLD = 0.50

# Common names emitted by Casanovo / common proteomics notation.
# Numeric matching is deliberately used only when a mass is available.
COMMON_MOD_MASSES = {
    "oxidation": 15.994915,
    "oxidation(m)": 15.994915,
    "oxidation of m": 15.994915,
    "carbamidomethyl": 57.021464,
    "carbamidomethylation": 57.021464,
    "acetyl": 42.010565,
    "acetylation": 42.010565,
    "deamidated": 0.984016,
    "deamidation": 0.984016,
    "gln->pyro-glu": -17.026549,
    "gln->pyro-glu (n-term q)": -17.026549,
    "pyro-glu": -17.026549,
    "met->hsl": -48.003371,
    "trp->kynurenin": 3.9949,
    "sulfo": 79.956815,
    "phospho": 79.966331,
    "methyl": 14.015650,
    "dimethyl": 28.031300,
    "formyl": 27.994915,
    "amidated": -0.984016,
    "amidation": -0.984016,
}

# Useful UNIMOD IDs for common modifications.
UNIMOD_MASSES = {
    1: 42.010565,   # Acetyl
    4: 57.021464,   # Carbamidomethyl
    7: 0.984016,    # Deamidated
    21: 79.966331,  # Phospho
    35: 15.994915,  # Oxidation
    40: 79.956815,  # Sulfo
    26: 14.015650,  # Methyl
    36: 28.031300,  # Dimethyl
    2: 14.015650,   # Amidated is not actually Unimod 2; retained only as fallback below
}


def read_tsv(path: Path | None):
    if not path or not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def first(path: Path | None, pattern: str):
    return next(path.glob(pattern), None) if path and path.exists() else None


def truth(value) -> bool:
    return str(value or "").strip().lower() in {"true", "1", "yes"}


def number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def int_number(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def pct(n, d):
    return 100.0 * n / d if d else 0.0


def median(values):
    values = sorted(v for v in values if v is not None)
    if not values:
        return 0.0
    n = len(values)
    return values[n // 2] if n % 2 else (values[n // 2 - 1] + values[n // 2]) / 2


def sample_name(path: Path, suffix: str) -> str:
    return path.name[:-len(suffix)] if suffix and path.name.endswith(suffix) else path.name


def strip_html(text: str) -> str:
    return html.unescape(re.sub(r"<[^>]*>", "", text or "")).strip()


def match_pct(label: str):
    m = re.search(r"\(([0-9.]+)%\s+match\)", label or "")
    return float(m.group(1)) if m else None


def slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", text).strip("_").lower() or "section"


def write_json(name, obj):
    Path(name).write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def mqc_table(name, section, description, data):
    write_json(name, {
        "id": slug(section),
        "section_name": section,
        "description": description,
        "plot_type": "table",
        "pconfig": {"id": slug(section) + "_table"},
        "data": data,
    })


def mqc_bar(name, section, title, ylab, data):
    write_json(name, {
        "id": slug(section),
        "section_name": section,
        "plot_type": "bargraph",
        "pconfig": {
            "id": slug(section) + "_plot",
            "title": title,
            "ylab": ylab,
            "cpswitch": True,
        },
        "data": data,
    })


def html_section(section_id, section_name, description, body):
    meta = "\n".join([
        "<!--",
        f"id: {section_id}",
        f"section_name: {section_name}",
        f"description: {description}",
        "-->",
    ])
    Path(f"{section_id}_mqc.html").write_text(meta + body, encoding="utf-8")


def relative_source(tool: str, directory: Path, file: Path | None = None) -> str:
    base = {"FragPipe": "fragpipe", "Casanovo": "casanovo", "AA_stat": "aa_stat"}[tool]
    target = directory.name
    if file:
        try:
            rel = file.relative_to(directory)
            return f"../{base}/{target}/{rel.as_posix()}"
        except ValueError:
            pass
    return f"../{base}/{target}/"


def source_cell(tool: str, path: str):
    label = f"{tool}: {path}"
    return f'<a href="{html.escape(path, quote=True)}">{html.escape(label)}</a>'


# ---------------------------------------------------------------------------
# Modification / peptide normalization
# ---------------------------------------------------------------------------

def normalize_name(name: str) -> str:
    x = str(name or "").strip().lower()
    x = x.replace("−", "-").replace("→", "->")
    x = re.sub(r"\s+", " ", x)
    return x


def parse_mod_token(token: str):
    """
    Return ('mass', mass_da) whenever the modification can be mapped to a
    numeric mass, otherwise ('name', normalized_name).

    This intentionally does not invent masses for unknown modifications.
    """
    token = str(token or "").strip()
    low = normalize_name(token)

    # ProForma UNIMOD:nn
    m = re.fullmatch(r"unimod:(\d+)", low)
    if m:
        uid = int(m.group(1))
        if uid in UNIMOD_MASSES:
            return ("mass", UNIMOD_MASSES[uid])

    # Pure numeric token, including +16, -17.0265, 57, etc.
    if re.fullmatch(r"[+-]?\d+(?:\.\d+)?", token):
        return ("mass", float(token))

    # Some FragPipe/ProForma strings have a leading +/- numeric value.
    m = re.fullmatch(r"([+-]\d+(?:\.\d+)?)", token)
    if m:
        return ("mass", float(m.group(1)))

    # Named common modification.
    if low in COMMON_MOD_MASSES:
        return ("mass", COMMON_MOD_MASSES[low])

    # Strip residue qualifiers in common names.
    simplified = re.sub(r"\([a-z]\)$", "", low)
    if simplified in COMMON_MOD_MASSES:
        return ("mass", COMMON_MOD_MASSES[simplified])

    return ("name", re.sub(r"\s+", "", low))


def mod_equal(a, b, tol=MOD_TOL_DA):
    ak, av = a
    bk, bv = b
    if ak in {"mass", "m"} and bk in {"mass", "m"}:
        return abs(float(av) - float(bv)) <= tol
    if ak in {"name", "n"} and bk in {"name", "n"}:
        return av == bv
    return False


def canonical_mod(mod):
    kind, value = mod
    if kind == "mass":
        # Keep enough precision for diagnostics. Equality is handled by
        # mod_equal(), not by exact tuple equality, because FragPipe often
        # reports integer modification masses whereas Casanovo ProForma may
        # report the corresponding precise mass.
        return ("m", round(float(value), 4))
    return ("n", str(value))


def parse_modified_peptide(text: str):
    """
    Parse FragPipe-style modified peptide:
        M[16]PEPTC[57]IDE
        [42]PEPTIDE
        PEPTIDE[80]

    Position 0 = peptide N-terminus.
    Residues are numbered 1..N.
    Position N+1 is reserved for an explicit C-terminal modification.

    Returns:
        stripped_sequence, tuple((position, canonical_mod), ...)
    """
    text = str(text or "").strip()
    if not text:
        return "", ()

    sequence = []
    mods = []
    i = 0
    pos = 0

    while i < len(text):
        if text[i] == "[":
            j = text.find("]", i + 1)
            if j < 0:
                break
            tok = text[i + 1:j]
            position = 0 if pos == 0 else pos + 1
            mods.append((position, canonical_mod(parse_mod_token(tok))))
            i = j + 1
            continue

        ch = text[i]
        if ch.isalpha() and ch.isupper():
            sequence.append(ch)
            pos += 1
            i += 1

            # Residue-attached modification.
            if i < len(text) and text[i] == "[":
                j = text.find("]", i + 1)
                if j >= 0:
                    tok = text[i + 1:j]
                    mods.append((pos, canonical_mod(parse_mod_token(tok))))
                    i = j + 1
            continue

        # ProForma separators or harmless punctuation.
        i += 1

    return "".join(sequence), tuple(sorted(mods))


def parse_proforma(text: str):
    """
    Parse the Casanovo ProForma column.

    Supports common forms:
        PEPTM[Oxidation]IDE
        [Acetyl]-PEPTIDE
        PEPTIDE-[Amidated]
        [UNIMOD:35]PEPTIDE

    This is intentionally conservative. Unknown modifications remain named
    rather than being guessed.
    """
    text = str(text or "").strip()
    if not text:
        return "", ()

    # Remove ProForma charge/adduct annotations that are not part of the
    # peptidoform sequence. We do not remove residue modifications.
    text = re.sub(r"/\d+$", "", text)

    sequence = []
    mods = []
    i = 0
    pos = 0

    while i < len(text):
        if text[i] == "[":
            j = text.find("]", i + 1)
            if j < 0:
                break
            tok = text[i + 1:j]
            position = 0 if pos == 0 else pos + 1
            mods.append((position, canonical_mod(parse_mod_token(tok))))
            i = j + 1
            continue

        ch = text[i]
        if ch.isalpha() and ch.isupper():
            sequence.append(ch)
            pos += 1
            i += 1

            if i < len(text) and text[i] == "[":
                j = text.find("]", i + 1)
                if j >= 0:
                    tok = text[i + 1:j]
                    mods.append((pos, canonical_mod(parse_mod_token(tok))))
                    i = j + 1
            continue

        i += 1

    return "".join(sequence), tuple(sorted(mods))


def strip_modifications(text: str):
    """
    Generic sequence stripper for ProForma / FragPipe representations.
    """
    seq, _ = parse_proforma(text)
    return seq


def il_key(sequence: str) -> str:
    # I/L are isobaric at MS1 and generally indistinguishable in DDA MS/MS.
    return sequence.upper().replace("I", "J").replace("L", "J")


def peptidoform_key(sequence: str, mods):
    return (sequence.upper(), tuple(mods))


def peptidoform_il_key(sequence: str, mods):
    return (il_key(sequence), tuple(mods))


def mods_match(mods_a, mods_b, tol=MOD_TOL_DA):
    """Position-aware, order-independent modification comparison."""
    if len(mods_a) != len(mods_b):
        return False

    unmatched = list(mods_b)
    for pos_a, mod_a in mods_a:
        found = None
        for i, (pos_b, mod_b) in enumerate(unmatched):
            if pos_a == pos_b and mod_equal(mod_a, mod_b, tol):
                found = i
                break
        if found is None:
            return False
        unmatched.pop(found)
    return not unmatched


def peptidoform_equal(a, b, il=False, tol=MOD_TOL_DA):
    seq_a, mods_a = a
    seq_b, mods_b = b
    if il:
        seq_a = il_key(seq_a)
        seq_b = il_key(seq_b)
    return seq_a == seq_b and mods_match(mods_a, mods_b, tol)


def count_peptidoform_matches(casanovo_forms, fragpipe_forms, il=False):
    """
    Count one-to-one peptidoform matches using position-aware modification
    matching and a mass tolerance. Exact duplicate peptidoforms are sets, so
    each FragPipe form can be consumed at most once.
    """
    remaining = list(fragpipe_forms)
    shared = 0

    # Match sequences first to avoid unnecessary pairwise comparisons.
    fp_by_seq = defaultdict(list)
    for form in remaining:
        seq = il_key(form[0]) if il else form[0]
        fp_by_seq[seq].append(form)

    for cform in casanovo_forms:
        seq = il_key(cform[0]) if il else cform[0]
        candidates = fp_by_seq.get(seq, [])
        hit = None
        for i, fform in enumerate(candidates):
            if peptidoform_equal(cform, fform, il=il):
                hit = i
                break
        if hit is not None:
            shared += 1
            candidates.pop(hit)

    return shared


# ---------------------------------------------------------------------------
# Spectrum keys
# ---------------------------------------------------------------------------

def fragpipe_spectrum_key(row):
    spectrum = (row.get("Spectrum") or "").strip()
    charge = str(row.get("Charge") or "").strip()
    if spectrum:
        return (spectrum, charge)

    sf = (row.get("Spectrum File") or "").strip()
    scan = str(row.get("Scan") or "").strip()
    if scan:
        return (sf, scan, charge)
    return None


def casanovo_spectrum_key(row):
    ref = (row.get("spectra_ref") or "").strip()
    charge = str(row.get("charge") or "").strip()

    # mzTab: ms_run[1]:index=123
    m = re.search(r"(?:index|scan)=(\d+)", ref, flags=re.I)
    if m:
        return ("scan", m.group(1), charge)

    return (ref, charge) if ref else None


def normalize_scan_number_from_fragpipe(row):
    spectrum = (row.get("Spectrum") or "").strip()
    # Common FragPipe format: file.scan.scan.charge
    nums = re.findall(r"\.(\d+)\.", spectrum)
    if nums:
        return nums[-1]
    m = re.search(r"\.(\d+)(?:\.\d+)?\.\d+$", spectrum)
    return m.group(1) if m else ""


# ---------------------------------------------------------------------------
# FragPipe parsing
# ---------------------------------------------------------------------------

def parse_fragpipe(directory: Path):
    psm_file = first(directory, "**/psm.tsv")
    pep_file = first(directory, "**/peptide.tsv")
    prot_file = first(directory, "**/protein.tsv")
    mod_file = first(directory, "**/global.modsummary.tsv")
    spectrum_file = first(directory, "**/spectrum_count.tsv")

    psm = read_tsv(psm_file)
    peptides = read_tsv(pep_file)
    proteins = read_tsv(prot_file)

    peptide_set = set()
    for r in peptides:
        pep = (r.get("Peptide") or "").strip().upper()
        if pep and not truth(r.get("Is Decoy")) and not truth(r.get("Is Contaminant")):
            peptide_set.add(pep)

    protein_set = {
        (r.get("Protein ID") or r.get("Protein") or "").strip()
        for r in proteins
        if (r.get("Protein ID") or r.get("Protein"))
        and not truth(r.get("Is Decoy"))
        and not truth(r.get("Is Contaminant"))
    }

    charge = defaultdict(int)
    missed = defaultdict(int)
    modified = contaminants = missed_gt0 = 0
    psm_spectrum_keys = set()
    psm_rows_target = []

    # Modified peptidoforms from target/non-contaminant PSMs.
    peptidoforms = set()
    peptidoforms_il = set()
    psm_by_sequence = defaultdict(int)

    for r in psm:
        charge[str(r.get("Charge") or "?").strip()] += 1

        mc = int_number(r.get("Number of Missed Cleavages"), 0)
        missed[mc] += 1
        missed_gt0 += mc > 0

        is_cont = truth(r.get("Is Contaminant"))
        is_decoy = truth(r.get("Is Decoy"))
        contaminants += is_cont

        is_modified = bool(
            (r.get("Assigned Modifications") or "").strip()
            or (r.get("Observed Modifications") or "").strip()
            or (r.get("Modified Peptide") or "").strip()
        )
        modified += is_modified

        if not is_decoy and not is_cont:
            key = fragpipe_spectrum_key(r)
            if key:
                psm_spectrum_keys.add(key)
            psm_rows_target.append(r)

            pep = (r.get("Peptide") or "").strip().upper()
            psm_by_sequence[pep] += 1

            modified_pep = (r.get("Modified Peptide") or "").strip()
            if modified_pep:
                seq, mods = parse_modified_peptide(modified_pep)
                if seq:
                    peptidoforms.add(peptidoform_key(seq, mods))
                    peptidoforms_il.add(peptidoform_il_key(seq, mods))
            elif pep:
                peptidoforms.add(peptidoform_key(pep, ()))
                peptidoforms_il.add(peptidoform_il_key(pep, ()))

    modifications = []
    if mod_file:
        for r in read_tsv(mod_file):
            name = (r.get("Modification") or "").strip()
            if not name or name.lower() == "none":
                continue
            psm_key = next((k for k in r if k.endswith("_PSMs")), None)
            pct_key = next((k for k in r if k.endswith("_percent_PSMs")), None)
            modifications.append({
                "name": name,
                "mass": number(r.get("Mass Shift")),
                "psm": int_number(r.get(psm_key), 0),
                "pct": number(r.get(pct_key), 0),
                "isotope": "isotop" in name.lower(),
                "source": relative_source("FragPipe", directory, mod_file),
            })

    total_spectra = ms2_spectra = 0
    if spectrum_file:
        rows = read_tsv(spectrum_file)
        if rows:
            total_spectra = int_number(rows[0].get("total_spectra"), 0)
            ms2_spectra = int_number(rows[0].get("ms2_spectra"), 0)

    return {
        "program": "FragPipe",
        "directory": directory,
        "psm_file": psm_file,
        "peptide_file": pep_file,
        "protein_file": prot_file,
        "mod_file": mod_file,
        "spectrum_file": spectrum_file,
        "psm": len(psm),
        "psm_target_rows": len(psm_rows_target),
        "psm_spectrum_keys": psm_spectrum_keys,
        "peptide": len(peptide_set),
        "protein": len(protein_set),
        "modified_psm_pct": pct(modified, len(psm)),
        "contaminant_psm_pct": pct(contaminants, len(psm)),
        "missed_cleavage_pct": pct(missed_gt0, len(psm)),
        "charge": dict(charge),
        "missed": dict(missed),
        "peptides": peptide_set,
        "proteins": protein_set,
        "peptidoforms": peptidoforms,
        "peptidoforms_il": peptidoforms_il,
        "psm_by_sequence": dict(psm_by_sequence),
        "modifications": modifications,
        "total_spectra": total_spectra,
        "ms2_spectra": ms2_spectra,
    }


# ---------------------------------------------------------------------------
# Casanovo parsing
# ---------------------------------------------------------------------------

def parse_casanovo(directory: Path):
    log = first(directory, "casanovo_*.log")
    mztab = first(directory, "casanovo_*.mztab")
    text = log.read_text(encoding="utf-8", errors="replace") if log else ""

    thresholds = {}
    for threshold in ("0.00", "0.50", "0.90", "0.95", "0.99"):
        pattern = rf"([0-9]+)\s+spectra\s+\(([0-9.]+)%\)\s+scored\s+≥\s+{re.escape(threshold)}"
        m = re.search(pattern, text)
        thresholds[threshold] = (
            (int(m.group(1)), float(m.group(2))) if m else (0, 0.0)
        )

    m = re.search(r"Sequenced\s+([0-9]+)\s+spectra", text)
    sequenced = int(m.group(1)) if m else 0

    rows = []
    if mztab:
        header = None
        with mztab.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n\r")
                if line.startswith("PSH\t"):
                    header = line.split("\t")[1:]
                    continue
                if line.startswith("PSM\t") and header:
                    fields = line.split("\t")[1:]
                    row = dict(zip(header, fields))
                    rows.append(row)

    # mzTab row count is the safest fallback when the log wording changes.
    sequenced = sequenced or len(rows)

    # If the log does not contain threshold summaries, derive them directly
    # from the mzTab PSM scores. This keeps the parser robust to Casanovo log
    # wording changes.
    for threshold in ("0.00", "0.50", "0.90", "0.95", "0.99"):
        if thresholds[threshold][0] == 0 and rows:
            t = float(threshold)
            n = sum(
                number(r.get("search_engine_score[1]"), float("-inf")) >= t
                for r in rows
            )
            thresholds[threshold] = (n, pct(n, sequenced))

    all_sequences = set()
    high_sequences = set()
    all_peptidoforms = set()
    high_peptidoforms = set()
    all_peptidoforms_il = set()
    high_peptidoforms_il = set()

    all_spectrum_keys = set()
    high_spectrum_keys = set()

    score_values = []

    for r in rows:
        score = number(r.get("search_engine_score[1]"), float("-inf"))
        if math.isfinite(score):
            score_values.append(score)

        seq = (r.get("sequence") or "").strip().upper()
        proforma = (r.get(
            "opt_global_cv_MS:1003169_proforma_peptidoform_sequence"
        ) or "").strip()

        # Casanovo 5.x provides an unmodified sequence and a ProForma sequence.
        # Prefer the explicit sequence column for sequence-only comparison.
        if not seq and proforma:
            seq, _ = parse_proforma(proforma)

        if seq:
            all_sequences.add(seq)

        if proforma:
            pseq, mods = parse_proforma(proforma)
        else:
            pseq, mods = seq, ()

        if pseq:
            all_peptidoforms.add(peptidoform_key(pseq, mods))
            all_peptidoforms_il.add(peptidoform_il_key(pseq, mods))

        skey = casanovo_spectrum_key(r)
        if skey:
            all_spectrum_keys.add(skey)

        if score >= DEFAULT_CASANOVO_THRESHOLD:
            if seq:
                high_sequences.add(seq)
            if pseq:
                high_peptidoforms.add(peptidoform_key(pseq, mods))
                high_peptidoforms_il.add(peptidoform_il_key(pseq, mods))
            if skey:
                high_spectrum_keys.add(skey)

    return {
        "program": "Casanovo",
        "directory": directory,
        "log": log,
        "mztab": mztab,
        "sequenced": sequenced,
        "rows": rows,
        "thresholds": thresholds,
        "sequences": all_sequences,
        "high_sequences": high_sequences,
        "peptidoforms": all_peptidoforms,
        "high_peptidoforms": high_peptidoforms,
        "peptidoforms_il": all_peptidoforms_il,
        "high_peptidoforms_il": high_peptidoforms_il,
        "spectrum_keys": all_spectrum_keys,
        "high_spectrum_keys": high_spectrum_keys,
        "score_values": score_values,
    }


# ---------------------------------------------------------------------------
# AA_stat
# ---------------------------------------------------------------------------

def parse_aastat(directory: Path):
    stats = first(directory, "aa_statistics_table.csv")
    inter = first(directory, "interpretations.json")
    loc = first(directory, "localization_statistics.csv")
    interpretations = (
        json.loads(inter.read_text(encoding="utf-8", errors="replace"))
        if inter else {}
    )

    shifts = []
    if stats:
        with stats.open(encoding="utf-8", errors="replace", newline="") as fh:
            for r in csv.DictReader(fh):
                mass = number(r.get("mass shift"), None)
                if mass is None or truth(r.get("is reference")):
                    continue
                npep = int_number(r.get("# peptides in bin"), 0)
                if npep <= 0:
                    continue

                isotope = truth(r.get("is isotope"))
                key = None
                if interpretations:
                    numeric_keys = []
                    for k in interpretations:
                        try:
                            numeric_keys.append((abs(float(k) - mass), k))
                        except ValueError:
                            pass
                    if numeric_keys:
                        key = min(numeric_keys)[1]

                annotations = []
                if key:
                    for item in interpretations.get(key, []):
                        if item.get("type") == "isotope":
                            continue
                        label = item.get("label", "")
                        urls = re.findall(r'href=["\']([^"\']+)', label)
                        annotations.append({
                            "label": strip_html(label),
                            "match": match_pct(label) or 0.0,
                            "url": urls[0] if urls else "",
                        })

                annotations.sort(key=lambda x: x["match"], reverse=True)

                shifts.append({
                    "mass": mass,
                    "npep": npep,
                    "isotope": isotope,
                    "annotation": annotations[0] if annotations else {
                        "label": "Isotope" if isotope else "Unassigned",
                        "match": 0.0,
                        "url": "",
                    },
                    "source": relative_source("AA_stat", directory, stats),
                })

    localized = 0
    if loc:
        with loc.open(encoding="utf-8", errors="replace", newline="") as fh:
            for r in csv.DictReader(fh):
                if truth(r.get("is isotope")):
                    continue
                if int_number(r.get("# peptides in bin"), 0) <= 0:
                    continue
                loc_text = r.get("localization", "")
                pairs = re.findall(
                    r"['\"]([^'\"]+)['\"]\s*:\s*(\d+)", loc_text
                )
                localized += any(
                    k.lower() != "non-localized" and int(v) > 0
                    for k, v in pairs
                )

    return {
        "program": "AA_stat",
        "directory": directory,
        "stats": stats,
        "interpretations": inter,
        "localizations": loc,
        "shifts": shifts,
        "localized": localized,
    }


# ---------------------------------------------------------------------------
# Comparison metrics
# ---------------------------------------------------------------------------

def set_metrics(a, b):
    shared = len(a & b)
    return {
        "a": len(a),
        "b": len(b),
        "shared": shared,
        "a_only": len(a - b),
        "b_only": len(b - a),
        "union": len(a | b),
        "jaccard_pct": pct(shared, len(a | b)),
        "a_recovered_by_b_pct": pct(shared, len(a)),
        "b_recovered_by_a_pct": pct(shared, len(b)),
    }


def modification_difference_categories(casanovo_forms, fragpipe_forms):
    """
    Separate amino-acid agreement from modification agreement.

    "Same sequence, different modification" means that the amino-acid sequence
    occurs in both tools, but no position-aware modification match exists after
    mass/name harmonization.
    """
    cas_by_seq = defaultdict(set)
    fp_by_seq = defaultdict(set)

    for seq, mods in casanovo_forms:
        cas_by_seq[seq].add(mods)
    for seq, mods in fragpipe_forms:
        fp_by_seq[seq].add(mods)

    same_seq = set(cas_by_seq) & set(fp_by_seq)
    same_peptidoform = 0
    same_seq_different_mod = 0

    for seq in same_seq:
        cforms = [(seq, mods) for mods in cas_by_seq[seq]]
        fforms = [(seq, mods) for mods in fp_by_seq[seq]]
        if count_peptidoform_matches(cforms, fforms) > 0:
            same_peptidoform += 1
        else:
            same_seq_different_mod += 1

    return {
        "same_sequence": len(same_seq),
        "same_peptidoform": same_peptidoform,
        "same_sequence_different_modification": same_seq_different_mod,
        "casanovo_high_conf_unique_sequence": len(cas_by_seq),
        "fragpipe_unique_sequence": len(fp_by_seq),
    }


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------

def classify_metrics(summary):
    psm_values = [v["fragpipe_psm_spectra"] for v in summary.values() if v["fragpipe_psm_spectra"]]
    id_values = [v["fragpipe_id_pct"] for v in summary.values() if v["fragpipe_id_pct"]]
    q50_values = [v["cas_q50_pct"] for v in summary.values() if v["cas_q50_pct"]]
    mc_values = [v["missed_cleavage_pct"] for v in summary.values() if v["missed_cleavage_pct"]]

    psm_med = median(psm_values)
    id_med = median(id_values)
    q50_med = median(q50_values)

    result = {}

    for sample, v in summary.items():
        flags = []

        if psm_med and v["fragpipe_psm_spectra"] < 0.75 * psm_med:
            flags.append("LOW PSM")

        if id_med and v["fragpipe_id_pct"] < 0.80 * id_med:
            flags.append("LOW ID RATE")

        if q50_med and v["cas_q50_pct"] < 0.80 * q50_med:
            flags.append("LOW DE NOVO")

        if v["contaminant_psm_pct"] > 2:
            flags.append("CONTAMINANTS")

        if v["missed_cleavage_pct"] > 25:
            flags.append("MISSED CLEAVAGE")

        # Per-metric status. This is intentionally descriptive rather than
        # a hard acceptance criterion.
        result[sample] = {
            "Overall": ", ".join(flags) if flags else "PASS",
            "PSM": "WARN" if "LOW PSM" in flags else "PASS",
            "ID rate": "WARN" if "LOW ID RATE" in flags else "PASS",
            "De novo": "WARN" if "LOW DE NOVO" in flags else "PASS",
            "Missed cleavage": "WARN" if "MISSED CLEAVAGE" in flags else "PASS",
            "Contaminants": "WARN" if "CONTAMINANTS" in flags else "PASS",
        }

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="+", required=True)
    ap.add_argument("--fragpipe_suffix", default=".FPv24")
    ap.add_argument("--casanovo_suffix", default=".DN")
    ap.add_argument("--aastat_suffix", default=".AA_statm")
    ap.add_argument("--casanovo_threshold", type=float, default=DEFAULT_CASANOVO_THRESHOLD)
    args = ap.parse_args()

    dirs = [Path(x) for x in args.dirs]

    frag_dirs = {
        sample_name(p, args.fragpipe_suffix): p
        for p in dirs if p.name.endswith(args.fragpipe_suffix)
    }
    cas_dirs = {
        sample_name(p, args.casanovo_suffix): p
        for p in dirs if p.name.endswith(args.casanovo_suffix)
    }
    aa_dirs = {
        sample_name(p, args.aastat_suffix): p
        for p in dirs if p.name.endswith(args.aastat_suffix)
    }

    frag = {s: parse_fragpipe(p) for s, p in frag_dirs.items()}
    cas = {s: parse_casanovo(p) for s, p in cas_dirs.items()}
    aa = {s: parse_aastat(p) for s, p in aa_dirs.items()}

    samples = sorted(set(frag) | set(cas) | set(aa))
    summary = {}

    for sample in samples:
        f = frag.get(sample, {})
        c = cas.get(sample, {})
        a = aa.get(sample, {})

        ms2 = f.get("ms2_spectra", 0)
        fp_spec = f.get("psm_spectrum_keys", set())
        cas_spec = c.get("high_spectrum_keys", set())

        # Spectrum-level overlap. This is kept conservative because the
        # two formats use different mzTab / FragPipe spectrum identifiers.
        # If scan keys can be harmonized, compare on scan+charge.
        fp_scan_charge = set()
        for row in read_tsv(f.get("psm_file")):
            if truth(row.get("Is Decoy")) or truth(row.get("Is Contaminant")):
                continue
            scan = normalize_scan_number_from_fragpipe(row)
            charge = str(row.get("Charge") or "").strip()
            if scan:
                fp_scan_charge.add(("scan", scan, charge))

        cas_scan_charge = set()
        for row in c.get("rows", []):
            score = number(row.get("search_engine_score[1]"), float("-inf"))
            if score < args.casanovo_threshold:
                continue
            sk = casanovo_spectrum_key(row)
            if sk and sk[0] == "scan":
                cas_scan_charge.add(sk)

        spectrum_overlap = set_metrics(fp_scan_charge, cas_scan_charge)

        # Primary peptide comparison: high-confidence Casanovo only.
        cas_seq = c.get("high_sequences", set())
        fp_seq = f.get("peptides", set())
        seq_exact = set_metrics(cas_seq, fp_seq)

        cas_il = {il_key(x) for x in cas_seq}
        fp_il = {il_key(x) for x in fp_seq}
        seq_il = set_metrics(cas_il, fp_il)

        # Modified-peptidoform comparison.
        cas_pf = c.get("high_peptidoforms", set())
        fp_pf = f.get("peptidoforms", set())
        pf_shared = count_peptidoform_matches(cas_pf, fp_pf, il=False)
        pf_exact = {
            "a": len(cas_pf),
            "b": len(fp_pf),
            "shared": pf_shared,
            "a_only": max(0, len(cas_pf) - pf_shared),
            "b_only": max(0, len(fp_pf) - pf_shared),
            "union": len(cas_pf) + len(fp_pf) - pf_shared,
        }
        pf_exact["jaccard_pct"] = pct(pf_shared, pf_exact["union"])
        pf_exact["a_recovered_by_b_pct"] = pct(pf_shared, len(cas_pf))
        pf_exact["b_recovered_by_a_pct"] = pct(pf_shared, len(fp_pf))

        cas_pf_il = c.get("high_peptidoforms_il", set())
        fp_pf_il = f.get("peptidoforms_il", set())
        pf_il_shared = count_peptidoform_matches(cas_pf_il, fp_pf_il, il=True)
        pf_il = {
            "a": len(cas_pf_il),
            "b": len(fp_pf_il),
            "shared": pf_il_shared,
        }

        mod_categories = modification_difference_categories(cas_pf, fp_pf)

        q = c.get("thresholds", {})
        q50_count, q50_log_pct = q.get("0.50", (0, 0.0))
        q90_count, q90_log_pct = q.get("0.90", (0, 0.0))
        q95_count, q95_log_pct = q.get("0.95", (0, 0.0))
        q99_count, q99_log_pct = q.get("0.99", (0, 0.0))

        # If log percentages are absent, derive them from counts / sequenced.
        sequenced = c.get("sequenced", 0)
        q50_pct = q50_log_pct if q50_count else pct(q50_count, sequenced)
        q90_pct = q90_log_pct if q90_count else pct(q90_count, sequenced)
        q95_pct = q95_log_pct if q95_count else pct(q95_count, sequenced)
        q99_pct = q99_log_pct if q99_count else pct(q99_count, sequenced)

        fp_id_pct = pct(len(fp_spec), ms2)

        # Confirmation metrics use spectrum-level overlap where scan+charge
        # harmonization succeeded.
        both = spectrum_overlap["shared"]
        cas_confirm_pct = pct(both, len(cas_scan_charge))
        fp_confirm_pct = pct(both, len(fp_scan_charge))

        # Protein reproducibility is added later across samples.
        summary[sample] = {
            "total_spectra": f.get("total_spectra", 0),
            "ms2_spectra": ms2,
            "fragpipe_psm_rows": f.get("psm", 0),
            "fragpipe_psm_spectra": len(fp_spec),
            "fragpipe_id_pct": fp_id_pct,
            "peptide": f.get("peptide", 0),
            "protein": f.get("protein", 0),
            "modified_psm_pct": f.get("modified_psm_pct", 0.0),
            "contaminant_psm_pct": f.get("contaminant_psm_pct", 0.0),
            "missed_cleavage_pct": f.get("missed_cleavage_pct", 0.0),
            "casanovo_sequenced": sequenced,
            "casanovo_sequenced_pct": pct(sequenced, ms2),
            "cas_q50_count": q50_count,
            "cas_q50_pct": pct(q50_count, ms2),
            "cas_q50_of_attempted_pct": q50_pct,
            "cas_q90_pct": pct(q90_count, ms2),
            "cas_q95_pct": pct(q95_count, ms2),
            "cas_q99_pct": pct(q99_count, ms2),
            "spectrum_both": both,
            "spectrum_fragpipe_only": spectrum_overlap["a_only"],
            "spectrum_casanovo_only": spectrum_overlap["b_only"],
            "spectrum_neither": max(
                0,
                ms2 - len(fp_scan_charge) - len(cas_scan_charge) + both,
            ),
            "casanovo_confirmation_pct": cas_confirm_pct,
            "fragpipe_confirmation_pct": fp_confirm_pct,
            "sequence_casanovo_high_conf": len(cas_seq),
            "sequence_fragpipe": len(fp_seq),
            "sequence_shared_exact": seq_exact["shared"],
            "sequence_overlap_pct": pct(seq_exact["shared"], len(cas_seq)),
            "sequence_shared_il": seq_il["shared"],
            "sequence_il_overlap_pct": pct(seq_il["shared"], len(cas_il)),
            "peptidoform_casanovo_high_conf": len(cas_pf),
            "peptidoform_fragpipe": len(fp_pf),
            "peptidoform_shared": pf_exact["shared"],
            "peptidoform_overlap_pct": pct(pf_exact["shared"], len(cas_pf)),
            "peptidoform_shared_il": pf_il["shared"],
            "peptidoform_il_overlap_pct": pct(pf_il["shared"], len(cas_pf_il)),
            "same_sequence_different_modification": mod_categories[
                "same_sequence_different_modification"
            ],
            "mass_shift_count": sum(
                not x["isotope"] for x in a.get("shifts", [])
            ),
            "localized_shift_count": a.get("localized", 0),
        }

    # Protein reproducibility.
    protein_sets = {s: v["proteins"] for s, v in frag.items()}
    protein_pairs = {}
    if len(protein_sets) > 1:
        names = sorted(protein_sets)
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                shared = protein_sets[a] & protein_sets[b]
                union = protein_sets[a] | protein_sets[b]
                smaller = min(len(protein_sets[a]), len(protein_sets[b]))
                protein_pairs[f"{a} ↔ {b}"] = {
                    "Shared proteins": len(shared),
                    "Union": len(union),
                    "Jaccard %": round(pct(len(shared), len(union)), 1),
                    "Recovery of smaller run %": round(pct(len(shared), smaller), 1),
                }

    # -----------------------------------------------------------------------
    # General statistics
    # -----------------------------------------------------------------------
    headers = {
        "total_spectra": {"title": "Total spectra", "format": "{:,.0f}"},
        "ms2_spectra": {"title": "MS2 spectra", "format": "{:,.0f}"},
        "fragpipe_psm_spectra": {"title": "FragPipe PSM spectra", "format": "{:,.0f}"},
        "fragpipe_id_pct": {"title": "FragPipe ID %", "suffix": "%", "format": "{:.1f}"},
        "peptide": {"title": "Peptides", "format": "{:,.0f}"},
        "protein": {"title": "Proteins", "format": "{:,.0f}"},
        "casanovo_sequenced": {"title": "Casanovo sequenced", "format": "{:,.0f}"},
        "casanovo_sequenced_pct": {"title": "Casanovo attempted %", "suffix": "%", "format": "{:.1f}"},
        "cas_q50_pct": {"title": "Casanovo ≥0.5 %", "suffix": "%", "format": "{:.1f}"},
        "cas_q50_of_attempted_pct": {"title": "≥0.5 / attempted %", "suffix": "%", "format": "{:.1f}"},
        "casanovo_confirmation_pct": {"title": "Casanovo confirmed %", "suffix": "%", "format": "{:.1f}"},
        "fragpipe_confirmation_pct": {"title": "FragPipe confirmed %", "suffix": "%", "format": "{:.1f}"},
        "modified_psm_pct": {"title": "Modified PSMs", "suffix": "%", "format": "{:.1f}"},
        "contaminant_psm_pct": {"title": "Contaminant PSMs", "suffix": "%", "format": "{:.2f}"},
        "missed_cleavage_pct": {"title": "Missed cleavage", "suffix": "%", "format": "{:.1f}"},
    }
    if aa:
        headers["mass_shift_count"] = {"title": "AA_stat shifts", "format": "{:,.0f}"}

    write_json("opensearch_generalstats_mqc.json", {
        "id": "opensearch_generalstats",
        "section_name": "OpenSearch General Statistics",
        "description": (
            "Integrated spectrum-, peptide-, protein-, de novo- and "
            "modification-level QC. Identification percentages use MS2 spectra."
        ),
        "plot_type": "generalstats",
        "headers": headers,
        "data": summary,
    })

    # Executive scorecard.
    scorecard = {}
    for s, v in summary.items():
        scorecard[s] = {
            "MS2": v["ms2_spectra"],
            "FragPipe ID %": round(v["fragpipe_id_pct"], 1),
            "Casanovo attempted %": round(v["casanovo_sequenced_pct"], 1),
            "Casanovo ≥0.5 %": round(v["cas_q50_pct"], 1),
            "≥0.5 / attempted %": round(v["cas_q50_of_attempted_pct"], 1),
            "Casanovo confirmed %": round(v["casanovo_confirmation_pct"], 1),
            "FragPipe confirmed %": round(v["fragpipe_confirmation_pct"], 1),
            "Missed cleavage %": round(v["missed_cleavage_pct"], 1),
            "Peptides": v["peptide"],
            "Proteins": v["protein"],
        }
    mqc_table(
        "opensearch_scorecard_mqc.json",
        "OpenSearch Run Scorecard",
        "Compact run-level interpretation. Confirmation metrics are spectrum-level and use scan+charge harmonization where available.",
        scorecard,
    )

    # PTM-Shepherd.
    ranked = defaultdict(int)
    for v in frag.values():
        for m in v["modifications"]:
            if not m["isotope"]:
                ranked[m["name"]] += m["psm"]

    top_names = [
        n for n, _ in sorted(ranked.items(), key=lambda x: x[1], reverse=True)[:15]
    ]
    ptm = {
        name: {
            s: next(
                (m["pct"] for m in frag.get(s, {}).get("modifications", [])
                 if m["name"] == name),
                0.0,
            )
            for s in frag
        }
        for name in top_names
    }
    if ptm:
        mqc_bar(
            "opensearch_ptm_shepherd_mqc.json",
            "PTM-Shepherd Modification Landscape",
            "Top non-isotope PTM-Shepherd modifications",
            "% of PSMs",
            ptm,
        )

    # Charge and missed cleavage.
    if frag:
        charge = {
            s: {
                f"z={k}": n
                for k, n in sorted(v["charge"].items(), key=lambda x: str(x[0]))
            }
            for s, v in frag.items()
        }
        missed = {
            s: {
                f"{k} missed": n
                for k, n in sorted(v["missed"].items())
            }
            for s, v in frag.items()
        }
        mqc_bar(
            "opensearch_charge_mqc.json",
            "Precursor Charge Distribution",
            "FragPipe precursor charge states",
            "PSMs",
            charge,
        )
        mqc_bar(
            "opensearch_missed_cleavage_mqc.json",
            "Missed Cleavage Distribution",
            "FragPipe missed-cleavage counts",
            "PSMs",
            missed,
        )

    # Identification funnel.
    funnel = {}
    for s, v in summary.items():
        funnel[s] = {
            "Casanovo sequenced": v["casanovo_sequenced"],
            "Casanovo ≥0.5": v["cas_q50_count"],
            "FragPipe PSM spectra": v["fragpipe_psm_spectra"],
            "Unique peptides": v["peptide"],
            "Proteins": v["protein"],
        }
    mqc_bar(
        "opensearch_funnel_mqc.json",
        "Identification Funnel",
        "Integrated identification funnel",
        "Count",
        funnel,
    )

    # Casanovo confidence.
    if cas:
        confidence = {}
        for s, c in cas.items():
            q = c["thresholds"]
            confidence[s] = {
                "≥0.00": pct(q.get("0.00", (0, 0))[0], summary[s]["ms2_spectra"]),
                "≥0.50": pct(q.get("0.50", (0, 0))[0], summary[s]["ms2_spectra"]),
                "≥0.90": pct(q.get("0.90", (0, 0))[0], summary[s]["ms2_spectra"]),
                "≥0.95": pct(q.get("0.95", (0, 0))[0], summary[s]["ms2_spectra"]),
                "≥0.99": pct(q.get("0.99", (0, 0))[0], summary[s]["ms2_spectra"]),
            }

        write_json("opensearch_casanovo_confidence_mqc.json", {
            "id": "opensearch_casanovo_confidence",
            "section_name": "Casanovo Confidence",
            "description": "Percentage of MS2 spectra above Casanovo confidence thresholds.",
            "plot_type": "linegraph",
            "pconfig": {
                "id": "opensearch_casanovo_confidence_plot",
                "title": "Casanovo Confidence Thresholds",
                "ylab": "% of MS2 spectra",
                "ymin": 0,
                "ymax": 100,
            },
            "data": confidence,
        })

    # Spectrum-level comparison.
    spectrum_table = {}
    for s, v in summary.items():
        spectrum_table[s] = {
            "Input MS2": v["ms2_spectra"],
            "FragPipe PSM spectra": v["fragpipe_psm_spectra"],
            "Casanovo ≥0.5 spectra": (
                len(cas[s]["high_spectrum_keys"]) if s in cas else 0
            ),
            "Both tools": v["spectrum_both"],
            "FragPipe only": v["spectrum_fragpipe_only"],
            "Casanovo only": v["spectrum_casanovo_only"],
            "Neither": v["spectrum_neither"],
        }

    mqc_table(
        "opensearch_spectrum_comparison_mqc.json",
        "Spectrum Identification Overview",
        "Spectrum-level comparison. Both/only categories use scan number plus precursor charge when both formats expose compatible scan identifiers.",
        spectrum_table,
    )

    # Peptide-level comparison: sequence and modification-aware.
    sequence_table = {}
    for s, v in summary.items():
        sequence_table[s] = {
            "Casanovo ≥0.5 sequences": v["sequence_casanovo_high_conf"],
            "FragPipe stripped sequences": v["sequence_fragpipe"],
            "Exact AA overlap": v["sequence_shared_exact"],
            "Exact AA overlap %": round(v["sequence_overlap_pct"], 2),
            "I/L-equivalent overlap": v["sequence_shared_il"],
            "I/L-equivalent overlap %": round(v["sequence_il_overlap_pct"], 2),
            "Same sequence, different modification": v["same_sequence_different_modification"],
            "Casanovo ≥0.5 peptidoforms": v["peptidoform_casanovo_high_conf"],
            "FragPipe peptidoforms": v["peptidoform_fragpipe"],
            "Exact peptidoform overlap": v["peptidoform_shared"],
            "Exact peptidoform overlap %": round(v["peptidoform_overlap_pct"], 2),
            "I/L-equivalent peptidoform overlap": v["peptidoform_shared_il"],
            "I/L-equivalent peptidoform overlap %": round(v["peptidoform_il_overlap_pct"], 2),
        }

    mqc_table(
        "opensearch_peptide_comparison_mqc.json",
        "Casanovo ↔ FragPipe Peptide Comparison",
        (
            "Primary comparison uses Casanovo PSMs with score ≥0.50. "
            "AA-sequence comparison ignores modifications; peptidoform comparison "
            "uses position-aware modification matching. Common modification names "
            "are converted to masses and compared within 0.05 Da."
        ),
        sequence_table,
    )

    # Modification harmonization diagnostics.
    mod_diag = {}
    for s in summary:
        if s not in cas or s not in frag:
            continue

        cas_pf = cas[s]["high_peptidoforms"]
        fp_pf = frag[s]["peptidoforms"]

        cat = modification_difference_categories(cas_pf, fp_pf)

        mod_diag[s] = {
            "High-conf Casanovo peptidoforms": len(cas_pf),
            "FragPipe peptidoforms": len(fp_pf),
            "Same AA sequence": cat["same_sequence"],
            "Same AA sequence + same mods": cat["same_peptidoform"],
            "Same AA sequence + different mods": cat["same_sequence_different_modification"],
        }

    if mod_diag:
        mqc_table(
            "opensearch_modification_harmonization_mqc.json",
            "Modification Harmonization Diagnostics",
            "Shows how much apparent disagreement remains after separating amino-acid sequence disagreement from modification disagreement.",
            mod_diag,
        )

    # AA_stat annotations.
    aa_rows = []
    for sample, v in aa.items():
        for x in sorted(
            (z for z in v["shifts"] if not z["isotope"]),
            key=lambda z: z["npep"],
            reverse=True,
        )[:15]:
            ann = x["annotation"]
            label = html.escape(ann["label"])
            if ann.get("url"):
                label = (
                    f'<a href="{html.escape(ann["url"], quote=True)}" '
                    f'target="_blank" rel="noopener">{label}</a>'
                )
            aa_rows.append(
                f"<tr><td>{html.escape(sample)}</td>"
                f"<td>{x['mass']:+.4f}</td>"
                f"<td>{x['npep']:,}</td>"
                f"<td>{label}</td>"
                f"<td>{ann['match']:.0f}%</td>"
                f"<td>{source_cell('AA_stat', x['source'])}</td></tr>"
            )

    if aa_rows:
        body = (
            '<table class="mqc-table"><thead><tr>'
            "<th>Sample</th><th>Mass shift</th><th>Peptides</th>"
            "<th>Annotation</th><th>Match</th><th>Source</th>"
            "</tr></thead><tbody>"
            + "".join(aa_rows)
            + "</tbody></table>"
        )
        html_section(
            "opensearch_aastat_annotations",
            "AA_stat Mass-Shift Annotations",
            "AA_stat annotations with source files and reported match percentages.",
            body,
        )

        # Cross-sample AA_stat table.
        shifts_by_mass = defaultdict(dict)
        for sample, v in aa.items():
            for x in v["shifts"]:
                if x["isotope"]:
                    continue
                key = round(x["mass"], 4)
                shifts_by_mass[key][sample] = x

        rows = []
        for mass in sorted(shifts_by_mass):
            per = shifts_by_mass[mass]
            annotation = next(
                (x["annotation"]["label"] for x in per.values()
                 if x["annotation"]["label"] != "Unassigned"),
                "Unassigned",
            )
            cells = [f"<td>{mass:+.4f}</td>", f"<td>{html.escape(annotation)}</td>"]
            for sample in sorted(aa):
                x = per.get(sample)
                if x:
                    cells.append(
                        f"<td>{x['npep']:,}</td><td>{x['annotation']['match']:.0f}%</td>"
                    )
                else:
                    cells.append("<td>0</td><td>0%</td>")
            rows.append("<tr>" + "".join(cells) + "</tr>")

        heads = "<th>Mass shift</th><th>Annotation</th>"
        for sample in sorted(aa):
            heads += f"<th>{html.escape(sample)} peptides</th><th>Match</th>"

        body = (
            '<table class="mqc-table"><thead><tr>'
            + heads
            + "</tr></thead><tbody>"
            + "".join(rows)
            + "</tbody></table>"
        )
        html_section(
            "opensearch_aastat_crosssample",
            "AA_stat Cross-Sample Modification Comparison",
            "Cross-sample view of the most prominent non-isotope AA_stat mass shifts.",
            body,
        )

    # Protein reproducibility.
    if protein_pairs:
        mqc_table(
            "opensearch_protein_reproducibility_mqc.json",
            "Protein-Level Reproducibility",
            "Pairwise non-decoy, non-contaminant protein overlap. Recovery of smaller run is shared proteins divided by the smaller protein set.",
            protein_pairs,
        )

    # QC flags.
    qc_status = classify_metrics(summary)
    qc = {}
    for s, v in summary.items():
        qc[s] = {
            **qc_status[s],
            "MS2 spectra": v["ms2_spectra"],
            "FragPipe ID %": round(v["fragpipe_id_pct"], 1),
            "Casanovo ≥0.5 %": round(v["cas_q50_pct"], 1),
            "≥0.5 / attempted %": round(v["cas_q50_of_attempted_pct"], 1),
            "Modified PSMs %": round(v["modified_psm_pct"], 1),
            "Contaminant PSMs %": round(v["contaminant_psm_pct"], 2),
            "Missed cleavage %": round(v["missed_cleavage_pct"], 1),
        }

    mqc_table(
        "opensearch_qc_flags_mqc.json",
        "Sample QC Flags",
        "Descriptive per-metric heuristic flags for relative run review; not acceptance criteria.",
        qc,
    )

    # Provenance.
    prov = []
    source_rows = []

    for s in samples:
        if s in frag:
            f = frag[s]
            for metric, path in (
                ("PSMs", f["psm_file"]),
                ("Peptides", f["peptide_file"]),
                ("Proteins", f["protein_file"]),
                ("PTM-Shepherd modifications", f["mod_file"]),
                ("Input spectrum count", f["spectrum_file"]),
            ):
                if path:
                    rel = relative_source("FragPipe", f["directory"], path)
                    prov.append({
                        "Sample": s,
                        "Program": "FragPipe/PTM-Shepherd"
                        if metric.startswith("PTM") else "FragPipe",
                        "Metric": metric,
                        "Source": rel,
                    })

        if s in cas:
            c = cas[s]
            for metric, path in (
                ("Casanovo scores", c["log"]),
                ("De novo sequences and ProForma", c["mztab"]),
            ):
                if path:
                    prov.append({
                        "Sample": s,
                        "Program": "Casanovo",
                        "Metric": metric,
                        "Source": relative_source("Casanovo", c["directory"], path),
                    })

        if s in aa:
            a = aa[s]
            for metric, path in (
                ("Mass shifts", a["stats"]),
                ("Annotations", a["interpretations"]),
                ("Localizations", a["localizations"]),
            ):
                if path:
                    prov.append({
                        "Sample": s,
                        "Program": "AA_stat",
                        "Metric": metric,
                        "Source": relative_source("AA_stat", a["directory"], path),
                    })

            report = first(a["directory"], "report.html")
            if report:
                source_rows.append(
                    f'<tr><td>{html.escape(s)}</td><td>AA_stat</td>'
                    f'<td><a href="{html.escape(relative_source("AA_stat", a["directory"], report), quote=True)}">'
                    "AA_stat report.html</a></td></tr>"
                )

        if s in cas:
            c = cas[s]
            for path, label in (
                (c["log"], "Casanovo log"),
                (c["mztab"], "Casanovo mzTab"),
            ):
                if path:
                    source_rows.append(
                        f'<tr><td>{html.escape(s)}</td><td>Casanovo</td>'
                        f'<td><a href="{html.escape(relative_source("Casanovo", c["directory"], path), quote=True)}">'
                        f"{label}</a></td></tr>"
                    )

        if s in frag:
            f = frag[s]
            for path, label in (
                (f["psm_file"], "FragPipe psm.tsv"),
                (f["peptide_file"], "FragPipe peptide.tsv"),
                (f["mod_file"], "PTM-Shepherd global.modsummary.tsv"),
                (f["spectrum_file"], "spectrum_count.tsv"),
            ):
                if path:
                    source_rows.append(
                        f'<tr><td>{html.escape(s)}</td>'
                        "<td>FragPipe/PTM-Shepherd</td>"
                        f'<td><a href="{html.escape(relative_source("FragPipe", f["directory"], path), quote=True)}">'
                        f"{label}</a></td></tr>"
                    )

    mqc_table(
        "opensearch_provenance_mqc.json",
        "OpenSearch Data Provenance",
        "Every integrated metric is derived from one of the source files listed below.",
        {str(i + 1): r for i, r in enumerate(prov)},
    )

    if source_rows:
        body = (
            "<p>The links below point to the original files published alongside the MultiQC report.</p>"
            '<table class="mqc-table"><thead><tr>'
            "<th>Sample</th><th>Program</th><th>Source report/file</th>"
            "</tr></thead><tbody>"
            + "".join(source_rows)
            + "</tbody></table>"
        )
        html_section(
            "opensearch_source_reports",
            "Detailed Source Reports and Files",
            "Direct links to published source outputs used by the integrated report.",
            body,
        )

    # Observations.
    notes = []

    if len(summary) > 1:
        best = max(summary, key=lambda s: summary[s]["fragpipe_psm_spectra"])
        low = min(summary, key=lambda s: summary[s]["fragpipe_psm_spectra"])
        notes.append(
            f"<li><strong>{html.escape(best)}</strong> has the highest FragPipe "
            f"PSM-spectrum yield ({summary[best]['fragpipe_psm_spectra']:,}); "
            f"<strong>{html.escape(low)}</strong> has the lowest "
            f"({summary[low]['fragpipe_psm_spectra']:,}).</li>"
        )

    for s, v in summary.items():
        if v["same_sequence_different_modification"] > 0:
            notes.append(
                f"<li><strong>{html.escape(s)}</strong>: "
                f"{v['same_sequence_different_modification']:,} high-confidence "
                "Casanovo/FragPipe peptide sequences agree at the amino-acid "
                "level but not at the reported modification level.</li>"
            )

    if summary:
        q50 = [v["cas_q50_pct"] for v in summary.values()]
        notes.append(
            f"<li>Casanovo ≥0.50 rates are "
            f"{min(q50):.1f}–{max(q50):.1f}% of MS2 spectra across the available runs.</li>"
        )

    if frag:
        mods = [
            m for v in frag.values()
            for m in v["modifications"]
            if not m["isotope"]
        ]
        if mods:
            top = max(mods, key=lambda m: m["psm"])
            notes.append(
                f"<li>The strongest non-isotope PTM-Shepherd signal is "
                f"<strong>{html.escape(top['name'])}</strong> "
                f"({top['pct']:.2f}% PSMs in its source sample).</li>"
            )

    html_section(
        "opensearch_observations",
        "OpenSearch Observations",
        "Descriptive observations generated from the integrated metrics.",
        "<ul>" + ("".join(notes) or "<li>No observations available.</li>") + "</ul>"
        + "<p><em>QC flags are heuristic. Casanovo-only sequences are candidates, "
          "not automatically novel peptides.</em></p>",
    )

    # Machine-readable summary.
    fields = ["sample"] + list(next(iter(summary.values())).keys()) if summary else ["sample"]
    with open("summary.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for sample, values in summary.items():
            writer.writerow({"sample": sample, **values})

    with open("provenance.tsv", "w", encoding="utf-8", newline="") as fh:
        fields = ["Sample", "Program", "Metric", "Source"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(prov)


if __name__ == "__main__":
    main()
