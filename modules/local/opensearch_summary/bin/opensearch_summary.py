#!/usr/bin/env python3
"""Create the integrated OpenSearch MultiQC report.

The script deliberately keeps one parser per upstream result type and derives all
report sections from one compact per-sample dictionary. It accepts any mixture of
FragPipe (.FPv24 by default), Casanovo (.DN) and AA_stat (.AA_statm) directories,
so disabled optional tools simply disappear from the report.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote


# ---------- small, reusable helpers ----------

def rows(path: Path):
    if not path or not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def first(path: Path, pattern: str):
    return next(path.glob(pattern), None) if path and path.exists() else None


def count_rows(path: Path) -> int:
    return len(rows(path))


def truth(value) -> bool:
    return str(value or "").strip().lower() == "true"


def number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def pct(n, d):
    return 100.0 * n / d if d else 0.0


def sample_name(path: Path, suffixes=()) -> str:
    for suffix in suffixes:
        if suffix and path.name.endswith(suffix):
            return path.name[: -len(suffix)]
    return path.name


def strip_html(text: str) -> str:
    return html.unescape(re.sub(r"<[^>]*>", "", text or "")).strip()


def match_pct(label: str) -> float | None:
    m = re.search(r"\(([0-9.]+)%\s+match\)", label or "")
    return float(m.group(1)) if m else None


def json_write(name: str, obj: dict):
    Path(name).write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def qid(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", text).strip("_").lower() or "x"


def mqc_table(name, section, description, data):
    json_write(name, {
        "id": qid(section), "section_name": section, "description": description,
        "plot_type": "table", "pconfig": {"id": qid(section) + "_table"},
        "data": data,
    })


def mqc_bar(name, section, title, ylab, data):
    json_write(name, {
        "id": qid(section), "section_name": section,
        "plot_type": "bargraph",
        "pconfig": {"id": qid(section) + "_plot", "title": title,
                    "ylab": ylab, "cpswitch": True},
        "data": data,
    })


# ---------- upstream parsers ----------

def parse_fragpipe(directory: Path):
    psm_file = first(directory, "**/psm.tsv")
    pep_file = first(directory, "**/peptide.tsv")
    prot_file = first(directory, "**/protein.tsv")
    mod_file = first(directory, "**/global.modsummary.tsv")

    psm = rows(psm_file)
    peptides = rows(pep_file)
    proteins = rows(prot_file)

    # Unmodified peptide sequences for Casanovo ↔ FragPipe overlap.
    peptide_set = {
        (r.get("Peptide") or "").strip().upper()
        for r in peptides
        if r.get("Peptide") and not truth(r.get("Is Decoy")) and not truth(r.get("Is Contaminant"))
    }
    protein_set = {
        (r.get("Protein ID") or r.get("Protein") or "").strip()
        for r in proteins
        if r.get("Protein ID") and not truth(r.get("Is Decoy")) and not truth(r.get("Is Contaminant"))
    }

    charges = defaultdict(int)
    missed = defaultdict(int)
    modified = contaminants = missed_gt0 = 0
    for r in psm:
        charge = str(r.get("Charge") or "?").strip()
        charges[charge] += 1
        mc = int(number(r.get("Number of Missed Cleavages"), 0))
        missed[mc] += 1
        missed_gt0 += mc > 0
        modified += bool((r.get("Assigned Modifications") or "").strip() or
                         (r.get("Observed Modifications") or "").strip())
        contaminants += truth(r.get("Is Contaminant"))

    modifications = []
    if mod_file:
        for r in rows(mod_file):
            name = (r.get("Modification") or "").strip()
            if not name or name.lower() == "none":
                continue
            modifications.append({
                "name": name,
                "mass": number(r.get("Mass Shift")),
                "psm": int(number(next((v for k, v in r.items() if k.endswith("_PSMs")), 0))),
                "pct": number(next((v for k, v in r.items() if k.endswith("_percent_PSMs")), 0)),
                "isotope": "isotop" in name.lower(),
            })

    return {
        "psm": len(psm), "peptide": len(peptides), "protein": len(proteins),
        "modified_psm_pct": pct(modified, len(psm)),
        "contaminant_psm_pct": pct(contaminants, len(psm)),
        "missed_cleavage_pct": pct(missed_gt0, len(psm)),
        "charge": dict(charges), "missed": dict(missed),
        "peptides": peptide_set, "proteins": protein_set,
        "modifications": modifications,
    }


def parse_casanovo(directory: Path):
    log = first(directory, "casanovo_*.log")
    mztab = first(directory, "casanovo_*.mztab")
    text = log.read_text(encoding="utf-8", errors="replace") if log else ""

    thresholds = {}
    for threshold in ("0.00", "0.50", "0.90", "0.95", "0.99"):
        pattern = rf"([0-9]+)\s+spectra\s+\(([0-9.]+)%\)\s+scored\s+≥\s+{re.escape(threshold)}"
        m = re.search(pattern, text)
        thresholds[threshold] = (int(m.group(1)), float(m.group(2))) if m else (0, 0.0)

    m = re.search(r"Sequenced\s+([0-9]+)\s+spectra", text)
    sequenced = int(m.group(1)) if m else 0
    sequences = set()
    mztab_psm = 0
    if mztab:
        with mztab.open(encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.startswith("PSM\t"):
                    continue
                mztab_psm += 1
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 1:
                    seq = re.sub(r"[^A-Z]", "", fields[1].upper())
                    if seq:
                        sequences.add(seq)
    if not sequenced:
        sequenced = mztab_psm

    return {"sequenced": sequenced, "thresholds": thresholds, "sequences": sequences}


def parse_aastat(directory: Path):
    stats = first(directory, "aa_statistics_table.csv")
    inter = first(directory, "interpretations.json")
    loc = first(directory, "localization_statistics.csv")
    interpretations = json.loads(inter.read_text(encoding="utf-8")) if inter else {}

    shifts = []
    if stats:
        with stats.open(encoding="utf-8", errors="replace", newline="") as f:
            for r in csv.DictReader(f):
                mass = number(r.get("mass shift"), None)
                if mass is None or truth(r.get("is reference")):
                    continue
                isotope = truth(r.get("is isotope"))
                candidates = r.get("unimod accessions", "")
                key = min(interpretations, key=lambda k: abs(number(k) - mass)) if interpretations else None
                anns = []
                if key:
                    for item in interpretations.get(key, []):
                        if item.get("type") == "isotope":
                            continue
                        label = item.get("label", "")
                        clean = strip_html(label)
                        mp = match_pct(label)
                        urls = re.findall(r'href=["\']([^"\']+)', label)
                        anns.append({"label": clean, "match": mp or 0.0, "url": urls[0] if urls else ""})
                anns.sort(key=lambda x: x["match"], reverse=True)
                shifts.append({
                    "mass": mass,
                    "npep": int(number(r.get("# peptides in bin"))),
                    "isotope": isotope,
                    "candidates": candidates,
                    "annotation": anns[0] if anns else {"label": "Isotope" if isotope else "Unassigned", "match": 0.0, "url": ""},
                })

    localized = 0
    if loc:
        with loc.open(encoding="utf-8", errors="replace", newline="") as f:
            for r in csv.DictReader(f):
                if truth(r.get("is isotope")):
                    continue
                loc_text = r.get("localization", "")
                pairs = re.findall(r"['\"]([^'\"]+)['\"]\s*:\s*(\d+)", loc_text)
                localized += any(k.lower() != "non-localized" and int(v) > 0 for k, v in pairs)

    return {"shifts": shifts, "localized": localized}


def classify(summary):
    flags = {}
    for s, v in summary.items():
        f = []
        if v["psm"] and v["psm"] < 0.75 * median([x["psm"] for x in summary.values() if x["psm"]]):
            f.append("LOW PSM")
        if v["q50_pct"] and v["q50_pct"] < 0.8 * median([x["q50_pct"] for x in summary.values() if x["q50_pct"]]):
            f.append("LOW DE NOVO")
        if v["contaminant_psm_pct"] > 2:
            f.append("CONTAMINANTS")
        if v["missed_cleavage_pct"] > 25:
            f.append("MISSED CLEAVAGE")
        flags[s] = ", ".join(f) if f else "PASS"
    return flags


def median(values):
    values = sorted(values)
    if not values:
        return 0.0
    n = len(values)
    return values[n // 2] if n % 2 else (values[n // 2 - 1] + values[n // 2]) / 2


# ---------- report ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="+", required=True, help="FragPipe, Casanovo and/or AA_stat result directories")
    ap.add_argument("--fragpipe_suffix", default=".FPv24")
    ap.add_argument("--casanovo_suffix", default=".DN")
    ap.add_argument("--aastat_suffix", default=".AA_statm")
    args = ap.parse_args()

    dirs = [Path(x) for x in args.dirs]
    suffixes = (args.fragpipe_suffix, args.casanovo_suffix, args.aastat_suffix)
    frag_dirs = {sample_name(p, suffixes): p for p in dirs if p.name.endswith(args.fragpipe_suffix)}
    cas_dirs = {sample_name(p, suffixes): p for p in dirs if p.name.endswith(args.casanovo_suffix)}
    aa_dirs = {sample_name(p, suffixes): p for p in dirs if p.name.endswith(args.aastat_suffix)}

    frag = {s: parse_fragpipe(p) for s, p in frag_dirs.items()}
    cas = {s: parse_casanovo(p) for s, p in cas_dirs.items()}
    aa = {s: parse_aastat(p) for s, p in aa_dirs.items()}
    samples = sorted(set(frag) | set(cas) | set(aa))

    summary = {}
    for s in samples:
        f, c, a = frag.get(s, {}), cas.get(s, {}), aa.get(s, {})
        cseq = c.get("sequences", set())
        fpep = f.get("peptides", set())
        overlap = len(cseq & fpep)
        sequenced = c.get("sequenced", 0)
        q = c.get("thresholds", {})
        summary[s] = {
            "psm": f.get("psm", 0), "peptide": f.get("peptide", 0), "protein": f.get("protein", 0),
            "casanovo": sequenced, "q50_pct": q.get("0.50", (0, 0))[1], "q90_pct": q.get("0.90", (0, 0))[1],
            "q95_pct": q.get("0.95", (0, 0))[1], "q99_pct": q.get("0.99", (0, 0))[1],
            "denovo_per_psm": sequenced / f["psm"] if f.get("psm") else 0.0,
            "cas_unique": len(cseq), "overlap": overlap,
            "cas_only": max(0, len(cseq) - overlap), "overlap_pct": pct(overlap, len(cseq)),
            "modified_psm_pct": f.get("modified_psm_pct", 0.0),
            "contaminant_psm_pct": f.get("contaminant_psm_pct", 0.0),
            "missed_cleavage_pct": f.get("missed_cleavage_pct", 0.0),
            "mass_shift_count": sum(not x["isotope"] for x in a.get("shifts", [])),
            "localized_shift_count": a.get("localized", 0),
        }

    # 1. Compact executive table. Keep it deliberately small; the detailed
    # metrics are exposed in the sections below.
    executive = {}
    for s, v in summary.items():
        executive[s] = {
            "PSMs": v["psm"], "Peptides": v["peptide"], "Proteins": v["protein"],
            "Modified PSMs %": round(v["modified_psm_pct"], 2),
            "Contaminants %": round(v["contaminant_psm_pct"], 2),
            "Missed cleavage %": round(v["missed_cleavage_pct"], 2),
        }
        if s in cas:
            executive[s].update({"Casanovo sequences": v["casanovo"], "De novo ≥0.5 %": round(v["q50_pct"], 2)})
    mqc_table("opensearch_executive_mqc.json", "OpenSearch Executive Summary",
               "Compact cross-tool view of the dataset. Detailed QC and interpretation follow below.", executive)

    # 2. Compact general statistics.
    headers = {
        "psm": {"title": "PSMs", "format": "{:,.0f}"},
        "peptide": {"title": "Peptides", "format": "{:,.0f}"},
        "protein": {"title": "Proteins", "format": "{:,.0f}"},
        "modified_psm_pct": {"title": "Modified PSMs", "suffix": "%", "format": "{:.1f}"},
        "contaminant_psm_pct": {"title": "Contaminants", "suffix": "%", "format": "{:.2f}"},
        "missed_cleavage_pct": {"title": "Missed cleavage", "suffix": "%", "format": "{:.1f}"},
    }
    if aa:
        headers["mass_shift_count"] = {"title": "AA_stat shifts", "format": "{:,.0f}"}
    if cas:
        headers.update({
            "casanovo": {"title": "Casanovo", "format": "{:,.0f}"},
            "q50_pct": {"title": "De novo ≥0.5", "suffix": "%", "format": "{:.1f}"},
            "q90_pct": {"title": "De novo ≥0.9", "suffix": "%", "format": "{:.2f}"},
        })
    json_write("opensearch_generalstats_mqc.json", {
        "id": "opensearch_generalstats", "section_name": "OpenSearch Sample QC",
        "description": "Integrated identification, digestion, contamination, modification and de novo metrics.",
        "plot_type": "generalstats", "headers": headers, "data": summary,
    })

    # 3. PTM-Shepherd: top modifications, excluding the unmodified row.
    top = sorted({m["name"] for v in frag.values() for m in v["modifications"]})
    ranked = defaultdict(int)
    for v in frag.values():
        for m in v["modifications"]:
            ranked[m["name"]] += m["psm"]
    top = [x for x, _ in sorted(ranked.items(), key=lambda kv: kv[1], reverse=True)[:15]]
    ptm = {name: {s: next((m["pct"] for m in frag.get(s, {}).get("modifications", []) if m["name"] == name), 0.0) for s in samples if s in frag} for name in top}
    if ptm:
        mqc_bar("opensearch_ptm_shepherd_mqc.json", "PTM-Shepherd Modification Landscape", "PTM-Shepherd: Top Modification Landscape", "% of PSMs", ptm)

    # 4. Charge state distribution.
    charge = {s: {f"z={k}": n for k, n in sorted(frag[s]["charge"].items(), key=lambda x: x[0])} for s in frag}
    if charge:
        mqc_bar("opensearch_charge_mqc.json", "Precursor Charge Distribution", "FragPipe PSM Charge States", "PSMs", charge)

    # 5. Missed cleavage distribution.
    missed = {s: {f"{k} missed": n for k, n in sorted(frag[s]["missed"].items())} for s in frag}
    if missed:
        mqc_bar("opensearch_missed_cleavage_mqc.json", "Missed Cleavage Distribution", "FragPipe Missed Cleavages", "PSMs", missed)

    # 6. Identification funnel.
    if cas:
        funnel = {s: {"Casanovo sequences": v["casanovo"], "FragPipe PSMs": v["psm"], "Unique peptides": v["peptide"], "Proteins": v["protein"]} for s, v in summary.items()}
    else:
        funnel = {s: {"FragPipe PSMs": v["psm"], "Unique peptides": v["peptide"], "Proteins": v["protein"]} for s, v in summary.items()}
    mqc_bar("opensearch_funnel_mqc.json", "Identification Funnel", "OpenSearch Identification Funnel", "Count", funnel)

    # 7. Casanovo confidence.
    confidence = {}
    for s in cas:
        t = cas[s]["thresholds"]
        confidence[s] = {f"≥{k}": v[1] for k, v in t.items()}
    if confidence:
        json_write("opensearch_casanovo_confidence_mqc.json", {
            "id": "opensearch_casanovo_confidence", "section_name": "Casanovo Confidence",
            "description": "Percentage of spectra above Casanovo confidence thresholds.",
            "plot_type": "linegraph", "pconfig": {"id": "opensearch_casanovo_confidence_plot", "title": "Casanovo Confidence Thresholds", "ylab": "% of sequenced spectra", "ymin": 0, "ymax": 100},
            "data": confidence,
        })

    # 8. Exact sequence overlap.
    overlap = {s: {"Casanovo unique": v["cas_unique"], "Found in FragPipe": v["overlap"], "Casanovo-only": v["cas_only"], "Overlap %": v["overlap_pct"]} for s, v in summary.items() if v["cas_unique"]}
    if overlap:
        mqc_table("opensearch_sequence_overlap_mqc.json", "Casanovo ↔ FragPipe Sequence Overlap", "Exact overlap of unmodified peptide sequences. Casanovo-only does not imply novelty.", overlap)

    # 9. AA_stat annotation HTML.
    aa_rows = []
    for s, v in aa.items():
        for x in sorted((z for z in v["shifts"] if not z["isotope"]), key=lambda z: z["npep"], reverse=True)[:15]:
            ann = x["annotation"]
            label = html.escape(ann["label"])
            if ann.get("url"):
                label = f'<a href="{html.escape(ann["url"], quote=True)}" target="_blank">{label}</a>'
            aa_rows.append(f'<tr><td>{html.escape(s)}</td><td>{x["mass"]:+.4f}</td><td>{x["npep"]:,}</td><td>{label}</td><td>{ann.get("match",0):.0f}%</td></tr>')
    if aa_rows:
        body = """<table class="mqc-table"><thead><tr><th>Sample</th><th>Mass shift</th><th>Peptides</th><th>AA_stat annotation</th><th>Match</th></tr></thead><tbody>""" + "".join(aa_rows) + "</tbody></table>"
        Path("opensearch_aastat_annotations_mqc.html").write_text("<!-- id: opensearch_aastat_annotations\nsection_name: AA_stat Mass-Shift Annotations\ndescription: AA_stat annotations and reported Unimod match percentages. -->" + body, encoding="utf-8")

    # 10. Protein-level reproducibility: Jaccard matrix.
    protein_sets = {s: frag[s]["proteins"] for s in frag}
    if len(protein_sets) > 1:
        pairs = {}
        names = sorted(protein_sets)
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                shared = protein_sets[a] & protein_sets[b]
                union = protein_sets[a] | protein_sets[b]
                pairs[f"{a} ↔ {b}"] = {
                    "Shared proteins": len(shared),
                    "Union": len(union),
                    "Jaccard %": round(pct(len(shared), len(union)), 1),
                }
        mqc_table("opensearch_protein_reproducibility_mqc.json", "Protein-Level Reproducibility",
                   "Pairwise overlap of non-decoy, non-contaminant protein IDs.", pairs)

    # 11. Sample QC flags.
    flags = classify(summary)
    flag_table = {}
    for s in samples:
        flag_table[s] = {
            "Status": flags[s], "PSMs": summary[s]["psm"],
            "Modified PSMs %": round(summary[s]["modified_psm_pct"], 1),
            "Contaminants %": round(summary[s]["contaminant_psm_pct"], 2),
            "Missed cleavage %": round(summary[s]["missed_cleavage_pct"], 1),
        }
        if s in cas:
            flag_table[s]["De novo ≥0.5 %"] = round(summary[s]["q50_pct"], 1)
    mqc_table("opensearch_qc_flags_mqc.json", "Sample QC Flags", "Rule-based descriptive QC flags; thresholds are heuristic and intended for run review, not acceptance criteria.", flag_table)

    # 12. Compact observations.
    notes = []
    if len(summary) > 1:
        best = max(summary, key=lambda s: summary[s]["psm"])
        low = min(summary, key=lambda s: summary[s]["psm"])
        notes.append(f"<li><strong>{html.escape(best)}</strong> has the highest PSM yield ({summary[best]['psm']:,}); <strong>{html.escape(low)}</strong> has the lowest ({summary[low]['psm']:,}).</li>")
    if cas:
        q90 = [summary[s]["q90_pct"] for s in cas if s in summary]
        if q90:
            notes.append(f"<li>Casanovo ≥0.90 rates span {min(q90):.2f}–{max(q90):.2f}% across the available runs.</li>")
    if frag:
        top_mod = max((m for v in frag.values() for m in v["modifications"] if not m["isotope"]), key=lambda m: m["psm"], default=None)
        if top_mod:
            notes.append(f"<li>The largest non-isotope PTM-Shepherd signal is <strong>{html.escape(top_mod['name'])}</strong> ({top_mod['pct']:.2f}% PSMs in its source sample).</li>")
    Path("opensearch_observations_mqc.html").write_text(
        "<!-- id: opensearch_observations\nsection_name: OpenSearch Observations\ndescription: Descriptive observations generated from integrated metrics. -->"
        + "<ul>" + ("".join(notes) or "<li>No observations available.</li>") + "</ul>"
        + "<p><em>QC flags and observations are descriptive. Casanovo-only sequences are candidates, not automatically novel peptides.</em></p>", encoding="utf-8")

    # Machine-readable summary.
    fields = ["sample"] + list(next(iter(summary.values())).keys()) if summary else ["sample"]
    with open("summary.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for s, v in summary.items():
            w.writerow({"sample": s, **v})


if __name__ == "__main__":
    main()
