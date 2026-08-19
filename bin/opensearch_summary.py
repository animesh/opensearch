#!/usr/bin/env python3
"""Build a compact, provenance-aware OpenSearch MultiQC report.

All report sections are derived from one per-sample data model.  Each metric
retains its source program and source file so the integrated report can be
traced back to FragPipe, PTM-Shepherd, Casanovo, or AA_stat outputs.
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


def read_tsv(path: Path | None):
    if not path or not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def first(path: Path | None, pattern: str):
    return next(path.glob(pattern), None) if path and path.exists() else None


def truth(value) -> bool:
    return str(value or "").strip().lower() == "true"


def number(value, default=0.0):
    try:
        return float(value)
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
        "pconfig": {"id": slug(section) + "_plot", "title": title, "ylab": ylab, "cpswitch": True},
        "data": data,
    })


def relative_source(tool: str, directory: Path, file: Path | None = None) -> str:
    """Path from results/multiqc/ to the published source result."""
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


def parse_fragpipe(directory: Path):
    psm_file = first(directory, "**/psm.tsv")
    pep_file = first(directory, "**/peptide.tsv")
    prot_file = first(directory, "**/protein.tsv")
    mod_file = first(directory, "**/global.modsummary.tsv")
    psm = read_tsv(psm_file)
    peptides = read_tsv(pep_file)
    proteins = read_tsv(prot_file)

    peptide_set = {
        (r.get("Peptide") or "").strip().upper()
        for r in peptides
        if r.get("Peptide") and not truth(r.get("Is Decoy")) and not truth(r.get("Is Contaminant"))
    }
    protein_set = {
        (r.get("Protein ID") or r.get("Protein") or "").strip()
        for r in proteins
        if (r.get("Protein ID") or r.get("Protein"))
        and not truth(r.get("Is Decoy")) and not truth(r.get("Is Contaminant"))
    }

    charge = defaultdict(int)
    missed = defaultdict(int)
    modified = contaminants = missed_gt0 = 0
    for r in psm:
        charge[str(r.get("Charge") or "?").strip()] += 1
        mc = int(number(r.get("Number of Missed Cleavages"), 0))
        missed[mc] += 1
        missed_gt0 += mc > 0
        modified += bool((r.get("Assigned Modifications") or "").strip() or (r.get("Observed Modifications") or "").strip())
        contaminants += truth(r.get("Is Contaminant"))

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
                "psm": int(number(r.get(psm_key), 0)),
                "pct": number(r.get(pct_key), 0),
                "isotope": "isotop" in name.lower(),
                "source": relative_source("FragPipe", directory, mod_file),
            })

    return {
        "program": "FragPipe",
        "directory": directory,
        "psm_file": psm_file,
        "peptide_file": pep_file,
        "protein_file": prot_file,
        "mod_file": mod_file,
        "psm": len(psm),
        "peptide": len(peptides),
        "protein": len(proteins),
        "modified_psm_pct": pct(modified, len(psm)),
        "contaminant_psm_pct": pct(contaminants, len(psm)),
        "missed_cleavage_pct": pct(missed_gt0, len(psm)),
        "charge": dict(charge),
        "missed": dict(missed),
        "peptides": peptide_set,
        "proteins": protein_set,
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
        with mztab.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith("PSM\t"):
                    continue
                mztab_psm += 1
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 1:
                    seq = re.sub(r"[^A-Z]", "", fields[1].upper())
                    if seq:
                        sequences.add(seq)
    sequenced = sequenced or mztab_psm
    return {
        "program": "Casanovo",
        "directory": directory,
        "log": log,
        "mztab": mztab,
        "sequenced": sequenced,
        "thresholds": thresholds,
        "sequences": sequences,
    }


def parse_aastat(directory: Path):
    stats = first(directory, "aa_statistics_table.csv")
    inter = first(directory, "interpretations.json")
    loc = first(directory, "localization_statistics.csv")
    interpretations = json.loads(inter.read_text(encoding="utf-8")) if inter else {}
    shifts = []
    if stats:
        with stats.open(encoding="utf-8", errors="replace", newline="") as fh:
            for r in csv.DictReader(fh):
                mass = number(r.get("mass shift"), None)
                if mass is None or truth(r.get("is reference")) or number(r.get("# peptides in bin"), 0) <= 0:
                    continue
                isotope = truth(r.get("is isotope"))
                key = min(interpretations, key=lambda k: abs(number(k) - mass)) if interpretations else None
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
                    "npep": int(number(r.get("# peptides in bin"), 0)),
                    "isotope": isotope,
                    "annotation": annotations[0] if annotations else {
                        "label": "Isotope" if isotope else "Unassigned", "match": 0.0, "url": ""
                    },
                    "source": relative_source("AA_stat", directory, stats),
                })

    localized = 0
    if loc:
        with loc.open(encoding="utf-8", errors="replace", newline="") as fh:
            for r in csv.DictReader(fh):
                if truth(r.get("is isotope")) or number(r.get("# peptides in bin"), 0) <= 0:
                    continue
                loc_text = r.get("localization", "")
                pairs = re.findall(r"['\"]([^'\"]+)['\"]\s*:\s*(\d+)", loc_text)
                localized += any(k.lower() != "non-localized" and int(v) > 0 for k, v in pairs)

    return {
        "program": "AA_stat",
        "directory": directory,
        "stats": stats,
        "interpretations": inter,
        "localizations": loc,
        "shifts": shifts,
        "localized": localized,
    }


def classify(summary):
    psm_med = median([v["psm"] for v in summary.values() if v["psm"]])
    q50_med = median([v["q50_pct"] for v in summary.values() if v["q50_pct"]])
    flags = {}
    for sample, v in summary.items():
        f = []
        if psm_med and v["psm"] < 0.75 * psm_med:
            f.append("LOW PSM")
        if q50_med and v["q50_pct"] < 0.8 * q50_med:
            f.append("LOW DE NOVO")
        if v["contaminant_psm_pct"] > 2:
            f.append("CONTAMINANTS")
        if v["missed_cleavage_pct"] > 25:
            f.append("MISSED CLEAVAGE")
        flags[sample] = ", ".join(f) if f else "PASS"
    return flags


def html_section(section_id, section_name, description, body):
    """Valid MultiQC custom HTML metadata followed by HTML content."""
    meta = "\n".join([
        "<!--",
        f"id: {section_id}",
        f"section_name: {section_name}",
        f"description: {description}",
        "-->",
    ])
    Path(f"{section_id}_mqc.html").write_text(meta + body, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="+", required=True)
    ap.add_argument("--fragpipe_suffix", default=".FPv24")
    ap.add_argument("--casanovo_suffix", default=".DN")
    ap.add_argument("--aastat_suffix", default=".AA_statm")
    args = ap.parse_args()

    dirs = [Path(x) for x in args.dirs]
    frag_dirs = {sample_name(p, args.fragpipe_suffix): p for p in dirs if p.name.endswith(args.fragpipe_suffix)}
    cas_dirs = {sample_name(p, args.casanovo_suffix): p for p in dirs if p.name.endswith(args.casanovo_suffix)}
    aa_dirs = {sample_name(p, args.aastat_suffix): p for p in dirs if p.name.endswith(args.aastat_suffix)}
    frag = {s: parse_fragpipe(p) for s, p in frag_dirs.items()}
    cas = {s: parse_casanovo(p) for s, p in cas_dirs.items()}
    aa = {s: parse_aastat(p) for s, p in aa_dirs.items()}
    samples = sorted(set(frag) | set(cas) | set(aa))

    summary = {}
    for sample in samples:
        f, c, a = frag.get(sample, {}), cas.get(sample, {}), aa.get(sample, {})
        cseq, fpep = c.get("sequences", set()), f.get("peptides", set())
        overlap = len(cseq & fpep)
        q = c.get("thresholds", {})
        summary[sample] = {
            "psm": f.get("psm", 0),
            "peptide": f.get("peptide", 0),
            "protein": f.get("protein", 0),
            "modified_psm_pct": f.get("modified_psm_pct", 0.0),
            "contaminant_psm_pct": f.get("contaminant_psm_pct", 0.0),
            "missed_cleavage_pct": f.get("missed_cleavage_pct", 0.0),
            "casanovo": c.get("sequenced", 0),
            "q50_pct": q.get("0.50", (0, 0.0))[1],
            "q90_pct": q.get("0.90", (0, 0.0))[1],
            "q95_pct": q.get("0.95", (0, 0.0))[1],
            "q99_pct": q.get("0.99", (0, 0.0))[1],
            "denovo_per_psm": c.get("sequenced", 0) / f["psm"] if f.get("psm") else 0.0,
            "cas_unique": len(cseq),
            "overlap": overlap,
            "cas_only": max(0, len(cseq) - overlap),
            "overlap_pct": pct(overlap, len(cseq)),
            "mass_shift_count": sum(not x["isotope"] for x in a.get("shifts", [])),
            "localized_shift_count": a.get("localized", 0),
        }

    # Executive summary.
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
    mqc_table("opensearch_executive_summary_mqc.json", "OpenSearch Executive Summary", "Integrated overview. Detailed sections retain program and source-file provenance.", executive)

    headers = {
        "psm": {"title": "PSMs", "format": "{:,.0f}"},
        "peptide": {"title": "Peptides", "format": "{:,.0f}"},
        "protein": {"title": "Proteins", "format": "{:,.0f}"},
        "modified_psm_pct": {"title": "Modified PSMs", "suffix": "%", "format": "{:.1f}"},
        "contaminant_psm_pct": {"title": "Contaminant PSMs", "suffix": "%", "format": "{:.2f}"},
        "missed_cleavage_pct": {"title": "Missed cleavage", "suffix": "%", "format": "{:.1f}"},
    }
    if cas:
        headers.update({
            "casanovo": {"title": "Casanovo", "format": "{:,.0f}"},
            "q50_pct": {"title": "De novo ≥0.5", "suffix": "%", "format": "{:.1f}"},
            "q90_pct": {"title": "De novo ≥0.9", "suffix": "%", "format": "{:.2f}"},
        })
    if aa:
        headers["mass_shift_count"] = {"title": "AA_stat shifts", "format": "{:,.0f}"}
    write_json("opensearch_generalstats_mqc.json", {
        "id": "opensearch_generalstats", "section_name": "OpenSearch Sample QC",
        "description": "Integrated FragPipe, PTM-Shepherd, Casanovo and AA_stat metrics.",
        "plot_type": "generalstats", "headers": headers, "data": summary,
    })

    # PTM-Shepherd modification landscape.
    ranked = defaultdict(int)
    for v in frag.values():
        for m in v["modifications"]:
            if not m["isotope"]:
                ranked[m["name"]] += m["psm"]
    top_names = [n for n, _ in sorted(ranked.items(), key=lambda x: x[1], reverse=True)[:15]]
    ptm = {
        name: {s: next((m["pct"] for m in frag.get(s, {}).get("modifications", []) if m["name"] == name), 0.0) for s in frag}
        for name in top_names
    }
    if ptm:
        mqc_bar("opensearch_ptm_shepherd_mqc.json", "PTM-Shepherd Modification Landscape", "Top non-isotope PTM-Shepherd modifications", "% of PSMs", ptm)

    # Charge and missed-cleavage distributions.
    if frag:
        charge = {s: {f"z={k}": n for k, n in sorted(v["charge"].items(), key=lambda x: str(x[0]))} for s, v in frag.items()}
        missed = {s: {f"{k} missed": n for k, n in sorted(v["missed"].items())} for s, v in frag.items()}
        mqc_bar("opensearch_charge_mqc.json", "Precursor Charge Distribution", "FragPipe precursor charge states", "PSMs", charge)
        mqc_bar("opensearch_missed_cleavage_mqc.json", "Missed Cleavage Distribution", "FragPipe missed-cleavage counts", "PSMs", missed)

    # Identification funnel.
    funnel = {}
    for s, v in summary.items():
        funnel[s] = ({"Casanovo sequences": v["casanovo"], "FragPipe PSMs": v["psm"], "Unique peptides": v["peptide"], "Proteins": v["protein"]}
                     if cas else {"FragPipe PSMs": v["psm"], "Unique peptides": v["peptide"], "Proteins": v["protein"]})
    mqc_bar("opensearch_funnel_mqc.json", "Identification Funnel", "Integrated identification funnel", "Count", funnel)

    # Casanovo confidence.
    if cas:
        confidence = {s: {f"≥{k}": v[1] for k, v in c["thresholds"].items()} for s, c in cas.items()}
        write_json("opensearch_casanovo_confidence_mqc.json", {
            "id": "opensearch_casanovo_confidence", "section_name": "Casanovo Confidence",
            "description": "Percentage of spectra above Casanovo confidence thresholds.", "plot_type": "linegraph",
            "pconfig": {"id": "opensearch_casanovo_confidence_plot", "title": "Casanovo Confidence Thresholds", "ylab": "% of sequenced spectra", "ymin": 0, "ymax": 100},
            "data": confidence,
        })

    # Sequence overlap.
    if cas:
        overlap_table = {
            s: {"Casanovo sequences": v["cas_unique"], "Found in FragPipe": v["overlap"], "Casanovo-only candidates": v["cas_only"], "Overlap %": round(v["overlap_pct"], 2)}
            for s, v in summary.items() if v["cas_unique"]
        }
        if overlap_table:
            mqc_table("opensearch_sequence_overlap_mqc.json", "Casanovo ↔ FragPipe Sequence Overlap", "Exact sequence overlap. Casanovo-only sequences are candidates, not automatically novel peptides.", overlap_table)

    # AA_stat annotations.
    aa_rows = []
    for sample, v in aa.items():
        for x in sorted((z for z in v["shifts"] if not z["isotope"]), key=lambda z: z["npep"], reverse=True)[:15]:
            ann = x["annotation"]
            label = html.escape(ann["label"])
            if ann.get("url"):
                label = f'<a href="{html.escape(ann["url"], quote=True)}" target="_blank" rel="noopener">{label}</a>'
            aa_rows.append(f'<tr><td>{html.escape(sample)}</td><td>{x["mass"]:+.4f}</td><td>{x["npep"]:,}</td><td>{label}</td><td>{ann["match"]:.0f}%</td><td>{source_cell("AA_stat", x["source"])}</td></tr>')
    if aa_rows:
        body = ('<table class="mqc-table"><thead><tr><th>Sample</th><th>Mass shift</th><th>Peptides</th><th>Annotation</th><th>Match</th><th>Source</th></tr></thead><tbody>'
                + "".join(aa_rows) + "</tbody></table>")
        html_section("opensearch_aastat_annotations", "AA_stat Mass-Shift Annotations", "AA_stat annotations with source files and reported match percentages.", body)

    # Protein-level reproducibility.
    protein_sets = {s: v["proteins"] for s, v in frag.items()}
    if len(protein_sets) > 1:
        pairs = {}
        names = sorted(protein_sets)
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                shared = protein_sets[a] & protein_sets[b]
                union = protein_sets[a] | protein_sets[b]
                pairs[f"{a} ↔ {b}"] = {"Shared proteins": len(shared), "Union": len(union), "Jaccard %": round(pct(len(shared), len(union)), 1)}
        mqc_table("opensearch_protein_reproducibility_mqc.json", "Protein-Level Reproducibility", "Pairwise non-decoy, non-contaminant protein overlap derived from FragPipe protein.tsv.", pairs)

    # QC flags.
    flags = classify(summary)
    qc = {}
    for s, v in summary.items():
        qc[s] = {
            "Status": flags[s], "PSMs": v["psm"], "Modified PSMs %": round(v["modified_psm_pct"], 1),
            "Contaminant PSMs %": round(v["contaminant_psm_pct"], 2), "Missed cleavage %": round(v["missed_cleavage_pct"], 1),
        }
        if s in cas:
            qc[s]["De novo ≥0.5 %"] = round(v["q50_pct"], 1)
    mqc_table("opensearch_qc_flags_mqc.json", "Sample QC Flags", "Descriptive heuristic flags for relative run review; not acceptance criteria.", qc)

    # Provenance map and original-source links.
    prov = []
    source_rows = []
    for s in samples:
        if s in frag:
            f = frag[s]
            for metric, path in (("PSMs", f["psm_file"]), ("Peptides", f["peptide_file"]), ("Proteins", f["protein_file"]), ("PTM-Shepherd modifications", f["mod_file"])):
                if path:
                    rel = relative_source("FragPipe", f["directory"], path)
                    prov.append({"Sample": s, "Program": "FragPipe/PTM-Shepherd" if metric.startswith("PTM") else "FragPipe", "Metric": metric, "Source": rel})
        if s in cas:
            c = cas[s]
            for metric, path in (("Casanovo scores", c["log"]), ("De novo sequences", c["mztab"])):
                if path:
                    prov.append({"Sample": s, "Program": "Casanovo", "Metric": metric, "Source": relative_source("Casanovo", c["directory"], path)})
        if s in aa:
            a = aa[s]
            for metric, path in (("Mass shifts", a["stats"]), ("Annotations", a["interpretations"]), ("Localizations", a["localizations"])):
                if path:
                    prov.append({"Sample": s, "Program": "AA_stat", "Metric": metric, "Source": relative_source("AA_stat", a["directory"], path)})

            report = first(a["directory"], "report.html")
            if report:
                source_rows.append(f'<tr><td>{html.escape(s)}</td><td>AA_stat</td><td><a href="{html.escape(relative_source("AA_stat", a["directory"], report), quote=True)}">AA_stat report.html</a></td></tr>')
        if s in cas:
            c = cas[s]
            for path, label in ((c["log"], "Casanovo log"), (c["mztab"], "Casanovo mzTab")):
                if path:
                    source_rows.append(f'<tr><td>{html.escape(s)}</td><td>Casanovo</td><td><a href="{html.escape(relative_source("Casanovo", c["directory"], path), quote=True)}">{label}</a></td></tr>')
        if s in frag:
            f = frag[s]
            for path, label in ((f["psm_file"], "FragPipe psm.tsv"), (f["mod_file"], "PTM-Shepherd global.modsummary.tsv")):
                if path:
                    source_rows.append(f'<tr><td>{html.escape(s)}</td><td>FragPipe/PTM-Shepherd</td><td><a href="{html.escape(relative_source("FragPipe", f["directory"], path), quote=True)}">{label}</a></td></tr>')

    mqc_table("opensearch_provenance_mqc.json", "OpenSearch Data Provenance", "Every integrated metric is derived from one of the source files listed below.", {str(i + 1): r for i, r in enumerate(prov)})
    if source_rows:
        body = ('<p>The links below point to the original files published alongside the MultiQC report.</p>'
                '<table class="mqc-table"><thead><tr><th>Sample</th><th>Program</th><th>Source report/file</th></tr></thead><tbody>'
                + "".join(source_rows) + "</tbody></table>")
        html_section("opensearch_source_reports", "Detailed Source Reports and Files", "Direct links to published source outputs used by the integrated report.", body)

    # Observations.
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
        mods = [m for v in frag.values() for m in v["modifications"] if not m["isotope"]]
        if mods:
            top = max(mods, key=lambda m: m["psm"])
            notes.append(f"<li>The strongest non-isotope PTM-Shepherd signal is <strong>{html.escape(top['name'])}</strong> ({top['pct']:.2f}% PSMs in its source sample).</li>")
    html_section("opensearch_observations", "OpenSearch Observations", "Descriptive observations generated from the integrated metrics.",
                 "<ul>" + ("".join(notes) or "<li>No observations available.</li>") + "</ul>"
                 + "<p><em>QC flags are heuristic. Casanovo-only sequences are candidates, not automatically novel peptides.</em></p>")

    # Machine-readable summary.
    fields = ["sample"] + list(next(iter(summary.values())).keys()) if summary else ["sample"]
    with open("summary.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for sample, values in summary.items():
            writer.writerow({"sample": sample, **values})

    # Machine-readable provenance, useful outside MultiQC too.
    with open("provenance.tsv", "w", encoding="utf-8", newline="") as fh:
        fields = ["Sample", "Program", "Metric", "Source"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(prov)


if __name__ == "__main__":
    main()
