#!/usr/bin/env python3
"""Build cross-tool OpenSearch summary as MultiQC custom-content files."""
from __future__ import annotations
import argparse, csv, html, json, re
from pathlib import Path
from collections import Counter, defaultdict


def sample_from_path(p: Path) -> str:
    n = p.name
    for suffix in ('.FPv24', '.DN', '.AA_statm'):
        if n.endswith(suffix):
            return n[:-len(suffix)]
    return n.replace('_fragpipe_mqc.tsv','').replace('_casanovo_mqc.tsv','').replace('_aastat_mqc.tsv','')


def find_dirs(paths, suffix):
    out = {}
    for p in paths:
        p = Path(p)
        if p.is_dir():
            sid = sample_from_path(p)
            if p.name.endswith(suffix) or any(x.is_file() for x in p.rglob('*')):
                out[sid] = p
    return out


def read_mqc_tsv(path):
    with open(path, encoding='utf-8') as fh:
        lines = [x.rstrip('\n') for x in fh if x.strip() and not x.startswith('#')]
    rows = list(csv.DictReader(lines, delimiter='\t'))
    return rows[0] if rows else {}


def parse_fragpipe(frag_dirs):
    result = {}
    for p in frag_dirs.values():
        sid = sample_from_path(p)
        # The pipeline's FPv24 directory contains one sample result directory.
        pep = next(p.rglob('peptide.tsv'), None)
        psm = next(p.rglob('psm.tsv'), None)
        prot = next(p.rglob('protein.tsv'), None)
        result[sid] = {
            'psm': count_data_rows(psm),
            'peptide': count_data_rows(pep),
            'protein': count_data_rows(prot),
            'peptides': read_peptides(pep),
        }
    return result


def count_data_rows(path):
    if not path or not path.exists(): return 0
    with open(path, encoding='utf-8', errors='replace') as fh:
        return max(0, sum(1 for _ in fh) - 1)


def read_peptides(path):
    vals = set()
    if not path or not path.exists(): return vals
    with open(path, encoding='utf-8', errors='replace', newline='') as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            pep = row.get('Peptide','').strip().upper()
            if pep and row.get('Is Decoy','false').lower() != 'true':
                vals.add(pep)
    return vals


def parse_casanovo(cas_dirs):
    result = {}
    for p in cas_dirs.values():
        sid = sample_from_path(p)
        log = next(p.glob('casanovo_*.log'), None)
        mztab = next(p.glob('casanovo_*.mztab'), None)
        text = log.read_text(encoding='utf-8', errors='replace') if log else ''
        def integer(pattern):
            m = re.search(pattern, text)
            return int(m.group(1)) if m else 0
        def percent(pattern):
            m = re.search(pattern, text)
            return float(m.group(1)) if m else 0.0
        sequenced = integer(r'Sequenced\s+(\d+)\s+spectra')
        q = {}
        for threshold in ('0.00','0.50','0.90','0.95','0.99'):
            key = threshold.replace('.','')
            q[key] = integer(r'(\d+)\s+spectra\s+\([^)]*%\)\s+scored\s+≥\s+' + re.escape(threshold))
        qp = {}
        for threshold in ('0.00','0.50','0.90','0.95','0.99'):
            key = threshold.replace('.','')
            qp[key] = percent(r'(\d+)\s+spectra\s+\(([0-9.]+)%\)\s+scored\s+≥\s+' + re.escape(threshold))
            # percent() above captures the first group; use explicit second capture below.
            m = re.search(r'\d+\s+spectra\s+\(([0-9.]+)%\)\s+scored\s+≥\s+' + re.escape(threshold), text)
            qp[key] = float(m.group(1)) if m else 0.0
        seqs = set()
        n_psm = 0
        if mztab:
            with open(mztab, encoding='utf-8', errors='replace') as fh:
                for line in fh:
                    if not line.startswith('PSM\t'): continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) >= 2:
                        seqs.add(re.sub(r'[^A-Z]', '', parts[1].upper()))
                        n_psm += 1
        result[sid] = {'sequenced': sequenced, 'q': q, 'qp': qp,
                       'sequences': {x for x in seqs if x}, 'mztab_psm': n_psm}
    return result


def strip_tags(s):
    s = re.sub(r'<[^>]+>', '', s)
    return html.unescape(s)


def parse_match_percent(label):
    m = re.search(r'\(([0-9.]+)%\s+match\)', label)
    return float(m.group(1)) if m else -1.0


def parse_aastat(aa_dirs):
    result = {}
    for p in aa_dirs.values():
        sid = sample_from_path(p)
        stats = next(p.glob('aa_statistics_table.csv'), None)
        inter = next(p.glob('interpretations.json'), None)
        loc = next(p.glob('localization_statistics.csv'), None)
        shifts = []
        if stats:
            with open(stats, encoding='utf-8', errors='replace', newline='') as fh:
                for row in csv.DictReader(fh):
                    try: mass = float(row['mass shift'])
                    except (ValueError, KeyError): continue
                    if abs(mass) < 1e-9 or str(row.get('is reference','')).lower() == 'true': continue
                    isotope = str(row.get('is isotope','')).lower() == 'true'
                    npep = int(float(row.get('# peptides in bin', 0) or 0))
                    shifts.append({'mass': mass, 'npep': npep, 'isotope': isotope})
        interpretations = json.loads(inter.read_text(encoding='utf-8')) if inter else {}
        for s in shifts:
            key = min(interpretations, key=lambda k: abs(float(k)-s['mass'])) if interpretations else None
            labels = []
            if key is not None:
                for item in interpretations.get(key, []):
                    if item.get('type') != 'isotope':
                        labels.append((parse_match_percent(item.get('label','')), strip_tags(item.get('label',''))))
            labels.sort(reverse=True)
            s['interpretation'] = labels[0][1] if labels else ('isotope' if s['isotope'] else 'Unassigned')
        localized_shifts = 0
        if loc:
            with open(loc, encoding='utf-8', errors='replace', newline='') as fh:
                for row in csv.DictReader(fh):
                    try: mass = float(row['mass shift'])
                    except (ValueError, KeyError): continue
                    if abs(mass) < 1e-9 or str(row.get('is isotope','')).lower() == 'true': continue
                    d = row.get('localization','')
                    if d and d not in ('{}','set()'):
                        try:
                            # The value is a Python-dict-like string. Count any entry except non-localized.
                            pairs = re.findall(r"'([^']+)':\s*(\d+)", d)
                            if any(k != 'non-localized' and int(v) > 0 for k,v in pairs):
                                localized_shifts += 1
                        except Exception: pass
        result[sid] = {'shifts': shifts, 'localized_shifts': localized_shifts}
    return result


def fmt(x):
    if isinstance(x, float): return f'{x:.2f}'
    return str(x)


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding='utf-8')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fragpipe', nargs='+', required=True)
    ap.add_argument('--casanovo', nargs='+', required=True)
    ap.add_argument('--aastat', nargs='+', required=True)
    args = ap.parse_args()
    frag_dirs = find_dirs(args.fragpipe, '.FPv24')
    cas_dirs = find_dirs(args.casanovo, '.DN')
    aa_dirs = find_dirs(args.aastat, '.AA_statm')
    frag, cas, aa = parse_fragpipe(frag_dirs), parse_casanovo(cas_dirs), parse_aastat(aa_dirs)
    samples = sorted(set(frag) | set(cas) | set(aa))

    summary = {}
    for s in samples:
        f, c, a = frag.get(s, {}), cas.get(s, {}), aa.get(s, {})
        fpep = f.get('peptides', set())
        cseq = c.get('sequences', set())
        overlap = len(fpep & cseq)
        cunq = len(cseq)
        summary[s] = {
            'psm': f.get('psm',0), 'peptide': f.get('peptide',0), 'protein': f.get('protein',0),
            'casanovo': c.get('sequenced',0), 'q50': c.get('q',{}).get('050',0), 'q90': c.get('q',{}).get('090',0),
            'q95': c.get('q',{}).get('095',0), 'q99': c.get('q',{}).get('099',0),
            'q00_pct': c.get('qp',{}).get('000',0.0), 'q50_pct': c.get('qp',{}).get('050',0.0), 'q90_pct': c.get('qp',{}).get('090',0.0),
            'cas_unique': cunq, 'overlap_unique': overlap,
            'casanovo_only_unique': max(0,cunq-overlap),
            'overlap_pct': (100*overlap/cunq) if cunq else 0.0,
            'denovo_per_psm': (c.get('sequenced',0)/f.get('psm',1)) if f.get('psm',0) else 0.0,
            'mass_shift_count': sum(1 for x in a.get('shifts',[]) if not x['isotope']),
            'localized_shift_count': a.get('localized_shifts',0),
        }

    # Main integrated table.
    write_json(Path('opensearch_executive_mqc.json'), {
      'id':'opensearch_executive','section_name':'OpenSearch Executive Summary',
      'description':'Integrated FragPipe, Casanovo and AA_stat metrics. Casanovo/FragPipe overlap is exact sequence overlap after removing modification annotations.',
      'plot_type':'table','pconfig':{'id':'opensearch_executive_table','title':'OpenSearch: Integrated Results'},
      'data': {s:{k:v for k,v in summary[s].items() if k not in ('cas_unique','overlap_unique','casanovo_only_unique')} for s in samples}
    })

    # General statistics.
    headers = [
      {'psm':{'title':'PSMs','description':'FragPipe PSM rows','min':0,'format':'{:,.0f}'}},
      {'peptide':{'title':'Peptides','min':0,'format':'{:,.0f}'}},
      {'protein':{'title':'Proteins','min':0,'format':'{:,.0f}'}},
      {'casanovo':{'title':'Casanovo','description':'Spectra sequenced by Casanovo','min':0,'format':'{:,.0f}'}},
      {'q50_pct':{'title':'Casanovo ≥0.5','suffix':'%','min':0,'max':100,'format':'{:.1f}'}},
      {'q90_pct':{'title':'Casanovo ≥0.9','suffix':'%','min':0,'max':100,'format':'{:.2f}'}},
      {'denovo_per_psm':{'title':'De novo / PSM','min':0,'format':'{:.2f}'}},
      {'overlap_pct':{'title':'Sequence overlap','suffix':'%','min':0,'max':100,'format':'{:.1f}'}},
      {'mass_shift_count':{'title':'Mass shifts','min':0,'format':'{:,.0f}'}},
      {'localized_shift_count':{'title':'Localized shifts','min':0,'format':'{:,.0f}'}},
    ]
    write_json(Path('opensearch_generalstats_mqc.json'), {
      'id':'opensearch_generalstats','plot_type':'generalstats','headers':headers,
      'data':summary
    })

    # Identification funnel.
    funnel={s:{'Casanovo sequences':v['casanovo'],'FragPipe PSMs':v['psm'],'Unique peptides':v['peptide'],'Proteins':v['protein']} for s,v in summary.items()}
    write_json(Path('opensearch_funnel_mqc.json'), {'id':'opensearch_funnel','section_name':'Identification Funnel','description':'Comparison of de novo sequence output with database-search identifications.','plot_type':'bargraph','pconfig':{'id':'opensearch_funnel_plot','title':'OpenSearch: Identification Funnel','ylab':'Count','cpswitch':True},'data':funnel})

    # Casanovo confidence thresholds.
    conf={}
    for s,v in summary.items():
        conf[s] = [['≥0.00',v['q00_pct']], ['≥0.50',v['q50_pct']],['≥0.90',v['q90_pct']],['≥0.95',100*v['q95']/v['casanovo'] if v['casanovo'] else 0],['≥0.99',100*v['q99']/v['casanovo'] if v['casanovo'] else 0]]
    # Replace ≥0.00 with the actual log percentage if available: currently inferred as 0 if absent.
    # The q50/q90 values are authoritative from Casanovo's log.
    write_json(Path('opensearch_casanovo_confidence_mqc.json'), {'id':'opensearch_casanovo_confidence','section_name':'Casanovo Confidence','description':'Fraction of sequenced spectra above Casanovo confidence thresholds, as reported in the Casanovo log.','plot_type':'linegraph','pconfig':{'id':'opensearch_casanovo_confidence_plot','title':'OpenSearch: Casanovo Confidence Thresholds','ylab':'Percent of sequenced spectra','ymin':0,'ymax':100},'data':conf})

    # Sequence overlap table.
    overlap={s:{'Casanovo unique sequences':v['cas_unique'],'Found in FragPipe peptides':v['overlap_unique'],'Casanovo-only sequences':v['casanovo_only_unique'],'Overlap %':v['overlap_pct']} for s,v in summary.items()}
    write_json(Path('opensearch_sequence_overlap_mqc.json'), {'id':'opensearch_sequence_overlap','section_name':'Casanovo ↔ FragPipe Sequence Overlap','description':'Exact overlap of unique unmodified peptide sequences. Casanovo modification annotations are removed for this comparison, so non-overlap is not evidence of novelty.','plot_type':'table','pconfig':{'id':'opensearch_sequence_overlap_table','title':'OpenSearch: De Novo / Database Sequence Overlap'},'data':overlap})

    # Mass-shift landscape, top non-isotope shifts per sample.
    mass_data={}
    mass_table={}
    for s,a in aa.items():
        rows=sorted([x for x in a.get('shifts',[]) if not x['isotope']], key=lambda x:x['npep'], reverse=True)[:10]
        mass_data[s]={f"{x['mass']:+.3f} ({x['interpretation']})":x['npep'] for x in rows}
        for x in rows:
            mass_table.setdefault(f"{x['mass']:+.3f}",{})[s]=x['npep']
    write_json(Path('opensearch_mass_shifts_mqc.json'), {'id':'opensearch_mass_shifts','section_name':'AA_stat Mass-Shift Landscape','description':'Top non-isotope mass shifts by peptide count. Interpretation is the highest reported AA_stat match percentage, not a definitive identification.','plot_type':'bargraph','pconfig':{'id':'opensearch_mass_shifts_plot','title':'OpenSearch: AA_stat Mass Shifts','ylab':'Peptides','cpswitch':True},'data':mass_data})

    # Human-readable observations as custom HTML.
    notes=[]
    if summary:
        best=max(summary,key=lambda s:summary[s]['psm'])
        notes.append(f'<li><strong>{html.escape(best)}</strong> has the highest FragPipe PSM count ({summary[best]["psm"]:,}).</li>')
        ratios=[v['denovo_per_psm'] for v in summary.values() if v['psm']]
        if ratios and max(ratios)-min(ratios) > 0.5:
            s=max(summary,key=lambda x:summary[x]['denovo_per_psm'])
            notes.append(f'<li><strong>{html.escape(s)}</strong> has the highest Casanovo-to-PSM ratio ({summary[s]["denovo_per_psm"]:.2f}).</li>')
        q90s=[v['q90_pct'] for v in summary.values()]
        if q90s and max(q90s)-min(q90s) < 0.5:
            notes.append(f'<li>Casanovo ≥0.90 rates are consistent across samples ({min(q90s):.2f}–{max(q90s):.2f}%).</li>')
        low=min(summary,key=lambda s:summary[s]['q50_pct'])
        notes.append(f'<li><strong>{html.escape(low)}</strong> has the lowest Casanovo ≥0.50 rate ({summary[low]["q50_pct"]:.2f}%).</li>')
    body=''.join(notes) or '<li>No integrated observations could be generated.</li>'
    html_content='<!--\nid: opensearch_observations\nsection_name: OpenSearch Observations\ndescription: Automatically generated observations from the integrated metrics.\n--><div class="mqc-opensearch-observations"><ul>'+body+'</ul><p><em>These are descriptive QC observations, not biological conclusions. Casanovo-only sequences require validation.</em></p></div>'
    Path('opensearch_observations_mqc.html').write_text(html_content, encoding='utf-8')

    with open('summary.tsv','w',encoding='utf-8',newline='') as fh:
        fields=['sample']+list(next(iter(summary.values())).keys()) if summary else ['sample']
        w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t'); w.writeheader()
        for s,v in summary.items(): w.writerow({'sample':s,**v})

if __name__=='__main__': main()
