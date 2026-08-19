#!/usr/bin/env python3
import json, tempfile, unittest
from pathlib import Path
from subprocess import run

SCRIPT = Path(__file__).resolve().parents[2] / 'bin' / 'opensearch_summary.py'

class TestOpenSearchSummary(unittest.TestCase):
    def test_integrated_summary(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            sample = 'S1'
            fp = root / f'{sample}.FPv24' / sample
            cn = root / f'{sample}.DN'
            aa = root / f'{sample}.AA_statm'
            fp.mkdir(parents=True); cn.mkdir(); aa.mkdir()

            (fp/'peptide.tsv').write_text('Peptide\tIs Decoy\nPEPTIDE\tfalse\nOTHER\tfalse\n')
            (fp/'psm.tsv').write_text('x\n1\n2\n')
            (fp/'protein.tsv').write_text('x\nP1\n')
            (cn/'casanovo_test.log').write_text(
                'Sequenced 100 spectra\n'
                '50 spectra (50.00%) scored ≥ 0.00\n'
                '25 spectra (25.00%) scored ≥ 0.50\n'
                '5 spectra (5.00%) scored ≥ 0.90\n'
                '0 spectra (0.00%) scored ≥ 0.95\n'
                '0 spectra (0.00%) scored ≥ 0.99\n')
            (cn/'casanovo_test.mztab').write_text(
                'PSM\tPEPTIDE\t1\tnull\tnull\tnull\tnull\tnull\t-1\tnull\tnull\t2\t1\t1\tscan=1\tnull\tnull\tnull\tnull\t0.9,0.9\tPEPTIDE\n'
                'PSM\tNOVEL\t1\tnull\tnull\tnull\tnull\tnull\t-1\tnull\tnull\t2\t1\t1\tscan=2\tnull\tnull\tnull\tnull\t0.8,0.8\tNOVEL\n')
            (aa/'aa_statistics_table.csv').write_text(
                'mass shift,# peptides in bin,is reference,is isotope\n'
                '0.0,10,True,False\n'
                '+57.0214,3,False,False\n'
                '+1.0030,4,False,True\n')
            (aa/'interpretations.json').write_text(json.dumps({'+57.0214':[{'label':'Carbamidomethyl (65% match)','type':'unimod'}]}))
            (aa/'localization_statistics.csv').write_text(
                'mass shift,# peptides in bin,is isotope,localization\n'
                "+57.0214,3,False,\"{'C_+57.0214': 2, 'non-localized': 1}\"\n")

            r = run(['python3', str(SCRIPT), '--dirs', str(fp.parent), str(cn), str(aa)], cwd=root, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            summary = (root/'summary.tsv').read_text()
            self.assertIn('S1', summary)
            self.assertIn('S1\t2\t2\t1\t0.0\t0.0\t0.0\t100\t25.0\t5.0', summary)
            self.assertTrue((root/'opensearch_sequence_overlap_mqc.json').exists())
            self.assertTrue((root/'opensearch_provenance_mqc.json').exists())
            self.assertTrue((root/'opensearch_source_reports_mqc.html').exists())
            self.assertTrue((root/'provenance.tsv').exists())
            html = (root/'opensearch_aastat_annotations_mqc.html').read_text()
            self.assertIn('id: opensearch_aastat_annotations', html)
            self.assertIn('section_name: AA_stat Mass-Shift Annotations', html)

if __name__ == '__main__':
    unittest.main()
