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
            fpdir = fp.parent
            cn = root / f'{sample}.DN'
            aa = root / f'{sample}.AA_statm'
            fp.mkdir(parents=True); cn.mkdir(); aa.mkdir()

            (fp/'peptide.tsv').write_text('Peptide\tIs Decoy\tIs Contaminant\nPEPTIDE\tfalse\tfalse\nOTHER\tfalse\tfalse\n')
            (fp/'psm.tsv').write_text(
                'Spectrum\tPeptide\tCharge\tIs Decoy\tIs Contaminant\tNumber of Missed Cleavages\n'
                'S.1.1.2\tPEPTIDE\t2\tfalse\tfalse\t0\n'
                'S.2.2.2\tOTHER\t2\tfalse\tfalse\t1\n')
            (fp/'protein.tsv').write_text('Protein ID\tIs Decoy\tIs Contaminant\nP1\tfalse\tfalse\n')
            (fpdir/'spectrum_count.tsv').write_text('sample\ttotal_spectra\tms2_spectra\nS1\t10\t8\n')
            (cn/'casanovo_test.log').write_text(
                'Sequenced 8 spectra\n'
                '50 spectra (50.00%) scored ≥ 0.00\n'
                '4 spectra (50.00%) scored ≥ 0.50\n'
                '1 spectra (12.50%) scored ≥ 0.90\n'
                '0 spectra (0.00%) scored ≥ 0.95\n'
                '0 spectra (0.00%) scored ≥ 0.99\n')
            (cn/'casanovo_test.mztab').write_text(
                'PSM\tPEPTIDE\t1\tnull\tnull\tnull\tnull\tnull\t-1\tnull\tnull\t2\t1\t1\tscan=1\tnull\tnull\tnull\tnull\t0.9,0.9\tPEPTIDE\n'
                'PSM\tNOVEL\t1\tnull\tnull\tnull\tnull\tnull\t-1\tnull\tnull\t2\t1\t1\tscan=2\tnull\tnull\tnull\tnull\t0.8,0.8\tNOVEL\n')
            (aa/'aa_statistics_table.csv').write_text(
                'mass shift,# peptides in bin,is reference,is isotope\n0.0,10,True,False\n+57.0214,3,False,False\n+1.0030,4,False,True\n')
            (aa/'interpretations.json').write_text(json.dumps({'+57.0214':[{'label':'Carbamidomethyl (65% match)','type':'unimod'}]}))
            (aa/'localization_statistics.csv').write_text('mass shift,# peptides in bin,is isotope,localization\n+57.0214,3,False,"{\'C_+57.0214\': 2, \'non-localized\': 1}"\n')

            r = run(['python3', str(SCRIPT), '--dirs', str(fpdir), str(cn), str(aa)], cwd=root, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            summary = (root/'summary.tsv').read_text()
            self.assertIn('S1\t10\t8\t2\t2\t2\t1', summary)
            data = json.loads((root/'opensearch_generalstats_mqc.json').read_text())
            self.assertNotIn('input_spectra', data['headers'])
            self.assertIn('total_spectra', data['headers'])
            self.assertIn('ms2_spectra', data['headers'])
            self.assertEqual(data['data']['S1']['spectrum_both'], 2)
            self.assertTrue((root/'opensearch_sequence_overlap_mqc.json').exists())
            self.assertTrue((root/'opensearch_provenance_mqc.json').exists())
            self.assertTrue((root/'opensearch_source_reports_mqc.html').exists())
            self.assertTrue((root/'provenance.tsv').exists())

if __name__ == '__main__':
    unittest.main()
