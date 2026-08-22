import sys
from pathlib import Path
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from count_mzml_spectra import count


def test_count_mzml():
    text = '''<?xml version="1.0"?><mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList count="3">
    <spectrum id="s1"><cvParam accession="MS:1000511" value="1"/></spectrum>
    <spectrum id="s2"><cvParam accession="MS:1000511" value="2"/></spectrum>
    <spectrum id="s3"><cvParam accession="MS:1000511" value="2"/></spectrum>
    </spectrumList></run></mzML>'''
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "x.mzML"
        p.write_text(text)
        assert count(p) == (3, 2)
