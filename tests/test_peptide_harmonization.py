import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))

from peptide_harmonization import mods_equal, normalize_sequence, parse_casanovo_mods, parse_fragpipe_mods


def test_sequence_normalization():
    assert normalize_sequence("M[147]PEP") == "MPEP"
    assert normalize_sequence("M[147]PEPI", True) == "MPEPL"


def test_fragpipe_and_casanovo_common_modifications_match():
    fp = parse_fragpipe_mods("MPEP", "1M(15.9949)")
    cas = parse_casanovo_mods("1-Oxidation (M):UNIMOD:35", "M[Oxidation]PEP", "MPEP")
    assert mods_equal(fp, cas, tolerance=0.05)


def test_nterm_acetyl_matches():
    fp = parse_fragpipe_mods("PEPTIDE", "N-term(42.0106)")
    cas = parse_casanovo_mods("0-Acetyl (N-term):UNIMOD:1", "[Acetyl]-PEPTIDE", "PEPTIDE")
    assert mods_equal(fp, cas, tolerance=0.05)


def test_unknown_modification_is_not_guessed():
    cas = parse_casanovo_mods("3-Something (K):UNIMOD:999999", "PEK[Something]T", "PEKT")
    assert isinstance(cas[0][1], str)
