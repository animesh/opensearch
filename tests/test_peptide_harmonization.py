#!/usr/bin/env python3
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "opensearch_summary.py"
spec = importlib.util.spec_from_file_location("osv3", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def assert_equal(a, b):
    assert a == b, (a, b)

# FragPipe integer-mass notation vs Casanovo named ProForma.
fseq, fmods = m.parse_modified_peptide("M[16]PEPTC[57]IDE")
cseq, cmods = m.parse_proforma("M[Oxidation]PEPTC[Carbamidomethyl]IDE")
assert_equal(fseq, cseq)
assert m.mods_match(fmods, cmods), (fmods, cmods)

# N-terminal modification.
fseq, fmods = m.parse_modified_peptide("[42]PEPTIDE")
cseq, cmods = m.parse_proforma("[Acetyl]-PEPTIDE")
assert_equal(fseq, cseq)
assert m.mods_match(fmods, cmods), (fmods, cmods)

# Sequence-only comparison ignores modifications.
assert_equal(m.strip_modifications("M[16]PEPTIDE"), "MPEPTIDE")
assert_equal(m.strip_modifications("M[Oxidation]PEPTIDE"), "MPEPTIDE")

# I/L-equivalent comparison.
assert_equal(m.il_key("PEPTIDEIL"), m.il_key("PEPTIDELI"))

# Different modification should remain different peptidoforms.
a = m.peptidoform_key("PEPTIDE", ((3, ("m", 16.0)),))
b = m.peptidoform_key("PEPTIDE", ((3, ("m", 15.0)),))
assert a != b

print("All OpenSearch v3 peptide harmonization tests passed.")
