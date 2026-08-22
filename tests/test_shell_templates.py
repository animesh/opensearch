from pathlib import Path


MODULES = [
    Path('modules/local/fragpipe.nf'),
    Path('modules/local/aa_stat.nf'),
    Path('modules/local/casanovo.nf'),
]


def test_fragpipe_always_copies_inputs():
    text = Path('modules/local/fragpipe.nf').read_text(encoding='utf-8')
    assert "stageInMode 'copy'" in text
    assert 'stageInMode params.fragpipe_stage_mode' not in text


def test_casanovo_uses_installed_default_model_and_absolute_mzml():
    text = Path('modules/local/casanovo.nf').read_text(encoding='utf-8')
    assert 'sequence' in text
    assert ' --model ' not in text
    assert 'realpath' in text
    assert 'default model' in text.lower()


def test_aastat_uses_absolute_mzml_and_pepxml():
    text = Path('modules/local/aa_stat.nf').read_text(encoding='utf-8')
    assert r'MZML="\$(realpath' in text
    assert r'PEPXML="\$(realpath' in text
    assert '--mzml "\\$MZML"' in text
    assert '--pepxml "\\$PEPXML"' in text


def test_versions_variables_are_defined_before_heredocs():
    expected = {
        Path('modules/local/fragpipe.nf'): 'fragpipe:',
        Path('modules/local/aa_stat.nf'): 'AA_stat:',
        Path('modules/local/casanovo.nf'): 'casanovo:',
    }
    for path, marker in expected.items():
        text = path.read_text(encoding='utf-8')
        assert 'cat > versions.yml <<END_VERSIONS' in text
        assert marker in text
