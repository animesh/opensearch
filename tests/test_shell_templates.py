from pathlib import Path


MODULES = [
    Path('modules/local/fragpipe.nf'),
    Path('modules/local/aa_stat.nf'),
    Path('modules/local/casanovo.nf'),
]


def test_versions_blocks_do_not_leave_nextflow_expressions_for_bash():
    for path in MODULES:
        text = path.read_text(encoding='utf-8')
        assert r'\${task.process}' not in text, f'Bash would see an unexpanded task.process in {path}'
        assert '"${task.process}":' in text, f'Nextflow task.process is not used in {path}'
        assert 'cat > versions.yml <<END_VERSIONS' in text


def test_version_shell_variables_are_defined_before_the_heredoc():
    expected = {
        Path('modules/local/fragpipe.nf'): 'FRAGPIPE_VERSION=',
        Path('modules/local/aa_stat.nf'): 'AA_STAT_VERSION=',
        Path('modules/local/casanovo.nf'): 'CASANOVO_VERSION=',
    }
    for path, marker in expected.items():
        text = path.read_text(encoding='utf-8')
        start = text.index(marker)
        end = text.index('END_VERSIONS', start)
        block = text[start:end]
        assert marker in block
        assert 'versions.yml' in block


def test_casanovo_uses_installed_default_model():
    text = Path('modules/local/casanovo.nf').read_text(encoding='utf-8')
    assert 'sequence ${mzml_file} --output_dir' in text
    assert '--model' not in text
    assert 'default model' in text.lower()
