from pathlib import Path


def test_fragpipe_copies_input_to_launch_directory_and_uses_absolute_manifest():
    text = Path('modules/local/fragpipe.nf').read_text(encoding='utf-8')
    assert 'cp -a' in text
    assert 'workflow.launchDir' in text
    assert 'INPUT_COPY=' in text
    assert 'GENERATED_MANIFEST=' in text
    assert 'printf \'%s\\\\t%s' in text


def test_casanovo_uses_installed_default_model_and_absolute_mzml():
    text = Path('modules/local/casanovo.nf').read_text(encoding='utf-8')
    assert 'sequence' in text
    assert '--model' not in text
    assert 'MZML="${mzml_path}"' in text
    assert 'default model' in text.lower()
    # Casanovo interprets --output_root relative to --output_dir; passing
    # OUTDIR/sample would make it create OUTDIR/OUTDIR/sample.log.
    assert 'OUTPUT_ROOT="${sample_id}"' in text
    assert r'LOG="\$OUTDIR/\${OUTPUT_ROOT}.log"' in text
    assert r'MZTAB="\$OUTDIR/\${OUTPUT_ROOT}.mztab"' in text


def test_aastat_uses_absolute_mzml_and_pepxml():
    text = Path('modules/local/aa_stat.nf').read_text(encoding='utf-8')
    assert 'MZML="${mzml_path}"' in text
    assert r'PEPXML="\$(realpath' in text
    assert '--mzml "\\$MZML"' in text
    assert '--pepxml "\\$PEPXML"' in text


def test_versions_variables_are_defined_before_heredocs():
    expected = [
        Path('modules/local/fragpipe.nf'),
        Path('modules/local/aa_stat.nf'),
        Path('modules/local/casanovo.nf'),
    ]
    for path in expected:
        text = path.read_text(encoding='utf-8')
        assert 'cat > versions.yml <<END_VERSIONS' in text
        assert 'tool_name' in text
        # task.process must be evaluated by Nextflow before the shell runs.
        # Escaping it as \$\{task.process\} leaves a Bash bad-substitution error.
        assert '"${task.process}":' in text
        assert '"\\${task.process}\":' not in text
