from pathlib import Path

from render_fragpipe_workflow import render


def test_render_tims_and_environment_overrides(tmp_path):
    source = tmp_path / "Open.workflow"
    dest = tmp_path / "generated.workflow"
    source.write_text(
        "# Workflow: Open\n"
        "database.db-path=/old/db.fas\n"
        "fragpipe-config.tools-folder=/old/tools\n"
        "fragpipe-config.bin-diann=/old/diann\n"
        "fragpipe-config.bin-python=/old/python\n"
        "crystalc.run-crystalc=true\n"
        "workflow.input.data-type.im-ms=false\n"
        "workflow.input.data-type.regular-ms=true\n"
        "msfragger.write_calibrated_mzml=false\n"
        "msfragger.run-msfragger=true\n"
        "database.decoy-tag=rev_\n"
    )
    import argparse

    render(argparse.Namespace(
        source=str(source), destination=str(dest), database="/new/db.fas",
        tools_folder="/new/tools", diann="/new/diann", python="/new/python",
        im_ms="true", crystalc="auto", write_calibrated=True,
    ))
    text = dest.read_text()
    assert text.count("database.db-path=") == 1
    assert "database.db-path=/new/db.fas" in text
    assert "workflow.input.data-type.im-ms=true" in text
    assert "workflow.input.data-type.regular-ms=false" in text
    assert "msfragger.write_calibrated_mzml=true" in text
    assert "crystalc.run-crystalc=false" in text
    assert "fragpipe-config.tools-folder=/new/tools" in text
    assert "database.decoy-tag=rev_" in text
