<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-opensearch_logo_dark.png">
    <img alt="nf-core/opensearch" src="docs/images/nf-core-opensearch_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/opensearch)
[![GitHub Actions CI Status](https://github.com/nf-core/opensearch/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/opensearch/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/opensearch/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/opensearch/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/opensearch/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/opensearch)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23opensearch-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/opensearch)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/opensearch** is a Nextflow pipeline for proteomics open-search analysis of DDA-mode Orbitrap raw files and timsTOF `.d` directories.

The workflow runs FragPipe first, then optionally runs Casanovo (de novo sequencing) and AA_stat (modification profiling) from FragPipe outputs, finally reports to [MultiQC](https://docs.seqera.io/multiqc/getting_started/installation). Inputs can be provided using a samplesheet (`ID`, `raw-file-name`) or by scanning an input directory using filename patterns. 

### Quick start

Install [Java](https://www.oracle.com/java/technologies/downloads/#java21), [Mono](https://www.mono-project.com/download/stable/), [Fragpipe](https://github.com/Nesvilab/FragPipe/releases#release-24.0), [MultiQC](https://docs.seqera.io/multiqc/getting_started/installation) 	
`pip install multiqc`. The repository includes `fp.dl.workflow.txt` as the canonical Open-search template for timsTOF. `--fragpipe_workflow` may point to another FragPipe workflow. The pipeline renders a per-sample workflow, overriding only environment paths and input-type-dependent settings. Install the optional tools because they provide complementary evidence from the FragPipe-generated `_calibrated.mzML` output produced by FragPipe:

```bash
pip install AA_stat casanovo
whereis AA_stat casanovo
```

For a local run, use `conf/local.config` a custom configuration and provide fasta-database/tool paths, something like

```bash
 ./nextflow run . -c conf/local.config --input_dir /mnt/f/tk/MSTK/ --raw_pattern '*TK12_CTR?_*.raw' --fragpipe_workflow fp.dl.workflow.txt --fragpipe_bin /root/fragpipe24v/bin/fragpipe --fragpipe_database /root/fragpipe24v/2024-06-01-decoys-contam-UP000005640.fas --fragpipe_tools_folder /root/fragpipe24v/tools --fragpipe_diann /root/fragpipe24v/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 --fragpipe_python /usr/bin/python3 --fragpipe_crystalc false --aa_stat_bin /root/miniforge3/bin/AA_stat --casanovo_bin /root/miniforge3/bin/casanovo --cpus 12 --ram_gb 36 --max_concurrent 1 -resume

 N E X T F L O W   ~  version 26.04.6

Launching `./main.nf` [soggy_sanger] revision: 91d87110a1

executor >  local (5)
[9b/99cd48] NFC…PENSEARCH:OPENSEARCH:FRAGPIPE (200313_SIRI_TK12_CTR2_20200323222228) [100%] 3 of 3, cached: 2 ✔
[85/6ec35a] NFC…PENSEARCH:OPENSEARCH:CASANOVO (200313_SIRI_TK12_CTR2_20200323222228) [100%] 3 of 3, cached: 2 ✔
[e1/3739f7] NFC…OPENSEARCH:OPENSEARCH:AA_STAT (200313_SIRI_TK12_CTR2_20200323222228) [100%] 3 of 3, cached: 2 ✔
[ab/dfbfdb] NFCORE_OPENSEARCH:OPENSEARCH:OPENSEARCH_SUMMARY (integrated summary)     [100%] 1 of 1 ✔
[a3/cdb429] NFCORE_OPENSEARCH:OPENSEARCH:MULTIQC (multiqc)                           [100%] 1 of 1 ✔
-[nf-core/opensearch] Pipeline completed successfully-
Completed at: 24-Aug-2026 14:10:02
Duration    : 20m 13s
CPU hours   : 14.0 (71.3% cached)
Succeeded   : 5
Cached      : 6
```

The workflow runs FragPipe first and can then run Casanovo (de novo sequencing) and AA_stat (mass-shift/modification profiling) from the calibrated `mzML` and FragPipe outputs. The final reporting layer combines the three analyses into an integrated OpenSearch summary and also preserves the original tool-specific reports and MultiQC output.

### Resource controls and failure isolation

The heavy processes use shared command-line resource parameters. Defaults are **20 CPUs and 40 GB RAM per heavy-tool task**:

```bash
./nextflow run . -c conf/local.config \
    --input_dir "$PWD" \
    --raw_pattern '*.raw' \
    --fragpipe_workflow "$PWD/fp.dl.workflow.txt" \
    --fragpipe_bin /home/ash022/fragpipe/bin/fragpipe \
    --aa_stat_bin /home/ash022/.local/bin/AA_stat \
    --casanovo_bin /home/ash022/.local/bin/casanovo \
    --cpus 20 \
    --ram_gb 40
```

`--cpus` controls the Nextflow CPU allocation, FragPipe `--threads`, and AA_stat worker count. `--ram_gb` controls the Nextflow memory allocation and FragPipe `--ram`. `--max_concurrent` limits the number of simultaneous heavy per-sample tasks.

FragPipe is the primary analysis, but per-sample failures are isolated by default so one problematic run does not cancel the remaining samples. The integrated report records the exact FragPipe status and message. Casanovo and AA_stat are optional and failure-isolated: a failed AA_stat run records the exit code and reason, produces a visible warning, and does not prevent Casanovo or the integrated report from running. The same applies in the opposite direction. Spectrum counting is QC metadata only; failure to parse a calibrated mzML cannot invalidate a successful FragPipe search.

The report distinguishes unavailable total-spectrum counts from MS2-only fallback counts rather than relabeling an MS2 denominator as total spectra.

The integrated report now keeps spectrum counts and identification counts distinct. For each sample it reports the total spectra and MS2 spectra from the calibrated mzML produced by FragPipe, unique spectra with target/non-contaminant FragPipe PSMs, total FragPipe PSM rows, Casanovo sequenced spectra, Casanovo spectra with score >=0.50, and the overlap between FragPipe PSM spectra and Casanovo score >=0.50 spectra. Casanovo is a de novo sequencer rather than a database-search PSM engine, so the report deliberately labels these as Casanovo sequences/spectra rather than PSMs. The Casanovo score-threshold counts are taken from its mzTab `search_engine_score[1]` field and cross-checked against the Casanovo log.

The report has one authoritative MultiQC General Statistics table. The individual FragPipe/Casanovo/AA_stat `*_mqc.tsv` files remain published for debugging/backwards compatibility but are not passed to MultiQC, preventing the duplicated PSM/peptide/protein columns that appeared in earlier reports.

The integrated report treats the three tools as complementary views of the same MS/MS data:

- **FragPipe:** what can be identified by database searching?
- **Casanovo:** what can be sequenced de novo without relying on the protein database?
- **AA_stat:** what unexplained precursor/peptide mass shifts and modification patterns are present?
- **OpenSearch Summary:** what do these analyses collectively say about the dataset?

NOTE: Casanovo needs GPU to be efficient, but it doesnt have to be the latest and greatest, RTX2070 via WSL is enough

![RTX2070 on WSL](<images/Screenshot 2026-08-15 153254.png>)


Default workflow steps:

1. Discover raw inputs from a samplesheet or directory pattern(s)
2. Run FragPipe in headless mode using the provided workflow and manifest templates
3. Run Casanovo from the FragPipe-generated calibrated `mzML` files (optional)
4. Run AA_stat from generated calibrated `mzML` and `pepXML` files (optional)
5. Produce standard nf-core pipeline metadata and reports

### Integrated report

The report adds a cross-tool analysis layer rather than simply presenting three independent tool summaries. Sections are generated conditionally, so disabling Casanovo or AA_stat does not break the report.

It includes:

1. **Single authoritative General Statistics table**
   - total spectra and MS2 spectra from the calibrated mzML in new runs
   - FragPipe PSMs, peptides and proteins
   - FragPipe identification rate
   - Casanovo sequences and sequence yield
   - spectra identified by both tools
   - modified PSM percentage, contaminant percentage and missed-cleavage percentage
   - no duplicate FragPipe/Casanovo general-stat columns

2. **Sample QC summary and flags**
   - compact run-level QC table
   - low-identification and low-de-novo-yield flags
   - contaminant and missed-cleavage warnings
   - flags are descriptive heuristics, not hard acceptance criteria

3. **PTM-Shepherd modification landscape**
   - top modifications from `global.modsummary.tsv`
   - percentage of PSMs carrying each modification
   - comparison across samples

4. **Spectrum Identification Overview**
   - input spectra
   - FragPipe PSM spectra
   - Casanovo spectra
   - spectra identified by both tools
   - FragPipe-only and Casanovo-only spectra
   - spectra receiving neither identification
   - both/only categories are matched by scan number plus precursor charge

5. **Casanovo confidence**
   - ≥0.00, ≥0.50, ≥0.90, ≥0.95 and ≥0.99 thresholds
   - percentage of sequenced spectra at each threshold
   - fallback to the number of PSM rows in the mzTab if the log does not report the total

6. **Casanovo ↔ FragPipe sequence overlap**
   - exact overlap of unmodified peptide sequences
   - Casanovo-only candidate sequences
   - overlap percentage

   Casanovo-only sequences are deliberately described as *Casanovo-only candidates*, not automatically as novel peptides.

7. **AA_stat mass-shift landscape and annotations**
   - observed mass shifts
   - peptide counts
   - AA_stat annotations
   - reported Unimod match percentages and links where available
   - isotope shifts are distinguished from other shifts

8. **Precursor charge distribution**
   - PSM counts by precursor charge state
   - useful as an MS2 quality and acquisition-consistency indicator

9. **Missed-cleavage distribution**
   - PSM counts by number of missed cleavages
   - overall missed-cleavage percentage in the general statistics

10. **Protein-level reproducibility**
    - pairwise shared protein counts
    - union size
    - Jaccard similarity between samples
    - decoys and contaminants excluded from this comparison

11. **Automatic observations**
    - run-to-run differences
    - identification efficiency
    - de novo sequencing quality
    - prominent modification signals

12. **Data provenance and source reports**
    - every integrated metric is mapped to its originating program and source file
    - `pipeline_info/provenance.tsv` provides the same mapping in machine-readable form
    - the MultiQC report provides relative links to the published FragPipe, PTM-Shepherd, Casanovo and AA_stat source files/reports

The report is intended to answer not only *how many identifications were obtained*, but also *what each analysis contributes beyond the others*. The detailed FragPipe, Casanovo and AA_stat reports remain available as the technical appendix.

### Reporting architecture

The reporting flow is:

```text
RAW / .d
   |
   +--------------------+
   |                    |
FragPipe             calibrated mzML
   |                    |
   |             +------+------+
   |             |             |
   |          Casanovo       AA_stat
   |             |             |
   +-------------+-------------+
                 |
        OpenSearch Summary
                 |
              MultiQC
                 |
        multiqc_report.html
```

The integrated summary is implemented as a pipeline-specific MultiQC custom-content layer. This keeps the standard MultiQC modules intact while adding the OpenSearch interpretation layer.


## More on Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

Prepare one of the two supported input modes.

1. Samplesheet mode (`--input`):

```csv
ID,raw-file-name
sample_01,/data/orbitrap/sample_01.raw
sample_02,/data/timstof/sample_02.d
```

2. Directory mode (`--input_dir`):

- Point to a parent directory containing raw files/directories.
- Optionally set `--raw_pattern` (comma-separated globs), e.g. `*.raw` or `*.d`.
- The default pattern is `*.d,*.raw,*.RAW,*.mzML,*.mzml`.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/opensearch \
   -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --fragpipe_workflow <PATH_TO_FRAGPIPE_WORKFLOW> \
  --outdir <OUTDIR>
```

Directory mode:

```bash
nextflow run nf-core/opensearch \
  -profile <docker/singularity/.../institute> \
  --input_dir /path/to/raw_inputs \
  --raw_pattern '*.raw,*.d' \
  --fragpipe_workflow <PATH_TO_FRAGPIPE_WORKFLOW> \
   --outdir <OUTDIR>
```

Orbitrap example:

```bash
nextflow run nf-core/opensearch \
  --input_dir /path/to/orbitrap_raws \
  --raw_pattern '*.raw' \
  --fragpipe_workflow "$PWD/fp.dl.workflow.txt" \
  --outdir results
```

timsTOF example:

```bash
nextflow run nf-core/opensearch \
  --input_dir /path/to/timstof_runs \
  --raw_pattern '*.d' \
  --fragpipe_workflow "$PWD/fp.dl.workflow.txt" \
  --outdir results
```

Provide the FragPipe workflow directly:

- `--fragpipe_workflow /path/to/FragPipe.workflow`

### Important parameters

- `--fragpipe_workflow`: optional path to a FragPipe workflow. The default is the repository `fp.dl.workflow.txt`.
- The manifest is generated automatically for every staged raw input. No manifest template is required.
- `--fragpipe_bin`: optional override if FragPipe is not in the default location.
- `--run_casanovo` and `--run_aa_stat`: optional downstream steps. If binaries are missing, these steps are skipped with warnings.

### Default values

Pipeline parameter defaults (`nextflow.config`):

- `outdir: ./results`
- `raw_pattern: *.d,*.raw,*.RAW,*.mzML,*.mzml`
- `run_casanovo: true`
- `run_aa_stat: true`
- `cpus: 20`
- `ram_gb: 40`
- `max_concurrent: 1`
- `Casanovo model: installed default`
- `FragPipe input staging: copy`
- `fragpipe_bin: $FRAGPIPE_BIN`
- `aa_stat_bin: $AA_STAT_BIN`
- `casanovo_bin: $CASANOVO_BIN`
- `fragpipe_database: $FRAGPIPE_DATABASE`
- `fragpipe_tools_folder: $FRAGPIPE_TOOLS_FOLDER`
- `fragpipe_diann: $FRAGPIPE_DIANN`
- Tool output directory names are derived from the executable basenames.

Default process resources (`conf/base.config`):

- `process_low`: `2 CPUs`, `8 GB`
- `process_medium`: `cpus`, `ram_gb`
- `process_high`: `cpus`, `ram_gb`

Default publish behavior (`conf/modules.config`):

- Tool outputs are copied to directories named from the executable paths, e.g. `results/fragpipe`, `results/casanovo`, and `results/AA_stat`.
- The integrated MultiQC report is copied to `results/multiqc/`.
- `versions.yml` is not copied to output process folders

### FragPipe workflow and manifest

The repository includes `fp.dl.workflow.txt`; a custom workflow can be supplied with `--fragpipe_workflow`. For each sample, OpenSearch copies the complete raw file/directory into the pipeline launch directory and generates a four-column, tab-separated FragPipe manifest from that known absolute path: LC-MS path, Experiment, Bioreplicate, and data type. The Bioreplicate field is blank because each FragPipe task processes one file as one experiment.

The generated workflow preserves the supplied FragPipe search settings and overrides only environment/input-dependent settings. The calibrated mzML path is deterministic: `<launch-directory>/<raw-stem>_calibrated.mzML`; no filesystem search is used. For timsTOF `.d` input it sets IM-MS mode and disables Crystal-C because Crystal-C currently does not support `.d`; for regular MS it sets Regular-MS mode. It also enables calibrated mzML writing so downstream tools have a standard spectrum file when calibration succeeds:

```text
database.db-path=<--fragpipe_database>
fragpipe-config.tools-folder=<--fragpipe_tools_folder>
fragpipe-config.bin-diann=<--fragpipe_diann>
fragpipe-config.bin-python=<--fragpipe_python>
crystalc.run-crystalc=<--fragpipe_crystalc>
```

### Monitoring progress

For live pipeline status, follow the Nextflow log:

```bash
tail -f .nextflow.log
```

For sample-level FragPipe progress, inspect the active task work directory and follow the process output:

```bash
find work -maxdepth 3 -type f -name .command.out
tail -f work/<hash>/<hash>/.command.out
```

If you want both stdout and stderr together, use:

```bash
tail -f work/<hash>/<hash>/.command.log
```

Useful optional Nextflow reports:

```bash
nextflow run nf-core/opensearch \
  --input_dir /path/to/raw_inputs \
  --fragpipe_workflow "$PWD/fp.dl.workflow.txt" \
  -with-report \
  -with-trace \
  -with-timeline \
  -with-dag flowchart.png
```

### Runtime notes

- Each raw `.raw` file or `.d` directory is copied into the pipeline launch directory before FragPipe is run.
- Casanovo receives the known absolute `<raw-stem>_calibrated.mzML` path directly as a value. AA_stat receives the same calibrated mzML plus the FragPipe pepXML.
- If FragPipe fails for a sample, the pipeline retains any usable artifacts and records the real exit code in that sample's `status.tsv`; other samples continue.
- For detailed process debugging, inspect `.command.sh`, `.command.out`, `.command.err`, and `.command.log` inside the relevant `work/` directory.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/opensearch/usage) and the [parameter documentation](https://nf-co.re/opensearch/parameters).

## Pipeline output

By default, outputs are written under `--outdir` (default: `results`) in process-specific subfolders:

- `results/<FragPipe executable name>/`
- `results/<Casanovo executable name>/` (if enabled and available)
- `results/<AA_stat executable name>/` (if enabled and available)
- `results/pipeline_info/summary.tsv` and `results/pipeline_info/provenance.tsv` (integrated machine-readable summary and provenance)
- `results/multiqc/multiqc_report.html` (integrated MultiQC report)

For FragPipe, each sample is published as a work directory named like:

- `results/<FragPipe executable name>/<sample>.<FragPipe executable name>/`

For example, PTM-Shepherd summary tables are typically found at:

- `results/<FragPipe executable name>/<sample>.<FragPipe executable name>/ptm-shepherd-output/global.modsummary.tsv`

Note that FragPipe creates many nested files. Each published FragPipe directory now also contains `spectrum_count.tsv`, recording total spectra and MS2 spectra in the calibrated mzML used by the downstream tools. Older results without this file fall back to Casanovo's sequenced + skipped spectrum counts when Casanovo was run.

For more details, please refer to the [output documentation](https://nf-co.re/opensearch/output).

## Development

```bash
python -m venv nf-core
source nf-core/bin/activate
pip install nf-core
nf-core pipelines lint .
nf-core modules lint .
```


## Credits

nf-core/opensearch is maintained by the nf-core community.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#opensearch` channel](https://nfcore.slack.com/channels/opensearch) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

[opensearch](https://github.com/animesh/opensearch) is created with great help from [github-copilot](https://github.com/copilot)

Please don't forget to cite what `opensearch` is really based upon, [Fragpipe](https://github.com/Nesvilab/FragPipe), [AA_stat](https://github.com/SimpleNumber/aa_stat), and  [Casanovo](https://github.com/Noble-Lab/casanovo)!

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

