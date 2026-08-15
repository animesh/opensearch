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

Install [Fragpipe](https://github.com/Nesvilab/FragPipe/releases#release-24.0), [MultiQC](https://docs.seqera.io/multiqc/getting_started/installation) 	
`pip install multiqc` and configure [fp.dl.workflow.txt](fp.dl.workflow.txt) specifically folliowing variables, note that `database.db-path` above points to a fasta file containing decoys and contaminants generated with Fragpipe itself and tools are downloaded via its GUI as well

```bash
database.db-path=/root/fragpipe24v/2024-06-01-decoys-contam-UP000005640.fas

fragpipe-config.tools-folder=/root/fragpipe24v/tools
fragpipe-config.bin-diann=/root/fragpipe24v/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8
fragpipe-config.bin-python=/usr/bin/python3
```

probably better to install optional tools as they can take in the `_calibrated.mzML` output from Fragpipe give a second opinion on mods like [AA_stat](https://pypi.org/project/AA-stat/) `pip install AA_stat` and [Casanovo](https://pypi.org/project/casanovo/) `pip install casanovo` which can run de-novo sequecing on the data, but need find there binaries `whereis AA_stat casanovo`to analyze `--input_dir`

```bash
nextflow run . --input_dir $PWD --scripts_dir $PWD --fragpipe_bin /root/fragpipe24v/bin/fragpipe --aa_stat_bin /root/miniforge3/bin/AA_stat --casanovo_bin /root/miniforge3/bin/casanovo -resume 

 N E X T F L O W   ~  version 26.04.6

Launching `./main.nf` [scruffy_ride] revision: e70a121c00

WARN: [nf-core/opensearch] You are attempting to run the pipeline without any custom configuration!

This will be dependent on your local compute environment but can be achieved via one or more of the following:
   (1) Using an existing pipeline profile e.g. `-profile docker` or `-profile singularity`
   (2) Using an existing nf-core/configs for your Institution e.g. `-profile crick` or `-profile uppmax`
   (3) Using your own local custom config e.g. `-c /path/to/your/custom.config`

Please refer to the quick start section and usage docs for the pipeline.
 
executor >  local (4)
[cd/a3a19d] NFCORE_OPENSEARCH:OPENSEARCH:FRAGPIPE (191107_SIRI_1_TK9_ctr1) [100%] 1 of 1 ✔
[65/26abc1] NFCORE_OPENSEARCH:OPENSEARCH:CASANOVO (191107_SIRI_1_TK9_ctr1) [100%] 1 of 1 ✔
[7c/627fdc] NFCORE_OPENSEARCH:OPENSEARCH:AA_STAT (191107_SIRI_1_TK9_ctr1)  [100%] 1 of 1 ✔
[a5/f048cd] NFCORE_OPENSEARCH:OPENSEARCH:MULTIQC (multiqc)                 [100%] 1 of 1 ✔
-[nf-core/opensearch] Pipeline completed successfully-
Completed at: 15-Aug-2026 18:47:20
Duration    : 26m 43s
CPU hours   : 2.8
Succeeded   : 4
```

NOTE: Casanovo needs GPU to be efficient, but it doesnt have to be the latest and greatest, RTX2070 via WSL is enough

![RTX2070 on WSL](<images/Screenshot 2026-08-15 153254.png>)


Default workflow steps:

1. Discover raw inputs from a samplesheet or directory pattern(s)
2. Run FragPipe in headless mode using the provided workflow and manifest templates
3. Run Casanovo from generated calibrated `mzML` files (optional)
4. Run AA_stat from generated calibrated `mzML` and `pepXML` files (optional)
5. Produce standard nf-core pipeline metadata and reports


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
  --scripts_dir <DIR_WITH_fp.manifest.txt_AND_fp.dl.workflow.txt> \
  --outdir <OUTDIR>
```

Directory mode:

```bash
nextflow run nf-core/opensearch \
  -profile <docker/singularity/.../institute> \
  --input_dir /path/to/raw_inputs \
  --raw_pattern '*.raw,*.d' \
  --scripts_dir <DIR_WITH_fp.manifest.txt_AND_fp.dl.workflow.txt> \
   --outdir <OUTDIR>
```

Orbitrap example:

```bash
nextflow run nf-core/opensearch \
  --input_dir /path/to/orbitrap_raws \
  --raw_pattern '*.raw' \
  --scripts_dir /path/to/scripts \
  --outdir results
```

timsTOF example:

```bash
nextflow run nf-core/opensearch \
  --input_dir /path/to/timstof_runs \
  --raw_pattern '*.d' \
  --scripts_dir /path/to/scripts \
  --outdir results
```

If you prefer explicit file paths instead of `--scripts_dir`, use:

- `--fragpipe_manifest /path/to/fp.manifest.txt`
- `--fragpipe_workflow /path/to/fp.dl.workflow.txt`

### Important parameters

- `--scripts_dir`: recommended way to provide FragPipe templates. Must contain `fp.manifest.txt` and `fp.dl.workflow.txt`.
- `--fragpipe_manifest` and `--fragpipe_workflow`: explicit alternative to `--scripts_dir`.
- `--fragpipe_bin`: optional override if FragPipe is not in the default location.
- `--run_casanovo` and `--run_aa_stat`: optional downstream steps. If binaries are missing, these steps are skipped with warnings.

### Default values

Pipeline parameter defaults (`nextflow.config`):

- `outdir: ./results`
- `raw_pattern: *.d,*.raw,*.RAW,*.mzML,*.mzml`
- `run_casanovo: true`
- `run_aa_stat: true`
- `fragpipe_threads: 8`
- `fragpipe_ram_gb: 32`
- `fragpipe_bin: $FRAGPIPE_BIN or $HOME/fragpipe/bin/fragpipe`
- `aa_stat_bin: $AA_STAT_BIN or $HOME/.local/bin/AA_stat`
- `casanovo_bin: $CASANOVO_BIN or $HOME/.local/bin/casanovo`

Default process resources (`conf/base.config`):

- `process_low`: `2 CPUs`, `8 GB`
- `process_medium`: `4 CPUs`, `32 GB`
- `process_high`: `8 CPUs`, `32 GB`

Default publish behavior (`conf/modules.config`):

- Outputs are copied to `results/<process_name>/` (lowercase process name)
- `versions.yml` is not copied to output process folders

FragPipe template defaults in the repository (`fp.dl.workflow.txt`):

```text
database.db-path=/root/fragpipe/2024-06-01-decoys-contam-UP000005640.fas
fragpipe-config.tools-folder=/root/fragpipe/tools
fragpipe-config.bin-diann=/root/fragpipe/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8
fragpipe-config.bin-python=/usr/bin/python3
```

FragPipe manifest template default (`fp.manifest.txt`):

```text
RAWDIR  RAWFILE  DDA
```

These FragPipe template values are not automatically discovered from your system. Update them for your environment, then pass them via `--scripts_dir` (or explicit `--fragpipe_manifest` + `--fragpipe_workflow`).

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
  --scripts_dir /path/to/scripts \
  -with-report \
  -with-trace \
  -with-timeline \
  -with-dag flowchart.png
```

### Runtime notes

- FragPipe inputs are staged into the task work directory before execution. This is intentional and helps avoid reader issues with mounted source paths.
- Casanovo and AA_stat start only after FragPipe produces the expected `mzML` and `pepXML` outputs.
- If FragPipe fails, the pipeline now stops with a real process error instead of reporting a misleading "successful with errors" summary.
- For detailed process debugging, inspect `.command.sh`, `.command.out`, `.command.err`, and `.command.log` inside the relevant `work/` directory.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/opensearch/usage) and the [parameter documentation](https://nf-co.re/opensearch/parameters).

## Pipeline output

By default, outputs are written under `--outdir` (default: `results`) in process-specific subfolders:

- `results/fragpipe/`
- `results/casanovo/` (if enabled and available)
- `results/aa_stat/` (if enabled and available)

For FragPipe, each sample is published as a work directory named like:

- `results/fragpipe/<sample>.FPv24hum/`

For example, PTM-Shepherd summary tables are typically found at:

- `results/fragpipe/20250909_CSF_13_b_Slot1-32_1_11095.FPv24hum/ptm-shepherd-output/global.modsummary.tsv`

Note that FragPipe creates many nested files; calibrated mzML and pepXML outputs are also produced and then used by optional Casanovo and AA_stat steps.

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

