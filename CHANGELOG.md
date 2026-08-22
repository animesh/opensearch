## [Unreleased]

### Fixed

- Corrected timsTOF `.d` processing to use IM-MS mode, disable unsupported Crystal-C, and enable calibrated mzML generation.
- Made FragPipe workflow rendering explicit and testable instead of relying on inline `sed` edits.
- Passed documented FragPipe headless configuration options for tools, DIA-NN and Python.
- Isolated per-sample FragPipe failures so one failed post-processing run does not cancel the remaining batch; partial MSFragger artifacts are retained.
- Staged FragPipe mzML/pepXML outputs correctly into Casanovo and AA_stat.
- Selected Casanovo's `timstof` model automatically for `.d` input and `orbitrap` for regular MS.
- Added explicit FragPipe/Casanovo/AA_stat status fields to the integrated summary.

### Added

- `bin/render_fragpipe_workflow.py` and unit coverage for input-type-aware workflow rendering.

### Added

- Shared `--cpus` and `--ram_gb` resource controls, defaulting to 20 CPUs and 40 GB RAM for heavy tools.
- Failure-isolated Casanovo and AA_stat execution with per-sample status files and visible report diagnostics.
- Standalone mzML spectrum counter with unit tests; spectrum-counting failures no longer fail FragPipe.
- Harmonized Casanovo/FragPipe peptide comparison using scan+charge, unmodified sequence identity, I/L-equivalent identity, and position-aware modification masses.
- Separate all-Casanovo and Casanovo-score>=0.50 sequence overlap metrics.

### Added

- Provenance-aware integrated MultiQC reporting for FragPipe, PTM-Shepherd, Casanovo and AA_stat.
- Direct source-file/report links and machine-readable `provenance.tsv`.
- Local `run_local.sh` wrapper and `conf/local.config` to avoid first-run `-resume` and no-custom-config warnings.
- Cleaner MultiQC custom-content ordering and disabled plot export to avoid non-fatal colour conversion warnings.

### Fixed

- AA_stat HTML custom-content metadata format.
- AA_stat localized-bin terminology.
- MultiQC ordering configuration that referred to custom sections as modules.

# nf-core/opensearch: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/opensearch, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

### v1.0.0dev integrated-report update

#### `Added`

- Integrated FragPipe/Casanovo/AA_stat MultiQC summary.
- Casanovo confidence-threshold plot from the Casanovo log.
- Casanovo-to-FragPipe exact peptide-sequence overlap analysis.
- AA_stat non-isotope mass-shift landscape with reported interpretation annotations.
- Machine-readable `summary.tsv` in `results/pipeline_info`.
- Unit test for the integrated summary parser.

#### `Fixed`

- Casanovo MultiQC parsing no longer expects the obsolete `dataset contains` log message.
- Casanovo report now uses the `Sequenced` and confidence-threshold values actually emitted by Casanovo.
