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
