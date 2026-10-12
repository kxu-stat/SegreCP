# Reproduce the 2018–2019 trade tensor analysis

This source-checkout-only directory reconstructs the paper's rank-three,
five-fold trade analysis: six 2018 exporter-level predictors are used to
predict each exporter's 2019 trade tensor. It contains download and
preprocessing code, source URLs, experiment settings, and SHA-256 references,
but no BACI, Gravity, WDI, WITS/TRAINS, or processed numerical data. The R
package build excludes the entire `reproducibility/` directory.

The frozen specification uses ordered trade-agreement scope, square-root 2018
trade weights, untransformed tariffs, and training-fold response centering.
Historical artifact identifiers remain inside the manifests and scripts to
preserve exact byte-level provenance; they are not user-facing method options.

For a GitHub release, publish the `software/SegreCP/` source tree only.
The parent research workspace contains large restricted raw and processed
data and must not be included in that repository.

## Run locally

Use a private, absolute data path outside the `SegreCP` package checkout.
The download and preprocessing steps may need several gigabytes of disk space.
From the package source directory:

```sh
Rscript reproducibility/real_data/run.R \
  --prepare --data-dir=/absolute/path/to/private-trade-work
```

For a fully downloaded private source tree, the same preprocessing can run
without network access:

```sh
Rscript reproducibility/real_data/run.R \
  --prepare --offline --data-dir=/absolute/path/to/private-trade-work
```

Inspect an existing private run with:

```sh
Rscript reproducibility/real_data/run.R \
  --check-inputs --data-dir=/absolute/path/to/private-trade-work
```

`--check-inputs` requires all 82 frozen official source files to match the
bundled source hashes and all four derived artifacts to match the frozen
byte hashes. `--check-derived` checks only the four processed files against
their frozen byte hashes. It can verify an existing processed bundle when
the large raw files are unavailable, but it does not establish that the
official source data can still be rebuilt. `--check-structure` checks source
hashes and derived dimensions,
keys, finite values, and fold balance while reporting any derived byte
differences. A structural pass alone does not establish an exact frozen replay.

Requirements: R with `data.table`, `digest`, `jsonlite`, and `xml2`.
Fresh downloads additionally require `bash`, `curl`, and `unzip`. The
real-data command-line workflow currently expects POSIX absolute paths and a
Unix-like shell (macOS, Linux, or Windows Subsystem for Linux); the installed
R package's tensor API is separate from these scripts. The
historical run used R 4.4.3 on Apple arm64; serialization and floating point
results may vary across environments. If a provider has changed a payload,
`--prepare` stops before data transformation and reports the first divergent
source paths.

## Official download paths and integrity

The complete per-file URLs and frozen SHA-256 values are in
`reference_source_manifest_sha256.csv`. The fetch scripts use:

| Source | Version or endpoint |
|---|---|
| CEPII BACI HS2017 | [`BACI_HS17_V202601.zip`](https://www.cepii.fr/DATA_DOWNLOAD/baci/data/BACI_HS17_V202601.zip) |
| CEPII Gravity | [`Gravity_rds_V202211.zip`](https://www.cepii.fr/DATA_DOWNLOAD/gravity/data/Gravity_rds_V202211.zip) |
| World Bank WDI and LPI | [World Bank API](https://api.worldbank.org/v2/) with the indicator, year, and query parameters recorded in the manifest |
| WITS / UNCTAD TRAINS | [WITS API](https://wits.worldbank.org/API/V1/wits/datasource/trn/country/ALL); exact reporter/year/type URLs are recorded in the manifest |

The CEPII archive hashes are checked before extraction. The run then checks
the frozen content of the remaining official inputs. WITS response JSON has
a generated `header.prepared` timestamp; the reference content hash
normalizes that field only. Other source files require exact byte and size
matches. The official Gravity V202211 ZIP and its extracted
`Gravity_V202211.rds` member were freshly recovered from CEPII and matched
the frozen reference. Their exact sizes and SHA-256 values are in the source
manifest. They are
two files for one source release, not two independently sourced datasets.
Availability and access through provider URLs can change; a fresh full
download of all sources, especially WITS API responses, has not been
revalidated as part of this package split. A provider may reject an automated
request. Follow
the provider's official access path if that occurs, place the same versioned
files in the private data directory, and run `--prepare --offline`.

Read the [WITS / UNCTAD TRAINS terms](https://wits.worldbank.org/data/public/UNCTAD_TermsofUseandDisclaimer.pdf)
before obtaining tariff data. They restrict redistribution of the data or
parts of it without written permission. This workflow has each user obtain
the official inputs locally and does not redistribute downloaded or derived
tariff values.

## Outputs

The private data directory receives `source_archives/`, `raw_official/`,
`metadata/`, `processed_official/`, and `analysis_official/`. The four
analysis inputs are the response tensor, six-feature exporter design, common
five-fold assignments, and preprocessing registry. Their exact relative paths
and hashes are in `reference_derived_sha256.csv`. A reader need not navigate
those paths manually:

```r
# From the SegreCP source directory, after preparing the private data.
source("reproducibility/real_data/read_prepared_trade_data.R")
trade <- read_prepared_trade_data("/absolute/path/to/private-trade-work")
dim(trade$X)   # 80 exporters x 6 features
dim(trade$Y)   # 80 x 25 x 21 x 8 response tensor
table(trade$fold)  # 16 exporters in each of five folds
```

The reader verifies the four frozen processed-file hashes, the registered
preprocessing choice, shapes, keys, feature order, and row alignment. It is
specific to the paper's prepared trade data; `SegreCP::segre_regress()` itself
accepts general user-supplied tensors and does not require these files. The
reader does **not** validate or download the 82 raw source files; use
`run.R --check-inputs` for that stronger provenance check.

The output tensor has 80 exporters, 25 importer partners, 21 HS sections,
and eight response channels. The exporter-level design has six features. Fold
standardization and response centering belong to the later fitting stage.
The historical scripts and settings are reproduced here; the `reference_derived_sha256.csv`
file permits an exact artifact comparison.

`--prepare` refuses to overwrite any of the four selected derived inputs
by default. After an interrupted run, use a fresh private directory or
inspect the existing files and explicitly pass `--overwrite-derived` to
rebuild them. Source downloads are checked before they replace any existing
source file.

The download and tensor-building code was adapted from the frozen official
2018-to-2019 analysis. Path handling now points to the private
`--data-dir`; the historical design and folds are rebuilt in a compact
script. In an isolated offline rebuild, all 82 official source files passed
their frozen content checks and all four newly generated analysis inputs
matched the historical files byte for byte and by SHA-256. This establishes
the local frozen-source-to-input transformation, but not that every provider
will still deliver the same file to a new user. Offline preparation rebuilds
the reporter-query metadata locally from the checked raw sources; no
pre-existing metadata files are required.

## Refit the authors' methods

Install `SegreCP` first, then use the four frozen processed inputs in a
private data directory:

```sh
Rscript reproducibility/real_data/fit_authors_methods.R \
  --data-dir=/absolute/path/to/private-trade-work --check-inputs

Rscript reproducibility/real_data/fit_authors_methods.R \
  --data-dir=/absolute/path/to/private-trade-work \
  --smoke --output-dir=/absolute/path/to/private-smoke-output

Rscript reproducibility/real_data/fit_authors_methods.R \
  --data-dir=/absolute/path/to/private-trade-work \
  --folds=1,2,3,4,5 --repetitions=5 \
  --output-dir=/absolute/path/to/private-full-output
```

The runner verifies the four input byte hashes, uses the common folds and
training-fold transformations, and fits rank-3 RGD and RGN through the
installed public package. `--smoke` runs one fold for one iteration and
is only a wiring check. The full command runs 5 folds x 5 timing repetitions
x 2 methods, saving per-fold accuracy, summaries, objective paths, input
hashes, and the R session. Accuracy is repeated on the same folds; timing
can vary by machine. A local 50-fit check agreed with the frozen
RGD/RGN RMSE and R-squared rows to a maximum absolute difference of
1.3e-14. Together with the byte-identical source-to-input rebuild above,
this validates the authors' two-method chain in stages. It is not a fresh
seven-method raw-source-to-table replay.

## Scope of the historical result

These scripts preserve the frozen trade-analysis sample selection and
preprocessing. The first WITS reporter selection uses **2019 tariff
availability** to pre-screen importers before the 2018 analysis; that screen
excludes KOR and MEX. Thus this exact frozen replay is not a strictly
2018-information-only prospective design. A revised 2018-only selection
would require a separately labeled sample, refit, and numerical table.

This directory rebuilds data inputs. It does not recreate the historical
seven-method table: that driver used older algorithm code outside this
public package. The SegreCP package exposes the current RGD/RGN and
initialization methods, and any new comparison table must be validated
against those public implementations before it is called a paper
reproduction.
