# Reproduce the trade-data analysis

This directory contains the data preparation and fitting scripts for the
2018-to-2019 trade analysis of the paper (Section 5 and Appendix D). Six exporter-level features from 2018
predict a 2019 response tensor. The scripts belong to the GitHub source
checkout, not the installed R package. Downloaded and processed data stay in
a private directory outside the checkout; neither is distributed here.

The analysis uses 80 exporters, 25 importers, 21 HS sections, eight response
channels, CP rank three, and five common folds of 16 exporters. The registered
preprocessing uses three-level trade-agreement scope, square-root 2018 trade
weights, untransformed tariffs, and training-fold feature standardization and
response centering.

## Data sources

| Input | Official source |
|:--|:--|
| Trade values and quantities | [CEPII BACI HS17 V202601](https://www.cepii.fr/DATA_DOWNLOAD/baci/data/BACI_HS17_V202601.zip) |
| Bilateral covariates | [CEPII Gravity V202211](https://www.cepii.fr/DATA_DOWNLOAD/gravity/data/Gravity_rds_V202211.zip) |
| GDP, population, and logistics | [World Bank WDI and LPI API](https://api.worldbank.org/v2/) |
| MFN tariffs | [WITS / UNCTAD TRAINS API](https://wits.worldbank.org/) |

The scripts specify the indicator, year, and reporter requests and check the
downloaded data against the versions used for the analysis. Provider access
and responses may change. If a download is rejected or differs from the
recorded version, obtain the versioned files through the provider's official
route and use the offline preparation command below. The CEPII Gravity archive
contains the `Gravity_V202211.rds` input used by the scripts.

Read the [UNCTAD TRAINS terms](https://wits.worldbank.org/data/public/UNCTAD_TermsofUseandDisclaimer.pdf)
before obtaining tariffs. They restrict redistribution without written
permission; this repository provides code and source links, not tariff values.

## Prepare the analysis inputs

Start in the `SegreCP` source checkout. You need R >= 4.2, a C++ build
toolchain, and the R packages `Rcpp`, `rTensor`, `data.table`, `digest`,
`jsonlite`, and `xml2`. Online preparation also uses `bash`, `curl`, and
`unzip`. The scripts expect a Unix-like shell and an absolute data path
outside the checkout. Allow several gigabytes of disk space.

```sh
Rscript -e 'install.packages(c("Rcpp", "rTensor", "data.table", "digest", "jsonlite", "xml2"), repos = "https://cloud.r-project.org")'
R CMD INSTALL .
TRADE_DATA_DIR=/absolute/path/to/private-trade-work
Rscript reproducibility/real_data/run.R \
  --prepare --data-dir="$TRADE_DATA_DIR"
Rscript reproducibility/real_data/run.R \
  --check-inputs --data-dir="$TRADE_DATA_DIR"
```

If the official source files are already present under the private directory
at the paths listed in the source manifest, replace `--prepare` with
`--prepare --offline`. Preparation does not overwrite existing analysis inputs
unless `--overwrite-derived` is explicitly requested; use a fresh directory
for an independent run. `--check-inputs` checks both the source files and the
four analysis inputs. If only the processed files remain, `--check-derived`
checks those inputs alone.

The four analysis inputs are the response tensor, six-feature exporter design,
fold assignments, and preprocessing registry. This helper loads and aligns
them:

```r
source("reproducibility/real_data/read_prepared_trade_data.R")
trade <- read_prepared_trade_data("/absolute/path/to/private-trade-work")
dim(trade$X)       # 80 x 6
dim(trade$Y)       # 80 x 25 x 21 x 8
table(trade$fold)  # 16 exporters per fold
```

## Fit the methods (Table 1 and Table 8)

From the repository root, with `TRADE_DATA_DIR` as above:

```sh
Rscript reproducibility/real_data/run_real_data.R \
  "$TRADE_DATA_DIR/analysis_official" reproducibility/output/real_data
Rscript reproducibility/real_data/make_tables.R reproducibility/output/real_data
```

`run_real_data.R` fits, in each of the five folds, CP-ALS, RGD and RGN from a
shared CPCA start, RRR, Tucker-RGD, Tucker-RGN and TT-IHT for 300 iterations,
and HOPLS with one to three latent blocks, and writes the held-out test RMSE,
test R² and training RMSE after 3, 5, 10, 30 and 300 iterations to
`accuracy.csv`. `make_tables.R` writes the rows of Table 1 (iteration 30) and
Table 8 (`table1_rows.tex`, `paths_rows.tex`). The run takes about 10 minutes
on one core. `reproducibility/output/real_data/` contains the `accuracy.csv`
and `session.txt` of the run reported in the paper.

Preprocessing within each fold: the six features are standardized with the
training means and standard deviations, and each response channel is centered
with its training mean. The CPCA start is computed from the training
cross-moment `X'Y/n` with the coefficient modes split into {feature, partner}
and {HS section, channel}. RGD uses the fixed step `1 / lambda_max(X'X)`, so
that `alpha * sigma_X^2 = 1` with `sigma_X^2` the largest eigenvalue of
`X'X/n`; RGD and RGN use no ridge. RRR is skipped if the `MultiwayRegression`
package is not installed.

## Sample selection

The selected sample used 2019 WITS tariff availability when screening
reporters, excluding KOR and MEX. Thus the analysis is not a strictly
2018-information-only prospective design; the screen did not use 2019 trade
outcomes.
