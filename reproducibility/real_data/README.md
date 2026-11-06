# Trade data analysis (Section 5 and Appendix D)

A tensor-on-tensor regression of 2019 export outcomes on 2018 covariates for
80 exporters: the response of exporter m is a 25 x 21 x 8 tensor (partner x
HS section x channel), and its design is a vector of 6 features.

The prepared data used in the paper (`response_Y.csv`, `design_X.csv`,
`folds.csv`) are in the supplementary material of the paper. With them, only
the last two steps below are needed. Steps 1 to 3 rebuild these files from
the public sources.

## Steps

Run from the repository root. `<data_dir>` is any directory outside the
repository; the downloads need about 2 GB.

| Step | Command | Time |
|---|---|---|
| 1 | `Rscript reproducibility/real_data/1_download_sources.R <data_dir>` | download |
| 2 | `Rscript reproducibility/real_data/2_download_tariffs.R <data_dir>` | download |
| 3 | `Rscript reproducibility/real_data/3_build_data.R <data_dir>` | 1 min |
| 4 | `Rscript reproducibility/real_data/run_real_data.R <data_dir>/prepared <out_dir>` | 10 min |
| 5 | `Rscript reproducibility/real_data/make_tables.R <out_dir>` | seconds |

Steps 1 to 3 use the R packages `data.table`, `jsonlite` and `xml2`. Step 4
writes `accuracy.csv` (test RMSE, test R² and training RMSE of every fold,
method and checkpoint) and step 5 the rows of Table 1 and of the fitting-path
table of Appendix D.

## Sources

| Data | Source | Version | License |
|---|---|---|---|
| Trade values and quantities by HS6 product | [CEPII BACI](https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37), HS 2017 | V202601 | Etalab 2.0 |
| Distance, trade agreements | [CEPII Gravity](https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=8) | V202211 | Etalab 2.0 |
| GDP, GDP per capita, population | [World Bank WDI](https://databank.worldbank.org/source/world-development-indicators) API | last updated 2026-07-13 | CC BY 4.0 |
| Logistics Performance Index | [World Bank LPI](https://lpi.worldbank.org/) API | 2018 | CC BY 4.0 |
| MFN applied tariffs by HS6 product | [WITS](https://wits.worldbank.org/), UNCTAD TRAINS | 2018 schedules | [TRAINS terms](https://wits.worldbank.org/data/public/UNCTAD_TermsofUseandDisclaimer.pdf) |

The TRAINS terms do not allow the tariff data to be redistributed, so the
tariff files are not included anywhere; step 2 downloads them.

## What the steps do

**Step 1** downloads BACI (2018 and 2019), Gravity, the World Bank indicators
and the WITS lists of available tariff schedules.

**Step 2** chooses 35 candidate partners: the largest importers of 2018 that
have WDI data for 2018 and 2019, the 2018 LPI, Gravity records for 2018 and
2019, and a WITS tariff schedule for 2019. (Korea and Mexico have no 2019
schedule in WITS and are therefore not partners; they remain exporters.) It
then downloads the 2018 schedule of each candidate. EU members share the
schedule of the EU. Specific duties are converted to ad valorem equivalents
where WITS provides them.

**Step 3** builds the data. HS6 products are grouped into the 21 HS sections,
and trade values and quantities are summed within each exporter, importer
and section.

- Partners: the 25 largest candidates by 2018 imports that have GDP,
  population, the LPI and a tariff in every section.
- Exporters: the 80 largest exporters by 2018 exports to the partners among
  those that trade with at least 13 partners in at least 11 sections, have
  GDP and population, and have a Gravity distance to every partner.
- Tariff of a partner in a section: the mean of its HS6 tariffs weighted by
  its 2018 imports of each product; the unweighted mean when the schedule is
  not in HS 2017 (Singapore) or when the weighted products cover less than
  half of the section's imports.
- Features of an exporter, partner and section cell (2018): tariff;
  `ntb_proxy`, minus the partner's LPI score (logistics friction); log GDP of
  the partner; log distance; `rta_scope_count`, 0 without a trade agreement,
  1 for an agreement on goods only and 2 for goods and services; and
  `world_demand`, log(1 + imports of the partner in the section).
- Design `design_X.csv`: for each exporter, the average of each feature over
  its 525 cells, weighted by the square root of the cell's 2018 trade value.
- Response `response_Y.csv` (2019), for trade value v and quantity q of the
  cell: log(1 + v); log(1 + q); log(v / q), or 0 if v or q is 0; v divided by
  the exporter's exports to all countries; v divided by the partner's imports
  in the section from all countries; the Balassa index of the exporter in the
  section (its share of the section in its exports divided by the section's
  share in world exports); 1 if v > 0, else 0; and log v, or 0 if v = 0.
- Folds `folds.csv`: five folds of 16 exporters, balanced over the quintiles
  of log(1 + 2018 exports to the partners), with a fixed seed.

The sample, the features and the weights use only 2018 data; the 2019 data
enter only the response. The screening in step 2 also uses the availability
of 2019 WDI data, Gravity records and WITS schedules, but no 2019 trade.

## Source versions

BACI and Gravity are versioned releases; when V202601 is no longer the current
BACI release, CEPII keeps it on its
[archive page](https://www.cepii.fr/DATA_DOWNLOAD/baci/doc/baci_webpage.html).
The World Bank revises the WDI over time. The WDI and LPI files used for the
paper are in the supplementary material; copied into `<data_dir>/sources`
before step 1, they are used instead of a new download. Values built on
different machines can differ in the last digits because of floating-point
rounding.
