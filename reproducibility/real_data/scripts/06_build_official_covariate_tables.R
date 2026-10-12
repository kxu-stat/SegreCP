#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

data_dir <- Sys.getenv("SEGRECP_REAL_DATA_DIR", unset = "")
if (!nzchar(data_dir)) stop("SEGRECP_REAL_DATA_DIR must name a private data directory")
root <- normalizePath(data_dir, mustWork = TRUE)
raw_dir <- file.path(root, "raw_official")
processed_dir <- file.path(root, "processed_official")
meta_dir <- file.path(root, "metadata")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

read_wb <- function(filename, value_name) {
  payload <- fromJSON(file.path(raw_dir, filename), simplifyDataFrame = TRUE)
  ans <- as.data.table(payload[[2L]])
  ans <- ans[
    as.integer(date) == 2018L & !is.na(countryiso3code) & nzchar(countryiso3code),
    .(iso3 = countryiso3code, value = as.numeric(value))
  ]
  setnames(ans, "value", value_name)
  unique(ans, by = "iso3")
}

message("Building official 2018 WDI macro table...")
macro <- Reduce(
  function(x, y) merge(x, y, by = "iso3", all = TRUE),
  list(
    read_wb("wdi_gdp_current_usd_2018_2019.json", "gdp_usd"),
    read_wb("wdi_gdp_pc_current_usd_2018_2019.json", "gdp_pc_usd"),
    read_wb("wdi_population_2018_2019.json", "pop")
  )
)
macro[, year := 2018L]
setcolorder(macro, c("year", "iso3", "gdp_usd", "gdp_pc_usd", "pop"))
fwrite(macro, file.path(processed_dir, "wdi_macro_2018.csv"))

message("Building official 2018 LPI logistics-friction table...")
lpi <- read_wb("world_bank_lpi_overall_2018.json", "lpi_score")
lpi[, `:=`(year = 2018L, logistics_friction_proxy = -lpi_score)]
setcolorder(lpi, c("year", "iso3", "lpi_score", "logistics_friction_proxy"))
fwrite(lpi, file.path(processed_dir, "lpi_logistics_friction_2018.csv"))

message("Building official 2018 CEPII Gravity table...")
gravity <- as.data.table(readRDS(file.path(raw_dir, "Gravity_V202211.rds")))
needed <- c(
  "year", "iso3_o", "iso3_d", "distw_harmonic", "contig",
  "comlang_off", "comcol", "rta_type"
)
if (length(setdiff(needed, names(gravity)))) {
  stop("CEPII Gravity is missing required columns")
}
gravity_2018 <- gravity[year == 2018L, ..needed]
setnames(gravity_2018, "distw_harmonic", "distw_h")
gravity_2018 <- unique(gravity_2018, by = c("iso3_o", "iso3_d"))
fwrite(gravity_2018, file.path(processed_dir, "cepii_gravity_2018.csv"))
rm(gravity)
gc()

message("Expanding official 2018 tariffs from 19 schedules to 35 importers...")
tariffs <- fread(
  file.path(processed_dir, "wits_mfn_tariffs_preferred_2018.csv"),
  colClasses = list(character = c("wits_source_iso3", "hs6"))
)
importers <- fread(
  file.path(meta_dir, "wits_reporters_2018.csv"),
  colClasses = list(character = c("imp_iso3", "wits_source_iso3"))
)
importers <- unique(importers[, .(
  imp_iso3,
  country_name,
  wits_source_iso3,
  wits_reporter_code,
  wits_nomenclature,
  wits_aveestimated_available
)])
tariffs_importer <- merge(
  importers,
  tariffs[, .(
    wits_source_iso3,
    hs6,
    tariff_avg_pct = tariff_simple_average_pct,
    nomenclature,
    source_datatype_request,
    source_file
  )],
  by = "wits_source_iso3",
  allow.cartesian = TRUE
)

hs_sections <- data.table(
  hs_section = sprintf("HS%02d", 1:21),
  hs2_lo = c(1, 6, 15, 16, 25, 28, 39, 41, 44, 47, 50, 64, 68, 71, 72, 84, 86, 90, 93, 94, 97),
  hs2_hi = c(5, 14, 15, 24, 27, 38, 40, 43, 46, 49, 63, 67, 70, 71, 83, 85, 89, 92, 93, 96, 97)
)
tariffs_importer[, hs2 := as.integer(substr(hs6, 1L, 2L))]
tariffs_importer[, hs_section := NA_character_]
for (row in seq_len(nrow(hs_sections))) {
  tariffs_importer[
    hs2 >= hs_sections$hs2_lo[row] & hs2 <= hs_sections$hs2_hi[row],
    hs_section := hs_sections$hs_section[row]
  ]
}
setorder(tariffs_importer, imp_iso3, hs6)
fwrite(
  tariffs_importer,
  file.path(processed_dir, "wits_mfn_tariffs_by_importer_2018.csv")
)

covariate_qa <- data.table(
  dataset = c("WDI macro", "LPI", "CEPII Gravity", "WITS tariff importer-HS6"),
  rows = c(nrow(macro), nrow(lpi), nrow(gravity_2018), nrow(tariffs_importer)),
  duplicate_keys = c(
    anyDuplicated(macro$iso3),
    anyDuplicated(lpi$iso3),
    gravity_2018[, anyDuplicated(paste(iso3_o, iso3_d))],
    tariffs_importer[, anyDuplicated(paste(imp_iso3, hs6))]
  ),
  primary_measure_missing = c(
    sum(!is.finite(macro$gdp_usd) | !is.finite(macro$gdp_pc_usd) | !is.finite(macro$pop)),
    sum(!is.finite(lpi$lpi_score)),
    sum(!is.finite(gravity_2018$distw_h)),
    sum(!is.finite(tariffs_importer$tariff_avg_pct))
  )
)
fwrite(covariate_qa, file.path(meta_dir, "official_covariate_qa_2018.csv"))

if (covariate_qa[, any(duplicate_keys != 0L)]) stop("Duplicate covariate keys remain")
if (tariffs_importer[, anyNA(hs_section)]) {
  stop("At least one WITS product could not be assigned to an HS section")
}
if (tariffs_importer[, uniqueN(imp_iso3)] != 35L) {
  stop("Tariff expansion did not cover all 35 selected importers")
}

message("Official covariate tables built successfully.")
print(covariate_qa)
