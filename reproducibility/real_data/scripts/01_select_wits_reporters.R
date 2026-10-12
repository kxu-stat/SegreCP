#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(xml2)
})

data_dir <- Sys.getenv("SEGRECP_REAL_DATA_DIR", unset = "")
if (!nzchar(data_dir)) stop("SEGRECP_REAL_DATA_DIR must name a private data directory")
root <- normalizePath(data_dir, mustWork = TRUE)
raw_dir <- file.path(root, "raw_official")
meta_dir <- file.path(root, "metadata")
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

read_wb <- function(filename) {
  payload <- jsonlite::fromJSON(file.path(raw_dir, filename), simplifyDataFrame = TRUE)
  ans <- as.data.table(payload[[2L]])
  ans[, .(
    iso3 = countryiso3code,
    year = as.integer(date),
    value = as.numeric(value)
  )]
}

message("Reading official BACI 2018 importer totals...")
baci_2018 <- fread(
  file.path(raw_dir, "BACI_HS17_Y2018_V202601.csv"),
  select = c("j", "v"),
  colClasses = list(integer = "j", numeric = "v"),
  showProgress = TRUE
)
imports <- baci_2018[
  is.finite(v) & v > 0,
  .(import_value_2018_kusd = sum(v)),
  by = j
]
rm(baci_2018)
gc()

countries <- fread(file.path(raw_dir, "country_codes_V202601.csv"))
setnames(countries, "country_code", "j")
countries[, baci_country_code := sprintf("%03d", as.integer(j))]
imports <- merge(
  imports,
  countries[, .(j, baci_country_code, country_name, imp_iso3 = country_iso3)],
  by = "j",
  all.x = TRUE
)
setorder(imports, -import_value_2018_kusd)
imports[, rank_by_2018_imports := seq_len(.N)]

gdp <- read_wb("wdi_gdp_current_usd_2018_2019.json")
gdp_pc <- read_wb("wdi_gdp_pc_current_usd_2018_2019.json")
pop <- read_wb("wdi_population_2018_2019.json")
lpi <- read_wb("world_bank_lpi_overall_2018.json")

complete_two_years <- function(dt) {
  dt[
    year %in% 2018:2019 & is.finite(value) & !is.na(iso3) & nzchar(iso3),
    .(n_years = uniqueN(year)),
    by = iso3
  ][n_years == 2L, iso3]
}

complete_wdi <- Reduce(
  intersect,
  list(
    complete_two_years(gdp),
    complete_two_years(gdp_pc),
    complete_two_years(pop)
  )
)
complete_lpi <- lpi[year == 2018L & is.finite(value), unique(iso3)]

message("Reading official WITS/TRAINS 2019 availability...")
wits_doc <- read_xml(file.path(raw_dir, "wits_trains_dataavailability_2019.xml"))
wits_nodes <- xml_find_all(wits_doc, ".//wits:reporter", xml_ns(wits_doc))
wits_availability <- data.table(
  wits_reporter_code = xml_attr(wits_nodes, "countrycode"),
  wits_source_iso3 = xml_attr(wits_nodes, "iso3Code"),
  wits_country_name = xml_text(xml_find_first(wits_nodes, "wits:name", xml_ns(wits_doc))),
  wits_year = as.integer(xml_text(xml_find_first(wits_nodes, "wits:year", xml_ns(wits_doc)))),
  wits_nomenclature = xml_attr(
    xml_find_first(wits_nodes, "wits:reporternernomenclature", xml_ns(wits_doc)),
    "reporternernomenclaturecode"
  ),
  wits_aveestimated_available = xml_text(
    xml_find_first(
      wits_nodes,
      "wits:isspecificdutyexpressionestimatedavailable",
      xml_ns(wits_doc)
    )
  ),
  wits_last_updated = xml_text(
    xml_find_first(wits_nodes, "wits:lastupdateddate", xml_ns(wits_doc))
  )
)
wits_availability <- unique(wits_availability[wits_year == 2019L])

# In 2019 the EU-28 applied a common external tariff. WITS reports that
# schedule under EUN/918 rather than separately for every member state.
eu28_2019 <- c(
  "AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA",
  "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD",
  "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE", "GBR"
)
imports[, wits_source_iso3 := fifelse(imp_iso3 %chin% eu28_2019, "EUN", imp_iso3)]
imports <- merge(imports, wits_availability, by = "wits_source_iso3", all.x = TRUE)
imports[, has_wits_tariff_2019 :=
  !is.na(wits_reporter_code) & nzchar(wits_reporter_code) &
  wits_year == 2019L & !is.na(wits_nomenclature) & nzchar(wits_nomenclature)]

message("Reading official CEPII Gravity coverage...")
gravity <- readRDS(file.path(raw_dir, "Gravity_V202211.rds"))
gravity_coverage <- as.data.table(gravity)[
  year %in% 2018:2019 & !is.na(iso3_d),
  .(n_gravity_years = uniqueN(year)),
  by = .(imp_iso3 = iso3_d)
]
rm(gravity)
gc()

imports[, has_wdi_2018_2019 := imp_iso3 %chin% complete_wdi]
imports[, has_lpi_2018 := imp_iso3 %chin% complete_lpi]
imports <- merge(imports, gravity_coverage, by = "imp_iso3", all.x = TRUE)
imports[is.na(n_gravity_years), n_gravity_years := 0L]
imports[, has_gravity_2018_2019 := n_gravity_years == 2L]
imports[, eligible_for_query :=
  !is.na(imp_iso3) & nzchar(imp_iso3) &
  has_wdi_2018_2019 & has_lpi_2018 & has_gravity_2018_2019 &
  has_wits_tariff_2019]

eligible <- imports[eligible_for_query == TRUE]
setorder(eligible, rank_by_2018_imports)
query_set <- eligible[seq_len(min(35L, .N))]
query_set[, query_order := seq_len(.N)]

fwrite(imports, file.path(meta_dir, "baci_2018_importer_coverage.csv"))
fwrite(
  query_set[, .(
    query_order,
    baci_country_code,
    imp_iso3,
    country_name,
    wits_reporter_code,
    wits_source_iso3,
    wits_country_name,
    wits_nomenclature,
    wits_aveestimated_available,
    wits_last_updated,
    import_value_2018_kusd,
    rank_by_2018_imports
  )],
  file.path(meta_dir, "wits_reporters_2019.csv")
)

wits_unique_queries <- unique(query_set[, .(
  wits_reporter_code,
  wits_source_iso3,
  wits_country_name,
  wits_nomenclature,
  wits_aveestimated_available,
  wits_last_updated
)])
setorder(wits_unique_queries, wits_source_iso3)
wits_unique_queries[, wits_query_order := seq_len(.N)]
setcolorder(wits_unique_queries, "wits_query_order")
fwrite(
  wits_unique_queries,
  file.path(meta_dir, "wits_unique_queries_2019.csv")
)

message("Selected ", nrow(query_set), " WITS reporter candidates.")
print(query_set[, .(
  query_order,
  baci_country_code,
  imp_iso3,
  country_name,
  wits_reporter_code,
  wits_source_iso3,
  wits_nomenclature,
  import_value_2018_kusd,
  rank_by_2018_imports
)])
message("These map to ", nrow(wits_unique_queries), " unique WITS schedules.")
