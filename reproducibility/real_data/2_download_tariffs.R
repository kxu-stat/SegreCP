# Step 2 of the trade data preparation: choose the candidate partner
# (importer) countries and download their 2018 tariff schedules.
#
#   Rscript reproducibility/real_data/2_download_tariffs.R <data_dir>
#
# Candidates: the 35 largest importers in 2018 (BACI) that have GDP, GDP per
# capita and population in 2018 and 2019 (WDI), the 2018 LPI, Gravity records
# as destination in 2018 and 2019, and a 2019 WITS tariff schedule. Step 3
# chooses the 25 partners among them.
#
# Writes <data_dir>/sources/tariff_schedules.csv (one row per candidate) and
# one WITS file per schedule in <data_dir>/sources/wits_tariffs/.

library(data.table)

data_dir <- commandArgs(trailingOnly = TRUE)[1]
source_dir <- file.path(data_dir, "sources")
tariff_dir <- file.path(source_dir, "wits_tariffs")
dir.create(tariff_dir, showWarnings = FALSE)

# Countries (ISO3 codes) with a finite value of a World Bank indicator in
# every year in `years`.
complete_countries <- function(file, years) {
  json <- jsonlite::fromJSON(file.path(source_dir, file))
  d <- data.table(iso3 = json[[2]]$countryiso3code,
                  year = as.integer(json[[2]]$date),
                  value = as.numeric(json[[2]]$value))
  d <- d[year %in% years & is.finite(value) & !is.na(iso3) & iso3 != ""]
  d[, .(n_years = uniqueN(year)), by = iso3][n_years == length(years), iso3]
}

# Tariff schedules available in WITS for one year: reporter code, ISO3 code,
# nomenclature (H5 = HS 2017) and whether WITS provides ad valorem
# equivalents of specific duties ("Yes"/"No").
wits_availability <- function(y) {
  doc <- xml2::read_xml(file.path(source_dir,
                                  paste0("wits_availability_", y, ".xml")))
  ns <- xml2::xml_ns(doc)
  nodes <- xml2::xml_find_all(doc, ".//wits:reporter", ns)
  field <- function(name) {
    xml2::xml_text(xml2::xml_find_first(nodes, paste0("wits:", name), ns))
  }
  nomenclature <- xml2::xml_find_first(nodes, "wits:reporternernomenclature", ns)
  d <- data.table(
    wits_iso3 = xml2::xml_attr(nodes, "iso3Code"),
    wits_code = xml2::xml_attr(nodes, "countrycode"),
    year = as.integer(field("year")),
    nomenclature = xml2::xml_attr(nomenclature, "reporternernomenclaturecode"),
    ave_available = field("isspecificdutyexpressionestimatedavailable")
  )
  unique(d[year == y])
}

# 2018 imports of every country (thousand USD), from positive BACI flows.
baci <- fread(file.path(source_dir, "BACI_HS17_Y2018_V202601.csv"),
              select = c("j", "v"))
imports <- baci[is.finite(v) & v > 0, .(imports_2018 = sum(v)), by = j]
codes <- fread(file.path(source_dir, "country_codes_V202601.csv"))
imports[, iso3 := codes$country_iso3[match(j, codes$country_code)]]
imports <- imports[order(-imports_2018)]
rm(baci)

# The EU members applied a common external tariff, reported by WITS under EUN.
eu28 <- c("AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA",
          "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD",
          "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE", "GBR")
imports[, wits_iso3 := ifelse(iso3 %in% eu28, "EUN", iso3)]

# Screening of the importers.
gravity <- as.data.table(readRDS(file.path(source_dir, "Gravity_V202211.rds")))
gravity_ok <- gravity[year %in% 2018:2019 & !is.na(iso3_d),
                      .(n_years = uniqueN(year)), by = iso3_d][n_years == 2, iso3_d]
rm(gravity)
wdi_ok <- Reduce(intersect, list(
  complete_countries("wdi_gdp_2018_2019.json", 2018:2019),
  complete_countries("wdi_gdp_per_capita_2018_2019.json", 2018:2019),
  complete_countries("wdi_population_2018_2019.json", 2018:2019)))
lpi_ok <- complete_countries("lpi_2018.json", 2018)
wits_2019 <- wits_availability(2019)
wits_ok <- wits_2019[!is.na(nomenclature) & nomenclature != "", wits_iso3]

candidates <- imports[!is.na(iso3) & iso3 %in% wdi_ok & iso3 %in% lpi_ok &
                        iso3 %in% gravity_ok & wits_iso3 %in% wits_ok]
candidates <- candidates[1:35, .(iso3, wits_iso3, imports_2018)]

# The 2018 schedule of each candidate.
wits_2018 <- wits_availability(2018)
candidates <- merge(candidates,
                    wits_2018[, .(wits_iso3, wits_code, nomenclature, ave_available)],
                    by = "wits_iso3", all.x = TRUE, sort = FALSE)
fwrite(candidates, file.path(source_dir, "tariff_schedules.csv"))

# Download each schedule once: the most-favoured-nation (MFN) applied tariff
# of every HS6 product, with specific duties converted to ad valorem
# equivalents where WITS provides them ("aveestimated"), and as reported
# otherwise.
schedules <- unique(candidates[, .(wits_iso3, wits_code, ave_available)])
for (i in seq_len(nrow(schedules))) {
  type <- if (schedules$ave_available[i] == "Yes") "aveestimated" else "reported"
  url <- paste0("https://wits.worldbank.org/API/V1/SDMX/V21/datasource/TRN/",
                "reporter/", schedules$wits_code[i],
                "/partner/000/product/all/year/2018/datatype/", type,
                "?format=JSON")
  file <- file.path(tariff_dir, sprintf("%s_%s_%s_2018.json", schedules$wits_iso3[i],
                                        schedules$wits_code[i], type))
  if (!file.exists(file)) download.file(url, file, mode = "wb")
}
message(nrow(candidates), " candidate importers, ", nrow(schedules),
        " tariff schedules.")
