# Step 3 of the trade data preparation: build the response tensor, the design
# and the cross-validation folds of the trade application.
#
#   Rscript reproducibility/real_data/3_build_data.R <data_dir>
#
# Reads the files downloaded by steps 1 and 2 from <data_dir>/sources and
# writes to <data_dir>/prepared:
#   response_Y.csv  one row per exporter x partner x HS section (80 x 25 x 21),
#                   with the 8 response channels of 2019
#   design_X.csv    one row per exporter, with the 6 features of 2018
#   folds.csv       the fold (1 to 5) of each exporter
#
# The sample (exporters and partners), the features and the weights use only
# 2018 data; 2019 data enter only the response. Trade values are in thousand
# USD and quantities in tonnes (BACI).

library(data.table)

data_dir <- commandArgs(trailingOnly = TRUE)[1]
source_dir <- file.path(data_dir, "sources")
out_dir <- file.path(data_dir, "prepared")
dir.create(out_dir, showWarnings = FALSE)

# The 21 sections of the Harmonized System, as ranges of HS2 chapters.
sections <- data.table(
  section = sprintf("HS%02d", 1:21),
  first_chapter = c(1, 6, 15, 16, 25, 28, 39, 41, 44, 47, 50,
                    64, 68, 71, 72, 84, 86, 90, 93, 94, 97)
)
section_of_chapter <- function(chapter) {
  index <- findInterval(chapter, c(sections$first_chapter, 98))
  ifelse(index >= 1 & index <= 21 & chapter <= 97, sections$section[index],
         NA_character_)
}

# ---------------------------------------------------------------------------
# 1. Tariffs of the candidate importers, per HS6 product (2018)
# ---------------------------------------------------------------------------

# One WITS file (SDMX-JSON): the MFN tariff of every HS6 product of one
# schedule. Each series is one product; its key "a:b:c:d" holds the indices
# of reporter, partner, product and data type, and only the product index
# varies within a file.
read_wits_tariffs <- function(file) {
  json <- jsonlite::fromJSON(file, simplifyVector = FALSE)
  dimensions <- json$structure$dimensions$series
  dimension_ids <- vapply(dimensions, `[[`, character(1), "id")
  products <- dimensions[[match("PRODUCTCODE", dimension_ids)]]$values
  attributes <- json$structure$attributes$observation
  attribute_ids <- vapply(attributes, `[[`, character(1), "id")
  # Value of an observation attribute (stored as an index into its values).
  attribute <- function(observation, id) {
    position <- match(id, attribute_ids)
    index <- observation[[position + 1]]
    attributes[[position]]$values[[index + 1]]$id
  }

  series <- json$dataSets[[1]]$series
  keys <- do.call(rbind, lapply(strsplit(names(series), ":"), as.integer))
  product_key <- which(apply(keys, 2, max) > 0)
  rows <- lapply(seq_along(series), function(i) {
    observation <- series[[i]]$observations[[1]]
    data.table(hs6 = products[[keys[i, product_key] + 1]]$id,
               tariff = as.numeric(observation[[1]]),
               nomenclature = attribute(observation, "NOMENCODE"))
  })
  rbindlist(rows)
}

candidates <- fread(file.path(source_dir, "tariff_schedules.csv"),
                    colClasses = list(character = "wits_code"))
schedule_files <- list.files(file.path(source_dir, "wits_tariffs"),
                             full.names = TRUE)
tariffs <- rbindlist(lapply(schedule_files, function(file) {
  d <- read_wits_tariffs(file)
  d[, wits_iso3 := sub("_.*", "", basename(file))]
  d[, type := ifelse(grepl("aveestimated", file), "aveestimated", "reported")]
  d
}))

# One row per candidate importer and product (EU members share the EU schedule).
tariffs <- merge(candidates[, .(iso3, wits_iso3)], tariffs, by = "wits_iso3",
                 allow.cartesian = TRUE)
tariffs[, hs6_code := as.integer(hs6)]
tariffs[, section := section_of_chapter(as.integer(substr(hs6, 1, 2)))]

# ---------------------------------------------------------------------------
# 2. BACI trade flows, summed over the HS6 products of each section
# ---------------------------------------------------------------------------

codes <- fread(file.path(source_dir, "country_codes_V202601.csv"))
codes <- codes[nchar(country_iso3) == 3]

# Flows of one year between different countries: value and quantity by
# exporter, importer and section. Also returns the imports of the
# `hs6_importers` by HS6 product, used to weight the tariffs.
read_baci <- function(year, hs6_importers = NULL) {
  d <- fread(file.path(source_dir, sprintf("BACI_HS17_Y%d_V202601.csv", year)),
             select = c("i", "j", "k", "v", "q"),
             colClasses = list(integer = c("i", "j", "k"), numeric = c("v", "q")))
  d[, exporter := codes$country_iso3[match(i, codes$country_code)]]
  d[, importer := codes$country_iso3[match(j, codes$country_code)]]
  d[, section := section_of_chapter(k %/% 10000L)]
  d <- d[!is.na(exporter) & !is.na(importer) & exporter != importer &
           is.finite(v) & v >= 0 & !is.na(section)]
  # A missing or negative quantity counts as zero.
  by_section <- d[, .(value = sum(v, na.rm = TRUE),
                      quantity = sum(fifelse(is.finite(q) & q >= 0, q, 0),
                                     na.rm = TRUE)),
                  by = .(exporter, importer, section)]
  by_hs6 <- d[importer %in% hs6_importers,
              .(value = sum(v, na.rm = TRUE)), by = .(importer, hs6_code = k, section)]
  list(by_section = by_section, by_hs6 = by_hs6)
}

baci_2018 <- read_baci(2018, hs6_importers = sort(unique(tariffs$iso3)))
trade_2018 <- baci_2018$by_section

# ---------------------------------------------------------------------------
# 3. Country covariates of 2018: WDI, LPI and CEPII Gravity
# ---------------------------------------------------------------------------

world_bank_2018 <- function(file) {
  json <- jsonlite::fromJSON(file.path(source_dir, file))
  d <- data.table(iso3 = json[[2]]$countryiso3code,
                  year = as.integer(json[[2]]$date),
                  value = as.numeric(json[[2]]$value))
  unique(d[year == 2018 & !is.na(iso3) & iso3 != "", .(iso3, value)], by = "iso3")
}
gdp <- world_bank_2018("wdi_gdp_2018_2019.json")
population <- world_bank_2018("wdi_population_2018_2019.json")
lpi <- world_bank_2018("lpi_2018.json")
has_gdp <- intersect(gdp[is.finite(value) & value > 0, iso3],
                     population[is.finite(value) & value > 0, iso3])
has_lpi <- lpi[is.finite(value), iso3]

# Gravity has several rows for some ISO3 pairs (historical entities, for
# example the former DEU, IDN and VNM). Keep one row per pair: one with a
# positive distance, preferring entities that still exist. The distance is
# the harmonic population-weighted one, or the arithmetic or simple distance
# when it is missing.
gravity <- as.data.table(readRDS(file.path(source_dir, "Gravity_V202211.rds")))
gravity <- gravity[year == 2018, .(iso3_o, iso3_d, country_exists_o,
                                   country_exists_d, distw_harmonic,
                                   distw_arithmetic, dist, rta_coverage)]
gravity[, distance := fcoalesce(as.numeric(distw_harmonic),
                                as.numeric(distw_arithmetic), as.numeric(dist))]
gravity[, valid_distance := is.finite(distance) & distance > 0]
gravity[, existing := (country_exists_o == 1) + (country_exists_d == 1)]
setorder(gravity, iso3_o, iso3_d, -valid_distance, -existing)
gravity <- gravity[!duplicated(paste(iso3_o, iso3_d))]

# ---------------------------------------------------------------------------
# 4. Partners: the 25 largest importers among the candidates with WDI, LPI
#    and a tariff in all 21 sections
# ---------------------------------------------------------------------------

full_tariff <- tariffs[is.finite(tariff), .(n = uniqueN(section)), by = iso3][n == 21, iso3]
partners <- trade_2018[importer %in% candidates$iso3,
                       .(imports = sum(value)), by = importer]
partners[, eligible := importer %in% has_gdp & importer %in% has_lpi &
           importer %in% full_tariff]
setorder(partners, -eligible, -imports, importer)
P <- partners[eligible == TRUE, importer][1:25]

# ---------------------------------------------------------------------------
# 5. Exporters: the 80 largest exporters to the partners that trade with at
#    least 13 partners in at least 11 sections, have WDI data and a Gravity
#    distance to all 25 partners
# ---------------------------------------------------------------------------

exporters <- trade_2018[importer %in% P & value > 0,
                        .(n_partners = uniqueN(importer),
                          n_sections = uniqueN(section),
                          exports = sum(value)), by = exporter]
with_distance <- gravity[iso3_d %in% P & valid_distance,
                         .(n = uniqueN(iso3_d)), by = iso3_o][n == 25, iso3_o]
exporters[, eligible := n_partners >= 13 & n_sections >= 11 &
            exporter %in% has_gdp & exporter %in% with_distance]
setorder(exporters, -eligible, -exports, exporter)
exporters <- exporters[eligible == TRUE][1:80]
E <- exporters$exporter

# ---------------------------------------------------------------------------
# 6. Tariff of each partner and section: the mean of the HS6 tariffs weighted
#    by the partner's 2018 imports. The unweighted mean is used for schedules
#    not in HS 2017 (nomenclature other than H5) and when the products with
#    imports cover less than half of the section's imports.
# ---------------------------------------------------------------------------

partner_tariffs <- tariffs[iso3 %in% P & is.finite(tariff) & tariff >= 0]
partner_tariffs <- merge(partner_tariffs,
                         baci_2018$by_hs6[, .(iso3 = importer, hs6_code, value)],
                         by = c("iso3", "hs6_code"), all.x = TRUE)
partner_tariffs[is.na(value), value := 0]
section_imports <- baci_2018$by_hs6[importer %in% P,
                                    .(section_imports = sum(value)),
                                    by = .(iso3 = importer, section)]
section_tariffs <- partner_tariffs[, .(
  nomenclature = paste(sort(unique(nomenclature)), collapse = "|"),
  mean = mean(tariff),
  weighted_mean = if (sum(value) > 0) weighted.mean(tariff, value) else NA_real_,
  matched_imports = sum(value)
), by = .(iso3, section)]
section_tariffs <- merge(section_tariffs, section_imports,
                         by = c("iso3", "section"), all.x = TRUE)
section_tariffs[, tariff := ifelse(
  nomenclature == "H5" & is.finite(weighted_mean) &
    matched_imports / pmax(section_imports, 1e-12) >= 0.5,
  weighted_mean, mean)]

# ---------------------------------------------------------------------------
# 7. One row per exporter, partner and section: 2018 trade and features
# ---------------------------------------------------------------------------

cells <- CJ(exporter = E, importer = P, section = sections$section,
            sorted = FALSE)
cells <- merge(cells, trade_2018[, .(exporter, importer, section,
                                     value_2018 = value)],
               by = c("exporter", "importer", "section"), all.x = TRUE)
cells[is.na(value_2018), value_2018 := 0]

# Features of the cell (all from 2018):
#   tariff           MFN tariff of the partner in the section (percent)
#   ntb_proxy        minus the partner's LPI score (logistics friction)
#   log_gdp_partner  log GDP of the partner (current USD)
#   log_distance     log distance between exporter and partner (km)
#   rta_scope_count  trade agreement in force: 0 none, 1 goods or services,
#                    2 goods and services
#   world_demand     log(1 + imports of the partner in the section)
demand <- trade_2018[importer %in% P, .(world_demand = log1p(sum(value))),
                     by = .(importer, section)]
cells <- merge(cells, demand, by = c("importer", "section"), all.x = TRUE)
cells <- merge(cells, section_tariffs[, .(importer = iso3, section, tariff)],
               by = c("importer", "section"), all.x = TRUE)
cells <- merge(cells, lpi[, .(importer = iso3, ntb_proxy = -value)],
               by = "importer", all.x = TRUE)
cells <- merge(cells, gdp[, .(importer = iso3, log_gdp_partner = log(value))],
               by = "importer", all.x = TRUE)
pairs <- gravity[iso3_o %in% E & iso3_d %in% P]
pairs[, coverage := trimws(as.character(rta_coverage))]
pairs[, rta_scope_count := as.numeric(coverage %in% c("Goods", "Goods & Services")) +
        as.numeric(coverage %in% c("Services", "Goods & Services"))]
cells <- merge(cells, pairs[, .(exporter = iso3_o, importer = iso3_d,
                                log_distance = log(distance), rta_scope_count)],
               by = c("exporter", "importer"), all.x = TRUE)

# Order the cells by exporter, then partner, then section.
cells <- cells[order(match(exporter, E), match(importer, P),
                     match(section, sections$section))]
features <- c("tariff", "ntb_proxy", "log_gdp_partner", "log_distance",
              "rta_scope_count", "world_demand")

# ---------------------------------------------------------------------------
# 8. Response (2019). Shares use the totals over all countries in BACI.
# ---------------------------------------------------------------------------

trade_2019 <- read_baci(2019)$by_section
cells <- merge(cells, trade_2019[, .(exporter, importer, section,
                                     value = value, quantity = quantity)],
               by = c("exporter", "importer", "section"), all.x = TRUE, sort = FALSE)
cells[is.na(value), value := 0]
cells[is.na(quantity), quantity := 0]

exporter_totals <- trade_2019[exporter %in% E, .(exporter_total = sum(value)),
                              by = exporter]
importer_section_totals <- trade_2019[importer %in% P,
                                      .(importer_section_total = sum(value)),
                                      by = .(importer, section)]
exporter_section_totals <- trade_2019[exporter %in% E,
                                      .(exporter_section_total = sum(value)),
                                      by = .(exporter, section)]
section_totals <- trade_2019[, .(section_total = sum(value)), by = section]
world_total <- trade_2019[, sum(value)]
cells <- merge(cells, exporter_totals, by = "exporter", all.x = TRUE, sort = FALSE)
cells <- merge(cells, importer_section_totals, by = c("importer", "section"),
               all.x = TRUE, sort = FALSE)
cells <- merge(cells, exporter_section_totals, by = c("exporter", "section"),
               all.x = TRUE, sort = FALSE)
cells <- merge(cells, section_totals, by = "section", all.x = TRUE, sort = FALSE)
cells[is.na(exporter_section_total), exporter_section_total := 0]
cells <- cells[order(match(exporter, E), match(importer, P),
                     match(section, sections$section))]

# Response channels of each cell (value v, quantity q of 2019):
#   log_export_value      log(1 + v)
#   log_export_quantity   log(1 + q)
#   log_unit_value        log(v / q) if v, q > 0, else 0
#   exporter_share        v / exports of the exporter to all countries
#   partner_market_share  v / imports of the partner in the section
#   rca_balassa           Balassa index of the exporter in the section:
#                         (exporter's exports in the section / its exports) /
#                         (world exports in the section / world exports)
#   extensive_margin      1 if v > 0, else 0
#   intensive_margin      log(v) if v > 0, else 0
eps <- 1e-12
cells[, log_export_value := log1p(value)]
cells[, log_export_quantity := log1p(quantity)]
cells[, log_unit_value := ifelse(value > 0 & quantity > 0, log(value / quantity), 0)]
cells[, exporter_share := value / pmax(exporter_total, eps)]
cells[, partner_market_share := value / pmax(importer_section_total, eps)]
cells[, rca_balassa := (exporter_section_total / pmax(exporter_total, eps)) /
        pmax(section_total / pmax(world_total, eps), eps)]
cells[, extensive_margin := as.numeric(value > 0)]
cells[, intensive_margin := ifelse(value > 0, log(pmax(value, eps)), 0)]
channels <- c("log_export_value", "log_export_quantity", "log_unit_value",
              "exporter_share", "partner_market_share", "rca_balassa",
              "extensive_margin", "intensive_margin")

# ---------------------------------------------------------------------------
# 9. Design: for each exporter, the average of each feature over its cells,
#    weighted by the square root of the cell's 2018 trade value
# ---------------------------------------------------------------------------

design <- cells[, lapply(.SD, function(z) {
  w <- sqrt(value_2018)
  if (sum(w) > 0) sum(z * w) / sum(w) else mean(z)
}), by = exporter, .SDcols = features]

# ---------------------------------------------------------------------------
# 10. Folds: five folds of 16 exporters, balanced over the quintiles of
#     log(1 + 2018 exports to the partners)
# ---------------------------------------------------------------------------

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(2549)
size <- log1p(exporters$exports)
quintile <- cut(size, breaks = quantile(size, probs = seq(0, 1, 0.2)),
                include.lowest = TRUE, labels = FALSE)
fold <- integer(80)
for (s in 1:5) {
  members <- which(quintile == s)
  fold[members] <- sample((seq_along(members) + s - 2) %% 5 + 1)
}

# ---------------------------------------------------------------------------
# 11. Write
# ---------------------------------------------------------------------------

response <- cells[, c("exporter", "importer", "section", channels), with = FALSE]
setnames(response, "importer", "partner")
for (channel in channels) {
  set(response, j = channel, value = sprintf("%.17g", response[[channel]]))
}
fwrite(response, file.path(out_dir, "response_Y.csv"), quote = FALSE)
write.csv(design, file.path(out_dir, "design_X.csv"), row.names = FALSE)
write.csv(data.frame(exporter = E, fold = fold),
          file.path(out_dir, "folds.csv"), row.names = FALSE)
message("Wrote ", out_dir, ": ", length(E), " exporters, ", length(P),
        " partners, 21 sections.")
