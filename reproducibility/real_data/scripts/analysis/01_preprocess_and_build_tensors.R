#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE)
code_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
data_dir <- Sys.getenv("SEGRECP_REAL_DATA_DIR", unset = "")
if (!nzchar(data_dir)) stop("SEGRECP_REAL_DATA_DIR must name a private data directory")
official_root <- normalizePath(data_dir, mustWork = TRUE)
analysis_root <- file.path(official_root, "analysis_official")
Sys.setenv(TT_ANALYSIS_ROOT = analysis_root)
source(file.path(code_root, "scripts", "_utils.R"))

stage <- "official_preprocess"
raw_official <- file.path(official_root, "raw_official")
processed_official <- file.path(official_root, "processed_official")
metadata_official <- file.path(official_root, "metadata")

save_dt <- function(x, target, extra = list()) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, target)
  .write_meta(target, c(extra, list(rows = nrow(x), cols = ncol(x))))
  invisible(x)
}

hs_lo <- HS_SECTIONS$hs2_lo
map_hs_section <- function(hs2) {
  idx <- findInterval(as.integer(hs2), c(hs_lo, 98L))
  out <- rep(NA_character_, length(idx))
  ok <- idx >= 1L & idx <= 21L & as.integer(hs2) <= 97L
  out[ok] <- HS_SECTIONS$hs1[idx[ok]]
  out
}

required <- c(
  file.path(raw_official, "BACI_HS17_Y2018_V202601.csv"),
  file.path(raw_official, "BACI_HS17_Y2019_V202601.csv"),
  file.path(raw_official, "country_codes_V202601.csv"),
  file.path(raw_official, "Gravity_V202211.rds"),
  file.path(processed_official, "wits_mfn_tariffs_by_importer_2018.csv"),
  file.path(processed_official, "wdi_macro_2018.csv"),
  file.path(processed_official, "lpi_logistics_friction_2018.csv"),
  file.path(processed_official, "cepii_gravity_2018.csv")
)
if (any(!file.exists(required))) {
  stop("Missing official input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

# Keep a compact, auditable raw layer without duplicating the 700+ MB BACI files.
manifest <- fread(file.path(metadata_official, "source_manifest_sha256.csv"))
needed_manifest <- manifest[local_path %in% c(
  "raw_official/BACI_HS17_Y2018_V202601.csv",
  "raw_official/BACI_HS17_Y2019_V202601.csv",
  "raw_official/country_codes_V202601.csv",
  "raw_official/Gravity_V202211.rds",
  "raw_official/wdi_gdp_current_usd_2018_2019.json",
  "raw_official/wdi_gdp_pc_current_usd_2018_2019.json",
  "raw_official/wdi_population_2018_2019.json",
  "raw_official/world_bank_lpi_overall_2018.json"
) | grepl("raw_official/wits_tariffs_2018_", local_path)]
needed_manifest[, resolved_path := file.path(official_root, local_path)]
needed_manifest[, exists_now := file.exists(resolved_path)]
needed_manifest[, size_matches_manifest :=
                  exists_now & file.info(resolved_path)$size == bytes]
save_dt(needed_manifest, file.path(paths$raw, "source_pointers.csv"),
        list(note = "No raw file is copied; resolved_path points to the official bundle"))
if (needed_manifest[, any(synthetic != FALSE | !exists_now | !size_matches_manifest)]) {
  stop("Official source provenance or size validation failed")
}

country_codes <- fread(file.path(raw_official, "country_codes_V202601.csv"))
country_codes <- country_codes[nchar(country_iso3) == 3L]

read_baci_section <- function(year, keep_exporters = NULL, keep_importers = NULL,
                              keep_hs6_importers = NULL) {
  input <- file.path(raw_official, sprintf("BACI_HS17_Y%d_V202601.csv", year))
  log_msg(sprintf("Reading official BACI %d (%s) ...", year, basename(input)), stage = stage)
  d <- fread(input, select = c("i", "j", "k", "v", "q"),
             colClasses = list(integer = c("i", "j", "k"),
                               numeric = c("v", "q")), showProgress = TRUE)
  n_in <- nrow(d)
  d[, `:=`(
    exp_iso3 = country_codes$country_iso3[match(i, country_codes$country_code)],
    imp_iso3 = country_codes$country_iso3[match(j, country_codes$country_code)],
    hs2 = k %/% 10000L
  )]
  d[, hs1 := map_hs_section(hs2)]
  n_unmapped_country <- d[is.na(exp_iso3) | is.na(imp_iso3), .N]
  n_self <- d[!is.na(exp_iso3) & exp_iso3 == imp_iso3, .N]
  n_bad_value <- d[!is.finite(v) | v < 0, .N]
  n_bad_section <- d[is.na(hs1), .N]
  n_quantity_missing <- d[!is.finite(q) | q < 0, .N]
  d <- d[!is.na(exp_iso3) & !is.na(imp_iso3) & exp_iso3 != imp_iso3 &
           is.finite(v) & v >= 0 & !is.na(hs1)]
  if (!is.null(keep_exporters)) d <- d[exp_iso3 %in% keep_exporters]
  if (!is.null(keep_importers)) d <- d[imp_iso3 %in% keep_importers]

  hs6_trade <- NULL
  if (!is.null(keep_hs6_importers)) {
    hs6_trade <- d[imp_iso3 %in% keep_hs6_importers,
                   .(trade_value_kusd = sum(v, na.rm = TRUE)),
                   by = .(imp_iso3, hs6_num = k, hs1)]
  }
  section <- d[, .(
    value_kusd = sum(v, na.rm = TRUE),
    quantity_t = sum(fifelse(is.finite(q) & q >= 0, q, 0), na.rm = TRUE),
    n_hs6_rows = .N,
    n_quantity_missing = sum(!is.finite(q) | q < 0)
  ), by = .(exp_iso3, imp_iso3, hs1)]
  qa <- data.table(
    year = year,
    rows_input = n_in,
    rows_used = nrow(d),
    unmapped_country_rows = n_unmapped_country,
    self_trade_rows = n_self,
    invalid_value_rows = n_bad_value,
    invalid_hs_section_rows = n_bad_section,
    missing_or_invalid_quantity_rows = n_quantity_missing,
    section_rows = nrow(section),
    total_value_kusd = section[, sum(value_kusd)]
  )
  rm(d); gc()
  list(section = section, hs6_trade = hs6_trade, qa = qa)
}

tariff_lines <- fread(file.path(processed_official,
                                "wits_mfn_tariffs_by_importer_2018.csv"),
                      colClasses = list(character = c("imp_iso3", "hs6", "hs_section")))
tariff_importers <- sort(unique(tariff_lines$imp_iso3))

b18 <- read_baci_section(2018L, keep_hs6_importers = tariff_importers)
trade18 <- b18$section
hs6_trade18 <- b18$hs6_trade

macro <- fread(file.path(processed_official, "wdi_macro_2018.csv"))
lpi <- fread(file.path(processed_official, "lpi_logistics_friction_2018.csv"))
# Rebuild the analysis gravity table from the official RDS. The earlier compact
# helper table used unique(..., by=ISO3 pair), which can retain an obsolete
# entity row with missing distance when historical country IDs share an ISO3
# code (notably DEU, IDN, and VNM). We select a finite-distance, currently
# existing entity row before ISO3-pair deduplication.
gravity_full <- as.data.table(readRDS(file.path(raw_official, "Gravity_V202211.rds")))
gravity <- gravity_full[year == 2018L, .(
  iso3_o, iso3_d, country_exists_o, country_exists_d,
  distw_harmonic, distw_arithmetic, dist, contig, comlang_off, comcol, rta_type
  , rta_coverage
)]
rm(gravity_full); gc()
gravity[, distw_h := fcoalesce(as.numeric(distw_harmonic),
                               as.numeric(distw_arithmetic), as.numeric(dist))]
gravity[, distance_source := fifelse(is.finite(distw_harmonic) & distw_harmonic > 0,
                                     "distw_harmonic",
                              fifelse(is.finite(distw_arithmetic) & distw_arithmetic > 0,
                                      "distw_arithmetic", "dist"))]
gravity[, valid_distance := is.finite(distw_h) & distw_h > 0]
gravity[, entity_score := fifelse(country_exists_o == 1, 1L, 0L) +
                            fifelse(country_exists_d == 1, 1L, 0L)]
gravity_duplicates_before <- gravity[, sum(duplicated(paste(iso3_o, iso3_d)))]
setorder(gravity, iso3_o, iso3_d, -valid_distance, -entity_score)
gravity <- gravity[!duplicated(paste(iso3_o, iso3_d))]
gravity_qa <- data.table(
  rows_after_year_filter = nrow(gravity),
  duplicate_iso3_rows_resolved = gravity_duplicates_before,
  missing_distance_after_resolution = gravity[, sum(!valid_distance)],
  finite_distance_rows = gravity[, sum(valid_distance)],
  harmonic_rows = gravity[, sum(valid_distance & distance_source == "distw_harmonic")],
  arithmetic_fallback_rows = gravity[, sum(valid_distance & distance_source == "distw_arithmetic")],
  simple_distance_fallback_rows = gravity[, sum(valid_distance & distance_source == "dist")]
)
save_dt(gravity[, .(year = 2018L, iso3_o, iso3_d, distw_h, distance_source,
                     contig, comlang_off, comcol, rta_type, rta_coverage)],
        file.path(paths$processed, "gravity_clean_for_analysis_2018.csv"),
        list(dedup_rule = "finite distance then current-entity score before ISO3 pair deduplication"))
save_dt(gravity_qa, file.path(paths$qc, "gravity_resolution_qa.csv"))
macro_ok <- macro[is.finite(gdp_usd) & gdp_usd > 0 & is.finite(pop) & pop > 0, iso3]
lpi_ok <- lpi[is.finite(lpi_score), iso3]
tariff_ok <- tariff_lines[is.finite(tariff_avg_pct),
                           .(n_sections = uniqueN(hs_section)), by = imp_iso3][n_sections == 21L, imp_iso3]

partner_log <- trade18[imp_iso3 %in% tariff_importers,
                       .(total_imports_kusd_2018 = sum(value_kusd)), by = imp_iso3]
partner_log[, `:=`(
  has_macro = imp_iso3 %in% macro_ok,
  has_lpi = imp_iso3 %in% lpi_ok,
  has_21_tariff_sections = imp_iso3 %in% tariff_ok
)]
partner_log[, eligible := has_macro & has_lpi & has_21_tariff_sections]
setorder(partner_log, -eligible, -total_imports_kusd_2018, imp_iso3)
eligible_partner_log <- partner_log[eligible == TRUE]
partners <- eligible_partner_log[seq_len(min(25L, nrow(eligible_partner_log)))]
if (nrow(partners) != 25L) stop("Fewer than 25 fully covered importer partners")
partners[, rank := seq_len(.N)]
P <- partners$imp_iso3

exp_log <- trade18[imp_iso3 %in% P & value_kusd > 0,
                   .(n_partners_covered = uniqueN(imp_iso3),
                     n_sections_covered = uniqueN(hs1),
                     n_positive_cells = .N,
                     total_export_kusd_2018 = sum(value_kusd)), by = exp_iso3]
grav_stats <- gravity[iso3_d %in% P & is.finite(distw_h) & distw_h > 0,
                      .(n_gravity_partners = uniqueN(iso3_d)), by = iso3_o]
exp_log <- merge(exp_log, grav_stats, by.x = "exp_iso3", by.y = "iso3_o", all.x = TRUE)
exp_log[is.na(n_gravity_partners), n_gravity_partners := 0L]
exp_log[, `:=`(
  has_macro = exp_iso3 %in% macro_ok,
  eligible_coverage = n_partners_covered >= 13L & n_sections_covered >= 11L,
  complete_gravity = n_gravity_partners == 25L
)]
exp_log[, eligible := has_macro & eligible_coverage & complete_gravity]
setorder(exp_log, -eligible, -total_export_kusd_2018, exp_iso3)
eligible_exporters <- exp_log[eligible == TRUE]
exporters <- eligible_exporters[seq_len(min(80L, nrow(eligible_exporters)))]
log_msg(sprintf(paste0("Exporter eligibility: %d candidates; %d pass trade coverage; ",
                       "%d have macro; %d have complete 25-partner gravity; %d pass all."),
                nrow(exp_log), exp_log[, sum(eligible_coverage)],
                exp_log[, sum(has_macro)], exp_log[, sum(complete_gravity)],
                nrow(eligible_exporters)), stage = stage)
if (nrow(exporters) < 30L) {
  save_dt(exp_log, file.path(paths$qc, "exporter_selection_diagnostic.csv"),
          list(status = "failed_before_tensor_build"))
  stop("Too few eligible exporters for five-fold analysis")
}
exporters[, rank := seq_len(.N)]
exp_log[, selected_top80 := exp_iso3 %in% exporters$exp_iso3]
E <- exporters$exp_iso3

save_dt(partner_log, file.path(paths$processed, "partner_selection_log.csv"),
        list(selection_year = 2018L, leakage_guard = "No 2019 outcome used"))
save_dt(partners[, .(imp_iso3, rank, total_imports_kusd_2018)],
        file.path(paths$processed, "partners.csv"),
        list(rule = "top 25 fully covered tariff-schedule importers by 2018 BACI imports"))
save_dt(exp_log, file.path(paths$processed, "exporter_selection_log.csv"),
        list(selection_year = 2018L, leakage_guard = "No 2019 outcome used"))
save_dt(exporters[, .(exp_iso3, rank, n_partners_covered, n_sections_covered,
                       n_positive_cells, n_gravity_partners, total_export_kusd_2018)],
        file.path(paths$processed, "exporters.csv"),
        list(rule = "coverage eligible, then top 80 by 2018 BACI exports"))

log_msg(sprintf("Selected %d partners and %d exporters from 2018 only.",
                length(P), length(E)), stage = stage)

# Tariffs are combined at the common HS-section level. H5 schedules use
# 2018 BACI HS6 trade weights when at least half of section trade is matched;
# H4 (Singapore) and low-match cells use the official simple mean.
tariff_lines[, hs6_num := as.integer(hs6)]
tariff_lines <- tariff_lines[imp_iso3 %in% P & hs_section %in% HS_SECTIONS$hs1 &
                               is.finite(tariff_avg_pct) & tariff_avg_pct >= 0]
tw <- merge(tariff_lines, hs6_trade18[, .(imp_iso3, hs6_num, trade_value_kusd)],
            by = c("imp_iso3", "hs6_num"), all.x = TRUE)
tw[is.na(trade_value_kusd), trade_value_kusd := 0]
section_trade_totals <- hs6_trade18[imp_iso3 %in% P,
                                    .(section_trade_value_kusd = sum(trade_value_kusd)),
                                    by = .(imp_iso3, hs_section = hs1)]
tariff_section <- tw[, .(
  nomenclature = paste(sort(unique(nomenclature)), collapse = "|"),
  source_datatype = paste(sort(unique(source_datatype_request)), collapse = "|"),
  n_tariff_lines = .N,
  n_trade_matched_lines = sum(trade_value_kusd > 0),
  simple_mean_tariff_pct = mean(tariff_avg_pct),
  weighted_mean_tariff_pct = if (sum(trade_value_kusd) > 0)
    weighted.mean(tariff_avg_pct, trade_value_kusd) else NA_real_,
  matched_trade_value_kusd = sum(trade_value_kusd),
  max_line_tariff_pct = max(tariff_avg_pct)
), by = .(imp_iso3, hs_section)]
tariff_section <- merge(tariff_section, section_trade_totals,
                        by = c("imp_iso3", "hs_section"), all.x = TRUE)
tariff_section[, matched_trade_share := matched_trade_value_kusd /
                  pmax(section_trade_value_kusd, 1e-12)]
tariff_section[, aggregation_rule := fifelse(
  nomenclature == "H5" & is.finite(weighted_mean_tariff_pct) & matched_trade_share >= 0.5,
  "2018_BACI_HS6_trade_weighted", "official_schedule_simple_mean_at_HS_section")]
tariff_section[, tariff_avg_pct := fifelse(
  aggregation_rule == "2018_BACI_HS6_trade_weighted",
  weighted_mean_tariff_pct, simple_mean_tariff_pct)]
if (tariff_section[, uniqueN(paste(imp_iso3, hs_section))] != 25L * 21L ||
    tariff_section[, any(!is.finite(tariff_avg_pct))]) {
  stop("Tariff section table is incomplete")
}
save_dt(tariff_section, file.path(paths$processed, "tariff_section_2018.csv"),
        list(tariff_year = 2018L,
             singapore_harmonization = "H4 and H5 combined only after HS2-to-section mapping"))

# The 2019 read occurs only after E and P have been frozen from 2018. We retain
# the full 2019 section table long enough to form standard world denominators
# for exporter shares, partner market shares, and Balassa RCA; only response
# numerators are then restricted to the selected E x P cube.
b19 <- read_baci_section(2019L)
trade19_all <- b19$section
trade19 <- trade19_all[exp_iso3 %in% E & imp_iso3 %in% P]
baci_qa <- rbind(b18$qa, b19$qa, fill = TRUE)
save_dt(baci_qa, file.path(paths$qc, "baci_processing_qa.csv"),
        list(note = "2019 rows_used retain the full cross-section for standard world response denominators; E x P response numerators are restricted afterward"))

# Dense selected cube, retaining structural zeros.
grid <- CJ(exp_iso3 = E, imp_iso3 = P, hs1 = HS_SECTIONS$hs1, sorted = FALSE)
lag_panel <- trade18[exp_iso3 %in% E & imp_iso3 %in% P,
                     .(exp_iso3, imp_iso3, hs1,
                       lag_value_kusd = value_kusd,
                       lag_quantity_t = quantity_t)]
out_panel <- trade19[, .(exp_iso3, imp_iso3, hs1,
                         value_kusd, quantity_t, n_quantity_missing)]
panel <- merge(grid, lag_panel, by = c("exp_iso3", "imp_iso3", "hs1"), all.x = TRUE)
panel <- merge(panel, out_panel, by = c("exp_iso3", "imp_iso3", "hs1"), all.x = TRUE)
for (v in c("lag_value_kusd", "lag_quantity_t", "value_kusd", "quantity_t",
            "n_quantity_missing")) set(panel, which(is.na(panel[[v]])), v, 0)

demand <- trade18[imp_iso3 %in% P,
                  .(world_demand = log1p(sum(value_kusd))), by = .(imp_iso3, hs1)]
tar_small <- tariff_section[, .(imp_iso3, hs1 = hs_section,
                                tariff = tariff_avg_pct, tariff_aggregation_rule = aggregation_rule)]
macro_small <- macro[iso3 %in% P, .(imp_iso3 = iso3, log_gdp_partner = log(gdp_usd))]
lpi_small <- lpi[iso3 %in% P,
                 .(imp_iso3 = iso3, lpi_score, logistics_friction_proxy)]
grav_small <- gravity[iso3_o %in% E & iso3_d %in% P,
                      .(exp_iso3 = iso3_o, imp_iso3 = iso3_d,
                        distw_h, log_distance = log(distw_h),
                        rta_type_label = trimws(as.character(rta_type)),
                        rta_coverage_label = trimws(as.character(rta_coverage)),
                        rta_indicator = as.numeric(!is.na(rta_type) &
                          nzchar(trimws(as.character(rta_type)))),
                        rta_goods = as.numeric(trimws(as.character(rta_coverage)) %in%
                          c("Goods", "Goods & Services")),
                        rta_services = as.numeric(trimws(as.character(rta_coverage)) %in%
                          c("Services", "Goods & Services")),
                        rta_scope_count = as.numeric(trimws(as.character(rta_coverage)) %in%
                          c("Goods", "Goods & Services")) +
                          as.numeric(trimws(as.character(rta_coverage)) %in%
                          c("Services", "Goods & Services")),
                        no_rta_recorded = is.na(rta_type) |
                          !nzchar(trimws(as.character(rta_type))))]
panel <- merge(panel, demand, by = c("imp_iso3", "hs1"), all.x = TRUE)
panel <- merge(panel, tar_small, by = c("imp_iso3", "hs1"), all.x = TRUE)
panel <- merge(panel, macro_small, by = "imp_iso3", all.x = TRUE)
panel <- merge(panel, lpi_small, by = "imp_iso3", all.x = TRUE)
panel <- merge(panel, grav_small, by = c("exp_iso3", "imp_iso3"), all.x = TRUE)
panel[, `:=`(exp_order = match(exp_iso3, E),
             imp_order = match(imp_iso3, P),
             hs_order = match(hs1, HS_SECTIONS$hs1))]
setorder(panel, exp_order, imp_order, hs_order)
panel[, c("exp_order", "imp_order", "hs_order") := NULL]

feature_cols <- c("tariff", "logistics_friction_proxy", "log_gdp_partner",
                  "log_distance", "rta_indicator", "world_demand")
missingness <- data.table(variable = names(panel),
                          n_missing = vapply(panel, function(x) sum(is.na(x)), integer(1L)),
                          pct_missing = 100 * vapply(panel, function(x) mean(is.na(x)), numeric(1L)))
save_dt(missingness, file.path(paths$qc, "panel_missingness.csv"))
if (panel[, any(!is.finite(as.matrix(.SD))), .SDcols = feature_cols]) {
  stop("Non-finite feature remains in the selected real panel")
}
if (anyDuplicated(panel[, paste(exp_iso3, imp_iso3, hs1)])) stop("Duplicate panel key")
if (nrow(panel) != length(E) * 25L * 21L) stop("Dense panel row-count mismatch")
save_dt(panel, file.path(paths$processed, "panel_real_2018_to_2019.csv"),
        list(selection_year = 2018L, covariate_year = 2018L, outcome_year = 2019L,
             leakage_guard = "2019 appears only in response columns value_kusd and quantity_t"))

# Construct tensors in the manuscript's mode order.
n <- length(E); J <- length(P); K <- 21L; F <- 6L; O <- 8L
idx_m <- match(panel$exp_iso3, E)
idx_j <- match(panel$imp_iso3, P)
idx_k <- match(panel$hs1, HS_SECTIONS$hs1)
lin <- idx_m + (idx_j - 1L) * n + (idx_k - 1L) * n * J
X <- array(NA_real_, dim = c(n, J, K, F),
           dimnames = list(exporter = E, partner = P, hs1 = HS_SECTIONS$hs1,
                           channel = INPUT_CHANNELS))
for (f in seq_len(F)) X[lin + (f - 1L) * n * J * K] <- panel[[feature_cols[f]]]
trade_weights_2018 <- array(0, dim = c(n, J, K))
value_arr <- array(0, dim = c(n, J, K))
qty_arr <- array(0, dim = c(n, J, K))
trade_weights_2018[lin] <- panel$lag_value_kusd
value_arr[lin] <- panel$value_kusd
qty_arr[lin] <- panel$quantity_t

# Standard 2019 response denominators from the full BACI cross-section.
exp_world_dt <- trade19_all[exp_iso3 %in% E,
                            .(export_value_world_kusd_2019 = sum(value_kusd)),
                            by = exp_iso3]
imp_sec_world_dt <- trade19_all[imp_iso3 %in% P,
                                .(import_value_world_kusd_2019 = sum(value_kusd)),
                                by = .(imp_iso3, hs1)]
exp_sec_world_dt <- trade19_all[exp_iso3 %in% E,
                                .(export_section_world_kusd_2019 = sum(value_kusd)),
                                by = .(exp_iso3, hs1)]
sec_world_dt <- trade19_all[, .(world_section_kusd_2019 = sum(value_kusd)), by = hs1]
world_total_2019 <- trade19_all[, sum(value_kusd)]
save_dt(exp_world_dt, file.path(paths$processed, "exporter_world_totals_2019.csv"))
save_dt(imp_sec_world_dt,
        file.path(paths$processed, "importer_section_world_totals_2019.csv"))
save_dt(exp_sec_world_dt,
        file.path(paths$processed, "exporter_section_world_totals_2019.csv"))
save_dt(sec_world_dt, file.path(paths$processed, "world_section_totals_2019.csv"),
        list(world_total_kusd_2019 = world_total_2019))

exp_total <- exp_world_dt$export_value_world_kusd_2019[match(E, exp_world_dt$exp_iso3)]
imp_sec_total <- matrix(0, J, K)
ii <- match(imp_sec_world_dt$imp_iso3, P); kk <- match(imp_sec_world_dt$hs1, HS_SECTIONS$hs1)
imp_sec_total[cbind(ii, kk)] <- imp_sec_world_dt$import_value_world_kusd_2019
exp_sec_total <- matrix(0, n, K)
mm <- match(exp_sec_world_dt$exp_iso3, E); kk <- match(exp_sec_world_dt$hs1, HS_SECTIONS$hs1)
exp_sec_total[cbind(mm, kk)] <- exp_sec_world_dt$export_section_world_kusd_2019
sec_total <- sec_world_dt$world_section_kusd_2019[match(HS_SECTIONS$hs1, sec_world_dt$hs1)]
world_total <- world_total_2019
if (any(!is.finite(c(exp_total, imp_sec_total, exp_sec_total, sec_total, world_total))) ||
    any(c(exp_total, imp_sec_total, sec_total, world_total) <= 0) ||
    any(exp_sec_total < 0)) {
  stop("Invalid full-2019 response denominator")
}

Y <- array(0, dim = c(n, J, K, O),
           dimnames = list(exporter = E, partner = P, hs1 = HS_SECTIONS$hs1,
                           channel = OUTPUT_CHANNELS))
eps <- 1e-12
positive_qty <- value_arr > 0 & qty_arr > 0
Y[,,,1L] <- log1p(value_arr)
Y[,,,2L] <- log1p(qty_arr)
Y[,,,3L] <- ifelse(positive_qty, log(value_arr / qty_arr), 0)
Y[,,,4L] <- value_arr / pmax(array(exp_total, dim = c(n, J, K)), eps)
Y[,,,5L] <- value_arr / pmax(array(rep(imp_sec_total, each = n), dim = c(n, J, K)), eps)
rca_num <- exp_sec_total / pmax(exp_total, eps)
rca_den <- sec_total / pmax(world_total, eps)
rca_mat <- sweep(rca_num, 2L, pmax(rca_den, eps), "/")
for (j in seq_len(J)) Y[,j,,6L] <- rca_mat
Y[,,,7L] <- as.numeric(value_arr > 0)
Y[,,,8L] <- ifelse(value_arr > 0, log(pmax(value_arr, eps)), 0)

assert_shape(X, c(n, J, K, F))
assert_shape(Y, c(n, J, K, O))
if (any(!is.finite(X)) || any(!is.finite(Y))) stop("Non-finite tensor entry")

X_in_raw <- matrix(NA_real_, nrow = n, ncol = F,
                   dimnames = list(exporter = E, channel = INPUT_CHANNELS))
for (f in seq_len(F)) {
  num <- apply(trade_weights_2018 * X[,,,f], 1L, sum)
  den <- apply(trade_weights_2018, 1L, sum)
  fallback <- apply(X[,,,f], 1L, mean)
  X_in_raw[,f] <- ifelse(den > 0, num / den, fallback)
}
Y_block_raw <- matrix(as.numeric(Y), nrow = n, ncol = J * K * O)
dims_full <- c(J, K, O, F)
d_0 <- 3L
sections <- HS_SECTIONS
channels <- list(input = INPUT_CHANNELS, output = OUTPUT_CHANNELS)
meta <- list(data_status = "official_real", selection_year = 2018L,
             covariate_year = 2018L, outcome_year = 2019L,
             n_samples = n, n_partners = J, n_sections = K,
             n_input_channels = F, n_output_channels = O,
             exporter_cap = 80L,
             response_denominators = "full 2019 BACI world totals, not the selected cube",
             leakage_guard = "2018 constructs sample/X/weights; 2019 constructs Y only")

save_rdata(file.path(paths$processed, "tensors_real.RData"),
           X, Y, trade_weights_2018, exporters, partners, sections, channels, meta,
           extra_meta = list(data_status = "official_real",
                             X_shape = paste(dim(X), collapse = "x"),
                             Y_shape = paste(dim(Y), collapse = "x")))
save_rdata(file.path(paths$processed, "regression_inputs_raw.RData"),
           X_in_raw, Y_block_raw, Y, exporters, partners, sections, channels,
           dims_full, d_0, meta,
           extra_meta = list(transform = "raw; scaling and centering are fold-specific",
                             X_shape = paste(dim(X_in_raw), collapse = "x"),
                             Y_shape = paste(dim(Y_block_raw), collapse = "x")))

channel_stats <- rbindlist(lapply(seq_len(O), function(o) {
  v <- Y[,,,o]
  data.table(channel = OUTPUT_CHANNELS[o], min = min(v), p01 = quantile(v, 0.01),
             p50 = median(v), p99 = quantile(v, 0.99), max = max(v),
             pct_zero = 100 * mean(v == 0))
}))
save_dt(channel_stats, file.path(paths$qc, "Y_channel_stats.csv"))
x_stats <- data.table(channel = INPUT_CHANNELS,
                      min = apply(X_in_raw, 2L, min),
                      p01 = apply(X_in_raw, 2L, quantile, probs = 0.01),
                      p50 = apply(X_in_raw, 2L, median),
                      p99 = apply(X_in_raw, 2L, quantile, probs = 0.99),
                      max = apply(X_in_raw, 2L, max),
                      sd = apply(X_in_raw, 2L, sd))
save_dt(x_stats, file.path(paths$qc, "X_exporter_summary_stats.csv"))

tensor_qa <- data.table(
  check = c("official_sources_only", "partner_count", "exporter_count",
            "section_count", "input_channels", "output_channels",
            "X_all_finite", "Y_all_finite", "panel_unique_keys",
            "2018_trade_weights_nonzero_each_exporter"),
  observed = c(
    as.character(needed_manifest[, all(synthetic == FALSE & exists_now & size_matches_manifest)]),
    as.character(J), as.character(n), as.character(K), as.character(F), as.character(O),
    as.character(all(is.finite(X))), as.character(all(is.finite(Y))),
    as.character(!anyDuplicated(panel[, paste(exp_iso3, imp_iso3, hs1)])),
    as.character(all(apply(trade_weights_2018, 1L, sum) > 0))
  ),
  passed = c(TRUE, J == 25L, n == 80L, K == 21L, F == 6L, O == 8L,
             all(is.finite(X)), all(is.finite(Y)),
             !anyDuplicated(panel[, paste(exp_iso3, imp_iso3, hs1)]),
             all(apply(trade_weights_2018, 1L, sum) > 0))
)
save_dt(tensor_qa, file.path(paths$qc, "tensor_validation.csv"))
if (tensor_qa[, any(!passed)]) stop("At least one tensor validation check failed")

steps <- data.table(
  step = 1:9,
  operation = c(
    "Validate official source manifest and byte sizes",
    "Map BACI numeric country identifiers with CEPII country codes",
    "Map HS2017 products to 21 HS sections from HS2 prefixes",
    "Rank and select 25 partners using 2018 BACI only",
    "Apply coverage rules and select top 80 exporters using 2018 BACI only",
    "Aggregate 2018 WITS schedules to HS sections; use lagged trade weights for H5",
    "Merge 2018 WDI, LPI, Gravity, tariffs, and lagged world demand",
    "Read 2019 BACI only after sample is frozen; use full-world denominators for shares and RCA",
    "Collapse cell-level X with 2018 trade weights; defer transforms to CV folds"
  ),
  leakage_control = c(
    "All manifest rows synthetic=FALSE", "Identifier mapping only", "Classification only",
    "No 2019 fields read", "No 2019 fields read", "Weights from 2018",
    "All covariates dated 2018", "2019 restricted to Y", "No 2019 response weights"
  )
)
save_dt(steps, file.path(paths$processed, "processing_steps.csv"),
        list(data_status = "official_real"))

log_msg(sprintf("Official preprocessing complete: X=%s, Y=%s, n=%d.",
                paste(dim(X), collapse = "x"), paste(dim(Y), collapse = "x"), n),
        stage = stage)
