#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

data_dir <- Sys.getenv("SEGRECP_REAL_DATA_DIR", unset = "")
if (!nzchar(data_dir)) stop("SEGRECP_REAL_DATA_DIR must name a private data directory")
root <- normalizePath(data_dir, mustWork = TRUE)
year <- as.integer(Sys.getenv("WITS_YEAR", "2019"))
if (is.na(year) || year < 1900L || year > 2100L) stop("Invalid WITS_YEAR")
raw_reported_dir <- file.path(
  root,
  "raw_official",
  sprintf("wits_tariffs_%d_by_reporter", year)
)
raw_aveestimated_dir <- file.path(
  root,
  "raw_official",
  sprintf("wits_tariffs_%d_aveestimated_by_reporter", year)
)
processed_dir <- file.path(root, "processed_official")
meta_dir <- file.path(root, "metadata")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

query_map <- fread(
  file.path(meta_dir, sprintf("wits_unique_queries_%d.csv", year)),
  colClasses = list(character = c("wits_reporter_code", "wits_source_iso3"))
)

attribute_value <- function(observation, attribute_position, attribute_definition) {
  index <- observation[[attribute_position + 1L]]
  if (is.null(index) || length(index) == 0L || is.na(index)) {
    return(NA_character_)
  }
  values <- attribute_definition$values
  if (length(values) < index + 1L) {
    stop("WITS observation attribute index is out of bounds")
  }
  values[[index + 1L]]$id
}

parse_one <- function(path) {
  payload <- fromJSON(path, simplifyVector = FALSE)
  series_definitions <- payload$structure$dimensions$series
  series_ids <- vapply(series_definitions, `[[`, character(1L), "id")
  reporter_position <- match("REPORTER", series_ids)
  partner_position <- match("PARTNER", series_ids)
  product_position <- match("PRODUCTCODE", series_ids)
  datatype_position <- match("DATATYPE", series_ids)

  if (anyNA(c(reporter_position, partner_position, product_position, datatype_position))) {
    stop("Required WITS series dimensions are missing in ", basename(path))
  }

  time_definition <- payload$structure$dimensions$observation[[1L]]
  year <- as.integer(time_definition$values[[1L]]$id)
  attribute_definitions <- payload$structure$attributes$observation
  attribute_ids <- vapply(attribute_definitions, `[[`, character(1L), "id")
  series <- payload$dataSets[[1L]]$series
  key_matrix <- do.call(
    rbind,
    lapply(strsplit(names(series), ":", fixed = TRUE), as.integer)
  )
  varying_key_positions <- which(apply(key_matrix, 2L, max) > 0L)
  if (length(varying_key_positions) != 1L) {
    stop("Expected exactly one varying WITS series-key component in ", basename(path))
  }
  product_key_position <- varying_key_positions[[1L]]

  # WITS currently serializes PRODUCTCODE in the third key component even
  # though its metadata advertises keyPosition=3 (zero-based fourth). The
  # response contains one reporter, partner, and datatype per file, so use
  # their sole dimension values and identify the product component by variation.
  sole_series_value <- function(position) {
    values <- series_definitions[[position]]$values
    if (length(values) != 1L) {
      stop("Expected a singleton WITS query dimension in ", basename(path))
    }
    values[[1L]]$id
  }
  reporter_value <- sole_series_value(reporter_position)
  partner_value <- sole_series_value(partner_position)
  datatype_value <- sole_series_value(datatype_position)
  product_values <- series_definitions[[product_position]]$values

  rows <- lapply(seq_along(series), function(row_index) {
    product_index <- key_matrix[row_index, product_key_position]
    observation <- series[[row_index]]$observations[[1L]]
    obs_attribute <- function(id) {
      position <- match(id, attribute_ids)
      if (is.na(position)) NA_character_ else
        attribute_value(observation, position, attribute_definitions[[position]])
    }

    data.table(
      wits_reporter_code = reporter_value,
      partner_code = partner_value,
      hs6 = product_values[[product_index + 1L]]$id,
      datatype = datatype_value,
      year = year,
      tariff_simple_average_pct = as.numeric(observation[[1L]]),
      nomenclature = obs_attribute("NOMENCODE"),
      tariff_type = obs_attribute("TARIFFTYPE"),
      sum_of_rates = as.numeric(obs_attribute("SUM_OF_RATES")),
      min_rate = as.numeric(obs_attribute("MIN_RATE")),
      max_rate = as.numeric(obs_attribute("MAX_RATE")),
      total_lines = as.integer(obs_attribute("TOTALNOOFLINES")),
      n_preferential_lines = as.integer(obs_attribute("NBR_PREF_LINES")),
      n_mfn_lines = as.integer(obs_attribute("NBR_MFN_LINES")),
      n_na_lines = as.integer(obs_attribute("NBR_NA_LINES")),
      observation_measure = obs_attribute("OBS_VALUE_MEASURE"),
      source_datatype_request = tolower(datatype_value),
      source_file = basename(path)
    )
  })
  rbindlist(rows)
}

reported_files <- sort(list.files(
  raw_reported_dir,
  pattern = sprintf("_%d[.]json$", year),
  full.names = TRUE
))
aveestimated_files <- sort(list.files(
  raw_aveestimated_dir,
  pattern = sprintf("_aveestimated_%d[.]json$", year),
  full.names = TRUE
))
expected_aveestimated <- query_map[wits_aveestimated_available == "Yes", .N]
if (length(reported_files) != nrow(query_map)) {
  stop("Expected ", nrow(query_map), " reported WITS files but found ", length(reported_files))
}
if (length(aveestimated_files) != expected_aveestimated) {
  stop(
    "Expected ", expected_aveestimated,
    " AVE-estimated WITS files but found ", length(aveestimated_files)
  )
}

message(
  "Parsing ", length(reported_files), " reported and ",
  length(aveestimated_files), " AVE-estimated official WITS schedules..."
)
reported_tariffs <- rbindlist(lapply(reported_files, parse_one), use.names = TRUE)
aveestimated_tariffs <- rbindlist(
  lapply(aveestimated_files, parse_one),
  use.names = TRUE
)
expected_ave_codes <- query_map[
  wits_aveestimated_available == "Yes",
  unique(wits_reporter_code)
]
if (!setequal(unique(aveestimated_tariffs$wits_reporter_code), expected_ave_codes)) {
  stop("AVE-estimated file coverage disagrees with WITS availability metadata")
}

# Prefer WITS/UNCTAD's AVE-estimated schedule when the official availability
# metadata says that specific-duty conversion exists. Otherwise retain the
# reported schedule; no locally imputed tariff values are introduced.
tariffs <- rbindlist(list(
  aveestimated_tariffs,
  reported_tariffs[!wits_reporter_code %chin% expected_ave_codes]
))
tariffs <- merge(
  tariffs,
  query_map[, .(
    wits_reporter_code,
    wits_source_iso3,
    expected_nomenclature = wits_nomenclature
  )],
  by = "wits_reporter_code",
  all.x = TRUE
)
setcolorder(
  tariffs,
  c("wits_source_iso3", "wits_reporter_code", "partner_code", "year", "hs6")
)
setorder(tariffs, wits_source_iso3, hs6)
requested_year <- year

qa <- tariffs[, .(
  n_rows = .N,
  n_unique_hs6 = uniqueN(hs6),
  n_duplicate_hs6 = .N - uniqueN(hs6),
  n_missing_tariff = sum(!is.finite(tariff_simple_average_pct)),
  share_missing_tariff = mean(!is.finite(tariff_simple_average_pct)),
  n_zero_tariff = sum(tariff_simple_average_pct == 0, na.rm = TRUE),
  share_zero_tariff = mean(tariff_simple_average_pct == 0, na.rm = TRUE),
  min_tariff_pct = min(tariff_simple_average_pct, na.rm = TRUE),
  median_tariff_pct = median(tariff_simple_average_pct, na.rm = TRUE),
  max_tariff_pct = max(tariff_simple_average_pct, na.rm = TRUE),
  n_negative_tariff = sum(tariff_simple_average_pct < 0, na.rm = TRUE),
  n_unexpected_nomenclature = sum(
    nomenclature != expected_nomenclature | is.na(nomenclature)
  ),
  n_non_mfn = sum(tariff_type != "MFN" | is.na(tariff_type)),
  n_non_requested_year = sum(year != requested_year | is.na(year)),
  n_non_world_partner = sum(partner_code != "000" | is.na(partner_code)),
  source_datatype_request = unique(source_datatype_request)
), by = .(wits_source_iso3, wits_reporter_code)]

summary_by_request <- function(dt, label) {
  dt[, .(
    n_rows = .N,
    n_zero_tariff = sum(tariff_simple_average_pct == 0, na.rm = TRUE),
    share_zero_tariff = mean(tariff_simple_average_pct == 0, na.rm = TRUE),
    median_tariff_pct = median(tariff_simple_average_pct, na.rm = TRUE),
    max_tariff_pct = max(tariff_simple_average_pct, na.rm = TRUE),
    request = label
  ), by = wits_reporter_code]
}
reported_vs_preferred <- merge(
  dcast(
    rbindlist(list(
      summary_by_request(reported_tariffs, "reported"),
      summary_by_request(tariffs, "preferred")
    )),
    wits_reporter_code ~ request,
    value.var = c(
      "n_rows", "n_zero_tariff", "share_zero_tariff",
      "median_tariff_pct", "max_tariff_pct"
    )
  ),
  query_map[, .(
    wits_reporter_code,
    wits_source_iso3,
    wits_aveestimated_available
  )],
  by = "wits_reporter_code",
  all.x = TRUE
)

if (anyNA(tariffs$wits_source_iso3)) stop("Unmapped WITS reporter codes remain")
if (tariffs[, anyDuplicated(paste(wits_source_iso3, hs6))] != 0L) {
  stop("Duplicate reporter-HS6 keys detected")
}
if (qa[, any(n_unexpected_nomenclature > 0L | n_non_mfn > 0L |
             n_non_requested_year > 0L |
             n_non_world_partner > 0L | n_negative_tariff > 0L)]) {
  stop("WITS validity check failed; inspect QA output")
}

fwrite(
  tariffs,
  file.path(processed_dir, sprintf("wits_mfn_tariffs_preferred_%d.csv", year))
)
fwrite(qa, file.path(meta_dir, sprintf("wits_tariff_qa_%d.csv", year)))
fwrite(
  reported_vs_preferred,
  file.path(
    meta_dir,
    sprintf("wits_tariff_reported_vs_preferred_%d.csv", year)
  )
)

message("Built ", nrow(tariffs), " preferred reporter-product rows.")
message(
  "All rows are ", year,
  ", MFN, world-partner observations with WITS-certified nomenclature."
)
message("No local tariff imputation was used.")
print(qa)
