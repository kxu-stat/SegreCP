#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

data_dir <- Sys.getenv("SEGRECP_REAL_DATA_DIR", unset = "")
if (!nzchar(data_dir)) stop("SEGRECP_REAL_DATA_DIR must name a private data directory")
root <- normalizePath(data_dir, mustWork = TRUE)

relative_path <- function(path) {
  resolved <- normalizePath(path, mustWork = TRUE)
  prefix <- paste0(root, "/")
  if (!startsWith(resolved, prefix)) stop("Source is outside the private data directory")
  substring(resolved, nchar(prefix) + 1L)
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE)
}

# WITS embeds the response-generation time in `header.prepared`.  That field
# changes on an otherwise identical download, so preserve the byte-level hash
# and also compute a stable content hash after replacing only this timestamp.
wits_content_identity <- function(path) {
  payload <- readChar(path, file.info(path)$size, useBytes = TRUE)
  prepared_pattern <- '"prepared":"[^"]+"'
  prepared_match <- regmatches(payload, regexpr(prepared_pattern, payload))
  prepared_at <- if (length(prepared_match) && nzchar(prepared_match)) {
    sub('^"prepared":"|"$', "", prepared_match)
  } else {
    NA_character_
  }
  normalized <- sub(prepared_pattern, '"prepared":"TIMESTAMP"', payload)
  list(
    content_sha256 = digest(normalized, algo = "sha256", serialize = FALSE),
    hash_rule = "WITS JSON with header.prepared replaced by TIMESTAMP",
    response_prepared_at = prepared_at
  )
}

core_files <- c(
  file.path(root, "source_archives", "BACI_HS17_V202601.zip"),
  file.path(root, "source_archives", "Gravity_rds_V202211.zip"),
  file.path(root, "raw_official", "BACI_HS17_Y2018_V202601.csv"),
  file.path(root, "raw_official", "BACI_HS17_Y2019_V202601.csv"),
  file.path(root, "raw_official", "country_codes_V202601.csv"),
  file.path(root, "raw_official", "product_codes_HS17_V202601.csv"),
  file.path(root, "raw_official", "Gravity_V202211.rds"),
  file.path(root, "raw_official", "wdi_gdp_current_usd_2018_2019.json"),
  file.path(root, "raw_official", "wdi_gdp_pc_current_usd_2018_2019.json"),
  file.path(root, "raw_official", "wdi_population_2018_2019.json"),
  file.path(root, "raw_official", "world_bank_lpi_overall_2018.json"),
  file.path(root, "raw_official", "wits_trains_country_metadata.xml"),
  file.path(root, "raw_official", "wits_trains_dataavailability_2018.xml"),
  file.path(root, "raw_official", "wits_trains_dataavailability_2019.xml")
)
wits_files <- list.files(
  file.path(root, "raw_official"),
  pattern = "[.]json$",
  recursive = TRUE,
  full.names = TRUE
)
wits_files <- wits_files[grepl("wits_tariffs_20(18|19)", wits_files)]
files <- sort(unique(c(core_files, wits_files)))
if (!all(file.exists(files))) stop("A manifest input file is missing")

source_info <- function(path) {
  rel <- relative_path(path)
  base <- basename(path)
  if (grepl("BACI", base) || grepl("V202601", base)) {
    return(list(
      source_organization = "CEPII",
      dataset = "BACI HS2017",
      source_url = "https://www.cepii.fr/DATA_DOWNLOAD/baci/data/BACI_HS17_V202601.zip"
    ))
  }
  if (grepl("Gravity|Countries_V202211|Label_", base)) {
    return(list(
      source_organization = "CEPII",
      dataset = "Gravity",
      source_url = "https://www.cepii.fr/DATA_DOWNLOAD/gravity/data/Gravity_rds_V202211.zip"
    ))
  }
  if (grepl("wdi_gdp_current", base)) {
    return(list(
      source_organization = "World Bank",
      dataset = "WDI NY.GDP.MKTP.CD",
      source_url = "https://api.worldbank.org/v2/country/all/indicator/NY.GDP.MKTP.CD?date=2018:2019&format=json&per_page=20000"
    ))
  }
  if (grepl("wdi_gdp_pc", base)) {
    return(list(
      source_organization = "World Bank",
      dataset = "WDI NY.GDP.PCAP.CD",
      source_url = "https://api.worldbank.org/v2/country/all/indicator/NY.GDP.PCAP.CD?date=2018:2019&format=json&per_page=20000"
    ))
  }
  if (grepl("wdi_population", base)) {
    return(list(
      source_organization = "World Bank",
      dataset = "WDI SP.POP.TOTL",
      source_url = "https://api.worldbank.org/v2/country/all/indicator/SP.POP.TOTL?date=2018:2019&format=json&per_page=20000"
    ))
  }
  if (grepl("lpi", base, ignore.case = TRUE)) {
    return(list(
      source_organization = "World Bank",
      dataset = "LPI LP.LPI.OVRL.XQ",
      source_url = "https://api.worldbank.org/v2/country/all/indicator/LP.LPI.OVRL.XQ?date=2018&format=json&per_page=20000"
    ))
  }
  if (grepl("wits_trains_country_metadata", base)) {
    return(list(
      source_organization = "World Bank WITS / UNCTAD TRAINS",
      dataset = "TRAINS country metadata",
      source_url = "https://wits.worldbank.org/API/V1/wits/datasource/trn/country/ALL"
    ))
  }
  if (grepl("wits_trains_dataavailability", base)) {
    year <- sub(".*_(20[0-9]{2})[.]xml$", "\\1", base)
    return(list(
      source_organization = "World Bank WITS / UNCTAD TRAINS",
      dataset = paste("TRAINS data availability", year),
      source_url = sprintf(
        "https://wits.worldbank.org/API/V1/wits/datasource/trn/dataavailability/country/all/year/%s",
        year
      )
    ))
  }
  if (grepl("wits_tariffs_", rel)) {
    year <- sub(".*_(20[0-9]{2})[.]json$", "\\1", base)
    code <- sub("^[A-Z]{3}_([0-9]{3})_.*$", "\\1", base)
    datatype <- if (grepl("aveestimated", base)) "aveestimated" else "reported"
    return(list(
      source_organization = "World Bank WITS / UNCTAD TRAINS",
      dataset = paste("TRAINS MFN tariff", datatype, year),
      source_url = sprintf(
        paste0(
          "https://wits.worldbank.org/API/V1/SDMX/V21/datasource/TRN/",
          "reporter/%s/partner/000/product/all/year/%s/datatype/%s?format=JSON"
        ),
        code, year, datatype
      )
    ))
  }
  stop("No source mapping for ", rel)
}

message("Computing SHA-256 for ", length(files), " official source files...")
manifest <- rbindlist(lapply(files, function(path) {
  info <- source_info(path)
  rel <- relative_path(path)
  byte_sha <- sha256_file(path)
  identity <- if (grepl("raw_official/wits_tariffs_", rel)) {
    wits_content_identity(path)
  } else {
    list(content_sha256 = byte_sha, hash_rule = "exact bytes",
         response_prepared_at = NA_character_)
  }
  data.table(
    local_path = rel,
    bytes = file.info(path)$size,
    sha256 = byte_sha,
    content_sha256 = identity$content_sha256,
    hash_rule = identity$hash_rule,
    response_prepared_at = identity$response_prepared_at,
    source_organization = info$source_organization,
    dataset = info$dataset,
    source_url = info$source_url,
    provenance_status = "official_real",
    synthetic = FALSE
  )
}))

fwrite(manifest, file.path(root, "metadata", "source_manifest_sha256.csv"))
message("Wrote ", nrow(manifest), " manifest rows; all are marked synthetic=FALSE.")
