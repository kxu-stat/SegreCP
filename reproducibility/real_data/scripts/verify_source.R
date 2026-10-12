#!/usr/bin/env Rscript

# Check a staged official download against the frozen trade-source manifest.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript verify_source.R PATH RELATIVE_SOURCE_PATH")
}
path <- args[1L]
relative_path <- args[2L]
if (!file.exists(path) || !requireNamespace("digest", quietly = TRUE)) {
  stop("Staged source is missing or R package 'digest' is unavailable")
}
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE)
reference <- utils::read.csv(
  file.path(dirname(script_path), "..", "reference_source_manifest_sha256.csv"),
  stringsAsFactors = FALSE
)
row <- reference[reference$local_path == relative_path, , drop = FALSE]
if (nrow(row) != 1L) stop("No frozen source reference for ", relative_path)
is_tariff <- grepl("^raw_official/wits_tariffs_", relative_path)
if (is_tariff) {
  payload <- readChar(path, file.info(path)$size, useBytes = TRUE)
  pattern <- '"prepared":"[^"]+"'
  if (!grepl(pattern, payload)) stop("WITS header.prepared is missing")
  normalized <- sub(pattern, '"prepared":"TIMESTAMP"', payload)
  actual <- digest::digest(normalized, algo = "sha256", serialize = FALSE)
  expected <- row$content_sha256
} else {
  actual <- digest::digest(path, algo = "sha256", file = TRUE)
  expected <- row$sha256
  if (file.info(path)$size != row$bytes) stop("Frozen source size differs")
}
if (!identical(actual, expected)) {
  stop("Frozen source hash differs for ", relative_path,
       "; keep the existing file and inspect provider version")
}
