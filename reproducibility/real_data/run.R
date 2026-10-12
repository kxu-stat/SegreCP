#!/usr/bin/env Rscript

# Source-checkout data preparation for the 2018–2019 trade analysis.
# All downloaded and derived data stay in --data-dir, outside this package.
args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Prepare or verify the frozen 2018–2019 trade tensor analysis inputs.",
  "Usage:",
  "  Rscript run.R --check-inputs --data-dir=/absolute/private/path",
  "  Rscript run.R --check-derived --data-dir=/absolute/private/path",
  "  Rscript run.R --check-structure --data-dir=/absolute/private/path",
  "  Rscript run.R --prepare --data-dir=/absolute/private/path",
  "    [--offline] [--overwrite-derived]",
  sep = "\n"
)
if (length(args) == 0L || "--help" %in% args) {
  cat(usage, "\n")
  quit(status = if (length(args) == 0L) 2L else 0L)
}
known <- grepl("^--data-dir=", args) |
  args %in% c("--check-inputs", "--check-derived", "--check-structure",
              "--prepare", "--offline", "--overwrite-derived")
if (any(!known)) stop("Unknown option(s): ", paste(args[!known], collapse = ", "))
mode <- intersect(args, c(
  "--check-inputs", "--check-derived", "--check-structure", "--prepare"
))
if (length(mode) != 1L) stop("Choose exactly one mode.\n", usage)
offline <- "--offline" %in% args
if (offline && mode != "--prepare") stop("--offline requires --prepare")
overwrite_derived <- "--overwrite-derived" %in% args
if (overwrite_derived && mode != "--prepare") {
  stop("--overwrite-derived requires --prepare")
}
data_arg <- args[grepl("^--data-dir=", args)]
if (length(data_arg) != 1L) stop("Specify one --data-dir.\n", usage)
data_dir <- sub("^--data-dir=", "", data_arg)
if (!nzchar(data_dir) || !startsWith(data_dir, "/")) {
  stop("--data-dir must be an absolute private path")
}

# Normalize ".." components and resolve any existing symlinked ancestor
# before creating a directory. Use the resolved path for every later write.
canonical_target <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
  stack <- character()
  for (part in parts) {
    if (!nzchar(part) || part == ".") next
    if (part == "..") {
      if (length(stack)) stack <- stack[-length(stack)]
    } else {
      stack <- c(stack, part)
    }
  }
  path <- paste0("/", paste(stack, collapse = "/"))
  cursor <- path
  suffix <- character()
  while (!file.exists(cursor)) {
    suffix <- c(basename(cursor), suffix)
    cursor <- dirname(cursor)
  }
  resolved <- normalizePath(cursor, mustWork = TRUE)
  for (part in suffix) resolved <- file.path(resolved, part)
  resolved
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE)
code_dir <- dirname(script_path)
package_dir <- normalizePath(file.path(code_dir, "..", ".."), mustWork = TRUE)
in_package <- function(path) {
  identical(path, package_dir) || startsWith(path, paste0(package_dir, "/"))
}
data_dir <- canonical_target(data_dir)
if (mode == "--prepare") {
  if (identical(data_dir, "/") || in_package(data_dir) ||
      startsWith(package_dir, paste0(data_dir, "/"))) {
    stop("--data-dir must be a dedicated directory outside the package source")
  }
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(data_dir)) stop("Data directory does not exist: ", data_dir)
data_dir <- normalizePath(data_dir, mustWork = TRUE)
if (mode == "--prepare" && in_package(data_dir)) {
  stop("--data-dir must be outside the SegreCP package source tree")
}
derived_rel <- c(
  "analysis_official/processed/tensors_real.RData",
  paste0("analysis_official/processed/preprocessing_sensitivity_inputs/",
         "v3_scope_sqrt_raw_center_X_raw.csv"),
  "analysis_official/outputs/tables/official_rank3_fold_assignments.csv",
  "analysis_official/outputs/tables/preprocessing_sensitivity_registry.csv"
)
derived_paths <- file.path(data_dir, derived_rel)
if (mode == "--prepare" && !overwrite_derived &&
    any(file.exists(derived_paths))) {
  stop("Existing derived artifact(s) would be overwritten:\n",
       paste(derived_rel[file.exists(derived_paths)], collapse = "\n"),
       "\nUse a new private --data-dir or pass --overwrite-derived explicitly.")
}
Sys.setenv(SEGRECP_REAL_DATA_DIR = data_dir)

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Install the R package 'digest' to check frozen source hashes")
}
reference <- utils::read.csv(
  file.path(code_dir, "reference_source_manifest_sha256.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
required_columns <- c(
  "local_path", "bytes", "sha256", "content_sha256", "hash_rule",
  "source_url", "provenance_status", "synthetic"
)
if (!all(required_columns %in% names(reference)) ||
    nrow(reference) != 82L ||
    anyDuplicated(reference$local_path) ||
    any(grepl("^/|(^|/)\\.\\.(/|$)", reference$local_path)) ||
    any(reference$synthetic != FALSE) ||
    any(reference$provenance_status != "official_real") ||
    any(!startsWith(reference$source_url, "https://"))) {
  stop("The bundled frozen source manifest is malformed")
}
is_tariff <- grepl("^raw_official/wits_tariffs_", reference$local_path)

# WITS adds a response timestamp to otherwise unchanged JSON. The frozen
# content hash replaces only header.prepared with a fixed token.
content_hash <- function(path, tariff) {
  if (!tariff) return(digest::digest(path, algo = "sha256", file = TRUE))
  payload <- readChar(path, file.info(path)$size, useBytes = TRUE)
  pattern <- '"prepared":"[^"]+"'
  if (!grepl(pattern, payload)) stop("No WITS header.prepared in ", path)
  normalized <- sub(pattern, '"prepared":"TIMESTAMP"', payload)
  digest::digest(normalized, algo = "sha256", serialize = FALSE)
}

check_sources <- function(rows, reject_extra_tariffs = FALSE) {
  missing <- rows$local_path[!file.exists(file.path(data_dir, rows$local_path))]
  if (length(missing)) {
    stop("Missing ", length(missing), " frozen official source file(s):\n",
         paste(utils::head(missing, 8L), collapse = "\n"))
  }
  mismatch <- character()
  for (i in seq_len(nrow(rows))) {
    rel <- rows$local_path[i]
    path <- file.path(data_dir, rel)
    tariff <- grepl("^raw_official/wits_tariffs_", rel)
    actual_hash <- content_hash(path, tariff)
    expected_hash <- if (tariff) rows$content_sha256[i] else rows$sha256[i]
    if (!identical(actual_hash, expected_hash) ||
        (!tariff && file.info(path)$size != rows$bytes[i])) {
      mismatch <- c(mismatch, rel)
    }
  }
  if (length(mismatch)) {
    stop("Frozen source content differs in ", length(mismatch),
         " file(s). The live provider may have changed its data; ",
         "the frozen trade inputs cannot be claimed to match.\n",
         paste(utils::head(mismatch, 8L), collapse = "\n"))
  }
  if (reject_extra_tariffs) {
    raw_dir <- file.path(data_dir, "raw_official")
    actual <- list.files(raw_dir, pattern = "[.]json$", recursive = TRUE)
    actual <- paste0("raw_official/", actual)
    actual <- actual[grepl("^raw_official/wits_tariffs_", actual)]
    extra <- setdiff(actual, reference$local_path[is_tariff])
    if (length(extra)) {
      stop("Unexpected WITS JSON file(s) could enter the parser:\n",
           paste(utils::head(extra, 8L), collapse = "\n"))
    }
  }
  message("Frozen official sources verified: ", nrow(rows), " files")
  invisible(TRUE)
}

check_derived <- function(strict_hash = TRUE) {
  base <- file.path(data_dir, "analysis_official")
  target <- c(
    file.path(base, "processed", "tensors_real.RData"),
    file.path(base, "processed", "preprocessing_sensitivity_inputs",
              "v3_scope_sqrt_raw_center_X_raw.csv"),
    file.path(base, "outputs", "tables", "official_rank3_fold_assignments.csv"),
    file.path(base, "outputs", "tables", "preprocessing_sensitivity_registry.csv")
  )
  missing <- target[!file.exists(target)]
  if (length(missing)) {
    stop("Missing selected trade-analysis input(s):\n",
         paste(missing, collapse = "\n"))
  }
  objects <- new.env(parent = baseenv())
  load(target[1L], envir = objects)
  needed <- c("X", "Y", "exporters")
  if (!all(needed %in% ls(objects)) ||
      !identical(as.integer(dim(objects$X)), c(80L, 25L, 21L, 6L)) ||
      !identical(as.integer(dim(objects$Y)), c(80L, 25L, 21L, 8L)) ||
      any(!is.finite(objects$X)) || any(!is.finite(objects$Y))) {
    stop("Derived tensors have unexpected dimensions or non-finite entries")
  }
  design <- utils::read.csv(target[2L], stringsAsFactors = FALSE)
  folds <- utils::read.csv(target[3L], stringsAsFactors = FALSE)
  registry <- utils::read.csv(target[4L], stringsAsFactors = FALSE)
  E <- as.character(objects$exporters$exp_iso3)
  if (nrow(design) != 80L || ncol(design) != 7L ||
      !identical(as.character(design$exp_iso3), E) ||
      any(!is.finite(as.matrix(design[-1L]))) ||
      !identical(as.character(folds$exporter), E) ||
      !identical(as.integer(table(factor(folds$fold, levels = 1:5))),
                 rep(16L, 5L))) {
    stop("Selected design or common folds fail structural checks")
  }
  selected <- registry[registry$variant_id == "v3_scope_sqrt_raw_center", ]
  if (nrow(selected) != 1L ||
      selected$rta_encoding != "scope_0_1_2" ||
      selected$aggregation_weight != "sqrt_value" ||
      selected$tariff_transform != "raw" ||
      selected$x_scaling != "mean_sd" ||
      selected$response_fit_scale != "center_only" ||
      selected$feature_count != 6L) {
    stop("The trade preprocessing registry differs from the frozen specification")
  }
  message("Selected trade inputs verified: Y=80x25x21x8, X=80x6, folds=16 each")
  frozen <- utils::read.csv(
    file.path(code_dir, "reference_derived_sha256.csv"),
    stringsAsFactors = FALSE
  )
  if (nrow(frozen) != 4L ||
      !identical(file.path(data_dir, frozen$local_path), target)) {
    stop("The bundled derived-artifact manifest is malformed")
  }
  current <- vapply(
    file.path(data_dir, frozen$local_path),
    digest::digest, character(1L), algo = "sha256", file = TRUE
  )
  exact <- current == frozen$sha256 &
    as.numeric(file.info(target)$size) == frozen$bytes
  message("Frozen derived byte hashes: ", sum(exact), "/", length(exact),
          " match")
  if (any(!exact)) {
    changed <- paste(frozen$local_path[!exact], collapse = ", ")
    if (strict_hash) {
      stop("Derived trade-analysis byte hashes differ from the frozen results: ", changed,
           ". Use --check-structure to inspect dimensions without claiming ",
           "exact artifact reproduction.")
    }
    message("Changed byte hashes: ", changed)
  }
  invisible(TRUE)
}

run_stage <- function(program, script) {
  message("Running ", basename(script))
  status <- system2(program, shQuote(script), stdout = "", stderr = "")
  if (!identical(as.integer(status), 0L)) {
    stop(basename(script), " failed with exit status ", status)
  }
}

if (mode %in% c("--check-inputs", "--check-derived", "--check-structure")) {
  if (mode != "--check-derived") {
    check_sources(reference, reject_extra_tariffs = TRUE)
  }
  check_derived(strict_hash = mode != "--check-structure")
  quit(status = 0L)
}

needed_packages <- c("data.table", "digest", "jsonlite", "xml2")
missing_packages <- needed_packages[
  !vapply(needed_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript is required")
if (!offline) {
  for (program in c("bash", "curl", "unzip")) {
    if (!nzchar(Sys.which(program))) stop(program, " is required for downloads")
  }
  run_stage("bash", file.path(code_dir, "scripts",
                              "00_download_core_official_sources.sh"))
  check_sources(reference[!is_tariff, , drop = FALSE])
  run_stage(rscript, file.path(code_dir, "scripts",
                               "01_select_wits_reporters.R"))
  run_stage(rscript, file.path(code_dir, "scripts",
                               "04_prepare_wits_2018_queries.R"))
  for (year in c("2018", "2019")) {
    for (datatype in c("reported", "aveestimated")) {
      Sys.setenv(WITS_YEAR = year, WITS_DATATYPE = datatype)
      run_stage("bash", file.path(code_dir, "scripts",
                                  "02_download_wits_tariffs.sh"))
    }
  }
}
# Never derive tensors from a changed source, even when a download completed.
check_sources(reference, reject_extra_tariffs = TRUE)
if (offline) {
  # Raw-only offline directories need the same deterministic query metadata as
  # online runs; it is rebuilt locally before the WITS responses are parsed.
  run_stage(rscript, file.path(code_dir, "scripts",
                               "01_select_wits_reporters.R"))
  run_stage(rscript, file.path(code_dir, "scripts",
                               "04_prepare_wits_2018_queries.R"))
}
for (year in c("2018", "2019")) {
  Sys.setenv(WITS_YEAR = year)
  run_stage(rscript, file.path(code_dir, "scripts",
                               "03_parse_and_validate_wits_tariffs.R"))
}
run_stage(rscript, file.path(code_dir, "scripts",
                             "05_write_source_manifest.R"))
run_stage(rscript, file.path(code_dir, "scripts",
                             "06_build_official_covariate_tables.R"))
run_stage(rscript, file.path(code_dir, "scripts", "analysis",
                             "01_preprocess_and_build_tensors.R"))
run_stage(rscript, file.path(code_dir, "scripts", "analysis",
                             "02_build_selected_design.R"))
check_derived()
message("Data preparation complete. Restricted data remain under ", data_dir)
