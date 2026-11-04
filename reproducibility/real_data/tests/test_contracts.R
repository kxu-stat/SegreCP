#!/usr/bin/env Rscript

# Fast, data-free checks for the repository-only trade workflow.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this file with Rscript, not source().", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
workflow_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
package_dir <- normalizePath(file.path(workflow_dir, "..", ".."), mustWork = TRUE)

check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
sha256 <- function(x) {
  is.character(x) && all(!is.na(x) & grepl("^[0-9a-f]{64}$", x))
}
relative_paths <- function(x) {
  is.character(x) && all(!is.na(x) & nzchar(x)) &&
    !any(grepl("^/|^[A-Za-z]:|\\\\|(^|/)\\.\\.(/|$)", x)) &&
    !anyDuplicated(x)
}

r_files <- list.files(workflow_dir, pattern = "[.]R$", recursive = TRUE,
                      full.names = TRUE)
check(length(r_files) >= 10L, "Expected preprocessing and fitting R scripts")
for (path in r_files) {
  tryCatch(parse(file = path), error = function(e) {
    stop("R syntax error in ", basename(path), ": ", conditionMessage(e),
         call. = FALSE)
  })
}

reader_env <- new.env(parent = baseenv())
base::source(file.path(workflow_dir, "read_prepared_trade_data.R"),
             local = reader_env)
check(is.function(reader_env$read_prepared_trade_data) &&
        identical(names(formals(reader_env$read_prepared_trade_data)),
                  "data_dir"),
      "The prepared-data reader must expose a one-argument source-only API")
missing_dir <- tempfile(pattern = "segrecp-missing-trade-")
check(!dir.exists(missing_dir), "The proposed missing-data test path exists")
missing_error <- tryCatch({
  reader_env$read_prepared_trade_data(missing_dir)
  ""
}, error = conditionMessage)
check(grepl("existing private data directory", missing_error, fixed = TRUE),
      "The prepared-data reader must explain missing private inputs")

source <- utils::read.csv(
  file.path(workflow_dir, "reference_source_manifest_sha256.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
required_source <- c("local_path", "bytes", "sha256", "content_sha256",
                     "hash_rule", "source_url", "provenance_status", "synthetic")
check(nrow(source) == 82L && all(required_source %in% names(source)),
      "The frozen source manifest must have 82 records and its required columns")
check(relative_paths(source$local_path),
      "Source paths must be unique, relative, and free of parent traversal")
check(is.numeric(source$bytes) && all(is.finite(source$bytes) & source$bytes > 0),
      "Source sizes must be positive finite numbers")
check(sha256(source$sha256) && sha256(source$content_sha256),
      "Source SHA-256 fields must contain lowercase 64-digit hashes")
check(all(!is.na(source$source_url) & startsWith(source$source_url, "https://")),
      "Each source must have an HTTPS provider URL")
check(all(source$provenance_status == "official_real") &&
        is.logical(source$synthetic) && all(!is.na(source$synthetic) &
                                              !source$synthetic),
      "The manifest must identify only official, nonsynthetic sources")
tariff <- grepl("^raw_official/wits_tariffs_", source$local_path)
check(all(grepl("header[.]prepared", source$hash_rule[tariff])) &&
        all(source$hash_rule[!tariff] == "exact bytes"),
      "Only WITS tariff responses may use the normalized timestamp hash")

derived <- utils::read.csv(
  file.path(workflow_dir, "reference_derived_sha256.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
check(identical(names(derived), c("local_path", "bytes", "sha256")) &&
        nrow(derived) == 4L,
      "The frozen derived manifest must describe exactly four artifacts")
check(relative_paths(derived$local_path) &&
        all(startsWith(derived$local_path, "analysis_official/")),
      "Derived paths must be unique and relative to the private analysis tree")
check(is.numeric(derived$bytes) && all(is.finite(derived$bytes) &
                                         derived$bytes > 0) && sha256(derived$sha256),
      "Derived sizes and SHA-256 hashes must be valid")

rscript <- Sys.which("Rscript")
check(nzchar(rscript), "Rscript is required for CLI checks")
run_cli <- function(script, args) {
  output <- suppressWarnings(system2(
    rscript, c(shQuote(script), vapply(args, shQuote, character(1L))),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  list(status = if (is.null(status)) 0L else as.integer(status),
       output = paste(output, collapse = "\n"))
}

driver <- file.path(workflow_dir, "run.R")
help <- run_cli(driver, "--help")
check(help$status == 0L && grepl("--prepare", help$output, fixed = TRUE) &&
        grepl("--check-inputs", help$output, fixed = TRUE),
      "run.R --help must succeed and describe the preparation/check modes")
fit_driver <- file.path(workflow_dir, "fit_authors_methods.R")
fit_help <- run_cli(fit_driver, "--help")
check(fit_help$status == 0L &&
        grepl("rank-three RGD/RGN", fit_help$output, fixed = TRUE) &&
        grepl("--output-dir", fit_help$output, fixed = TRUE) &&
        grepl("--folds", fit_help$output, fixed = TRUE) &&
        grepl("--smoke", fit_help$output, fixed = TRUE),
      "The authors-method CLI help must describe its purpose and fit options")

# Evaluate the runner's actual control and initialization expressions on a
# small numeric fixture; no official data or solver fit is needed.
fit_expressions <- parse(file = fit_driver)
assignment_value <- function(name) {
  values <- list()
  visit <- function(x) {
    if (missing(x)) return(invisible(NULL))
    if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
        identical(x[[2L]], as.name(name))) {
      values[[length(values) + 1L]] <<- x[[3L]]
    }
    if (is.call(x) || is.expression(x)) lapply(as.list(x), visit)
    invisible(NULL)
  }
  visit(fit_expressions)
  check(length(values) == 1L, paste("Expected one runner assignment for", name))
  values[[1L]]
}
check(requireNamespace("SegreCP", quietly = TRUE),
      "Install SegreCP before checking the fitting protocol")
protocol_env <- new.env(parent = globalenv())
protocol_env$iterations <- 30L
protocol_env$smoke <- FALSE
control <- eval(assignment_value("control"), protocol_env)
check(identical(control$ridge, 0) &&
        identical(control$normalize_regression, "none") &&
        identical(control$component_order, "decreasing_weight") &&
        identical(control$step_rule, "backtracking") && control$tolerance == 0,
      "The fresh trade fit must retain zero ridge and the original solver scale")
check(identical(eval(assignment_value("checkpoints"), protocol_env),
                c(3L, 5L, 10L, 30L)),
      "The full trade run must retain iterations 3, 5, 10 and 30")
set.seed(193L)
protocol_env$Xtr <- matrix(rnorm(12L * 6L), 12L, 6L)
protocol_env$Ytr_matrix <- matrix(rnorm(12L * 25L * 21L * 8L), 12L)
protocol_env$moment <- eval(assignment_value("moment"), protocol_env)
expected_moment <- array(crossprod(protocol_env$Xtr, protocol_env$Ytr_matrix) /
                          nrow(protocol_env$Xtr), c(6L, 25L, 21L, 8L))
check(isTRUE(all.equal(protocol_env$moment, expected_moment, tolerance = 1e-12)),
      "CPCA must use the training mean cross-moment, not the unnormalized sum")
initial <- eval(assignment_value("initial"), protocol_env)
expected_initial <- SegreCP::cp_initialize(expected_moment, rank = 3L,
  method = "cpca", split = c(1L, 2L), seed = 2549L)
check(isTRUE(all.equal(SegreCP::reconstruct_cp(initial$comps),
                       SegreCP::reconstruct_cp(expected_initial$comps),
                       tolerance = 1e-10)),
      "The shared trade initializer must preserve the normalized CPCA configuration")

# These requests must fail before they create a package-local directory.
guard_dir <- tempfile(pattern = ".segrecp-write-guard-", tmpdir = package_dir)
check(!file.exists(guard_dir), "The proposed write-guard path already exists")
prepare <- run_cli(driver, c("--prepare", "--offline",
                             paste0("--data-dir=", guard_dir)))
check(prepare$status != 0L && grepl("outside the package", prepare$output) &&
        !file.exists(guard_dir),
      "run.R must refuse a private-data directory inside the package source")

fit <- run_cli(fit_driver,
               c(paste0("--data-dir=", tempdir()),
                 paste0("--output-dir=", guard_dir)))
check(fit$status != 0L && grepl("outside the SegreCP package", fit$output) &&
        !file.exists(guard_dir),
      "The fitting driver must refuse outputs inside the package source")

cat(sprintf("PASS: parsed %d R scripts; validated 82 source and 4 derived ",
            length(r_files)))
cat("manifest rows; checked CLI help, zero-ridge controls, normalized CPCA, ",
    "checkpoints and package write guards.\n")
