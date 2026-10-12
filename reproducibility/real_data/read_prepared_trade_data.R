# Read the frozen paper inputs from a private data directory. This source-only
# helper is not part of the general SegreCP estimator API and distributes no data.
.trade_reader_dir <- local({
  files <- unlist(lapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
  }), use.names = FALSE)
  files <- files[!is.na(files) & nzchar(files)]
  if (!length(files)) {
    stop("Source read_prepared_trade_data.R from the SegreCP checkout.",
         call. = FALSE)
  }
  dirname(normalizePath(utils::tail(files, 1L), mustWork = TRUE))
})

#' Read the prepared 2018–2019 trade analysis inputs
#'
#' Supply the private directory created by `run.R --prepare`, not the
#' `analysis_official` subdirectory. The function verifies all four frozen
#' processed-file SHA-256 hashes and sizes before loading the data. It then
#' checks the registered preprocessing choice, tensor shape, exporter keys,
#' feature order, and balanced fold assignments. It does not download data,
#' fit an estimator, or verify the 82 official raw source files.
#'
#' @param data_dir Existing private data directory passed to `run.R`.
#' @return A list with `X` (80-by-6 numeric design matrix), `Y` (80-by-25-by-21-
#'   by-8 response array), `fold` (five aligned fold IDs), `exporters` (aligned
#'   ISO3 keys), `feature_names`, `metadata`, and `input_hashes`.
read_prepared_trade_data <- function(data_dir) {
  if (!is.character(data_dir) || length(data_dir) != 1L ||
      is.na(data_dir) || !nzchar(data_dir) || !dir.exists(data_dir)) {
    stop("`data_dir` must be an existing private data directory.", call. = FALSE)
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Install the 'digest' package to verify the prepared inputs.",
         call. = FALSE)
  }
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  relative <- c(
    tensor = "analysis_official/processed/tensors_real.RData",
    design = paste0("analysis_official/processed/",
                    "preprocessing_sensitivity_inputs/",
                    "v3_scope_sqrt_raw_center_X_raw.csv"),
    folds = "analysis_official/outputs/tables/official_rank3_fold_assignments.csv",
    registry = "analysis_official/outputs/tables/preprocessing_sensitivity_registry.csv"
  )
  paths <- stats::setNames(file.path(data_dir, relative), names(relative))
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing prepared trade input(s): ", paste(missing, collapse = ", "),
         ". Run the data preparation workflow first.", call. = FALSE)
  }

  manifest <- file.path(.trade_reader_dir, "reference_derived_sha256.csv")
  if (!file.exists(manifest)) {
    stop("The frozen processed-input manifest is missing from the checkout.",
         call. = FALSE)
  }
  expected <- utils::read.csv(manifest, stringsAsFactors = FALSE,
                              check.names = FALSE)
  if (!identical(names(expected), c("local_path", "bytes", "sha256")) ||
      nrow(expected) != 4L || anyDuplicated(expected$local_path) ||
      !setequal(expected$local_path, unname(relative))) {
    stop("The frozen processed-input manifest has an unexpected schema or path set.",
         call. = FALSE)
  }
  expected <- expected[match(relative, expected$local_path), , drop = FALSE]
  actual_size <- as.numeric(file.info(paths)$size)
  actual_hash <- vapply(paths, digest::digest, character(1L),
                        algo = "sha256", file = TRUE)
  if (anyNA(actual_size) || any(actual_size != expected$bytes) ||
      any(tolower(expected$sha256) != unname(actual_hash))) {
    stop("Prepared trade inputs differ from the frozen paper data. ",
         "Run `run.R --check-structure` to diagnose a changed source; ",
         "do not describe these inputs as an exact paper replay.",
         call. = FALSE)
  }

  tensor_env <- new.env(parent = baseenv())
  load(paths[["tensor"]], envir = tensor_env)
  if (!all(vapply(c("Y", "exporters", "meta"), exists, logical(1L),
                  envir = tensor_env, inherits = FALSE))) {
    stop("The response file lacks Y, exporters, or metadata.", call. = FALSE)
  }
  Y <- tensor_env$Y
  exporter_table <- tensor_env$exporters
  metadata <- tensor_env$meta
  if (!is.numeric(Y) || !identical(dim(Y), c(80L, 25L, 21L, 8L)) ||
      any(!is.finite(Y))) {
    stop("The response must be a finite 80 x 25 x 21 x 8 array.",
         call. = FALSE)
  }
  if (!is.data.frame(exporter_table) ||
      !"exp_iso3" %in% names(exporter_table)) {
    stop("The response file lacks exporter ISO3 keys.", call. = FALSE)
  }
  exporters <- as.character(exporter_table$exp_iso3)
  if (length(exporters) != 80L || anyNA(exporters) ||
      any(!nzchar(exporters)) || anyDuplicated(exporters)) {
    stop("Exporter ISO3 keys must contain 80 distinct nonempty values.",
         call. = FALSE)
  }
  if (!is.list(metadata) || !identical(metadata$data_status, "official_real") ||
      !identical(as.integer(metadata$selection_year), 2018L) ||
      !identical(as.integer(metadata$covariate_year), 2018L) ||
      !identical(as.integer(metadata$outcome_year), 2019L)) {
    stop("Response metadata do not describe the frozen 2018–2019 analysis.",
         call. = FALSE)
  }

  feature_names <- c("tariff", "ntb_proxy", "log_gdp_partner",
                     "log_distance", "rta_scope_count", "world_demand")
  design <- utils::read.csv(paths[["design"]], check.names = FALSE)
  if (!identical(names(design), c("exp_iso3", feature_names)) ||
      nrow(design) != 80L || anyNA(design$exp_iso3) ||
      anyDuplicated(design$exp_iso3) ||
      !setequal(as.character(design$exp_iso3), exporters) ||
      !all(vapply(design[feature_names], is.numeric, logical(1L)))) {
    stop("The design must have six numeric features and one row per exporter.",
         call. = FALSE)
  }
  design <- design[match(exporters, design$exp_iso3), , drop = FALSE]
  X <- as.matrix(design[feature_names])
  storage.mode(X) <- "double"
  if (any(!is.finite(X))) {
    stop("The design contains nonfinite values.", call. = FALSE)
  }

  fold_table <- utils::read.csv(paths[["folds"]], check.names = FALSE)
  if (!all(c("exporter", "fold") %in% names(fold_table)) ||
      nrow(fold_table) != 80L || anyNA(fold_table$exporter) ||
      anyDuplicated(fold_table$exporter) ||
      !setequal(as.character(fold_table$exporter), exporters) ||
      !is.numeric(fold_table$fold) || anyNA(fold_table$fold) ||
      any(fold_table$fold != as.integer(fold_table$fold)) ||
      !identical(as.integer(table(factor(fold_table$fold, levels = 1:5))),
                 rep(16L, 5L))) {
    stop("The fold file must assign 80 exporters to five balanced folds.",
         call. = FALSE)
  }
  fold_table <- fold_table[match(exporters, fold_table$exporter), , drop = FALSE]

  registry <- utils::read.csv(paths[["registry"]], check.names = FALSE)
  required <- c("variant_id", "rta_encoding", "aggregation_weight",
                "tariff_transform", "x_scaling", "response_fit_scale",
                "feature_count")
  if (!all(required %in% names(registry))) {
    stop("The preprocessing registry lacks required specification columns.",
         call. = FALSE)
  }
  spec <- registry[registry$variant_id == "v3_scope_sqrt_raw_center",
                   , drop = FALSE]
  choices <- c(rta_encoding = "scope_0_1_2",
               aggregation_weight = "sqrt_value", tariff_transform = "raw",
               x_scaling = "mean_sd", response_fit_scale = "center_only")
  if (nrow(spec) != 1L ||
      !identical(as.integer(spec$feature_count), 6L) ||
      !all(vapply(names(choices), function(name) {
        identical(as.character(spec[[name]]), choices[[name]])
      }, logical(1L)))) {
    stop("The preprocessing registry differs from the frozen paper design.",
         call. = FALSE)
  }

  list(X = X, Y = Y, fold = as.integer(fold_table$fold),
       exporters = exporters, feature_names = feature_names,
       metadata = metadata, input_hashes = actual_hash)
}
