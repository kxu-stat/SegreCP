#!/usr/bin/env Rscript

# Recreate the paper's six-feature trade design and common rank-three folds
# from the official panel. The five-fold split uses 2018 exporter size.
# The v3_scope_sqrt_raw_center basename below is a frozen artifact key, not
# an initialization or estimator name; see reproducibility/real_data/README.md.
suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE)
code_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
source(file.path(code_root, "scripts", "_utils.R"))

panel <- fread(file.path(paths$processed, "panel_real_2018_to_2019.csv"))
tensor_env <- new.env(parent = emptyenv())
load(file.path(paths$processed, "tensors_real.RData"), envir = tensor_env)
exporters <- tensor_env$exporters
E <- as.character(exporters$exp_iso3)

if (!identical(as.integer(dim(tensor_env$Y)), c(80L, 25L, 21L, 8L)) ||
    anyDuplicated(E) || length(E) != 80L) {
  stop("Unexpected official tensor dimensions or exporter identifiers")
}

# Keep the frozen RNG implementation and quintile allocation rule.
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
        sample.kind = "Rejection")
set.seed(2549L)
size_proxy <- log1p(exporters$total_export_kusd_2018)
breaks <- unique(as.numeric(quantile(
  size_proxy, probs = seq(0, 1, length.out = 6L), na.rm = TRUE
)))
strata <- as.integer(cut(
  size_proxy, breaks = breaks, include.lowest = TRUE, labels = FALSE
))
fold_id <- integer(length(E))
strata_levels <- sort(unique(strata))
for (s_pos in seq_along(strata_levels)) {
  idx <- which(strata == strata_levels[s_pos])
  labels <- ((seq_along(idx) - 1L + s_pos - 1L) %% 5L) + 1L
  fold_id[idx] <- sample(labels)
}
if (!identical(as.integer(table(factor(fold_id, levels = 1:5))),
               rep(16L, 5L))) {
  stop("The frozen five-fold allocation is not balanced")
}
fold_assign <- data.frame(
  exporter = E,
  fold = fold_id,
  stratum_2018_volume = strata,
  total_export_kusd_2018 = exporters$total_export_kusd_2018
)
save_csv(
  fold_assign,
  file.path(paths$tables, "official_rank3_fold_assignments.csv"),
  extra_meta = list(seed = 2549L, rule = "2018 exporter-volume quintiles")
)

# Use raw tariff, ordered trade-agreement scope, and square-root 2018 trade weights.
required <- c(
  "exp_iso3", "lag_value_kusd", "tariff", "logistics_friction_proxy",
  "log_gdp_partner", "log_distance", "rta_scope_count", "world_demand"
)
missing <- setdiff(required, names(panel))
if (length(missing)) stop("Missing panel columns: ", paste(missing, collapse = ", "))
if (any(!is.finite(panel$lag_value_kusd)) || any(panel$lag_value_kusd < 0)) {
  stop("2018 trade weights must be finite and nonnegative")
}
panel[, analysis_weight := sqrt(lag_value_kusd)]
columns <- c(
  "tariff", "logistics_friction_proxy", "log_gdp_partner",
  "log_distance", "rta_scope_count", "world_demand"
)
X <- panel[, {
  weights <- analysis_weight
  denominator <- sum(weights)
  as.list(vapply(.SD, function(z) {
    if (denominator > 0) sum(z * weights) / denominator else mean(z)
  }, numeric(1L)))
}, by = exp_iso3, .SDcols = columns]
setnames(X, columns, c(
  "tariff", "ntb_proxy", "log_gdp_partner", "log_distance",
  "rta_scope_count", "world_demand"
))
X <- X[match(E, exp_iso3)]
feature_cols <- setdiff(names(X), "exp_iso3")
if (anyNA(X$exp_iso3) ||
    !identical(as.character(X$exp_iso3), E) ||
    any(!is.finite(as.matrix(X[, ..feature_cols])))) {
  stop("The selected trade design is incomplete or non-finite")
}
design_dir <- file.path(paths$processed, "preprocessing_sensitivity_inputs")
dir.create(design_dir, recursive = TRUE, showWarnings = FALSE)
save_csv(
  X, file.path(design_dir, "v3_scope_sqrt_raw_center_X_raw.csv"),
  extra_meta = list(
    variant_id = "v3_scope_sqrt_raw_center",
    aggregation_weight = "sqrt_value",
    tariff_transform = "raw",
    x_scaling = "mean_sd",
    response_fit_scale = "center_only"
  )
)

# The registry is experiment configuration, not third-party data. Keep its
# frozen 13-row record so the selected archival row retains provenance.
registry_source <- file.path(code_root, "reference_preprocessing_registry.csv")
registry_target <- file.path(paths$tables, "preprocessing_sensitivity_registry.csv")
if (!file.copy(registry_source, registry_target, overwrite = TRUE)) {
  stop("Could not copy the frozen preprocessing registry")
}
.write_meta(registry_target, list(selection_policy = "frozen 13-row registry"))
message("Selected trade design and common five-fold assignment saved under ",
        dirname(paths$processed))
