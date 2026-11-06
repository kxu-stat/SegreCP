# Trade-data analysis: Table 1 and the fitting-path table of Appendix D.
#
# Run from the repository root after preparing the data (README.md):
#   Rscript reproducibility/real_data/run_real_data.R <data_dir>/analysis_official <out_dir> [rank]
# rank defaults to 3.
#
# Data (written by run.R --prepare):
#   Y: processed/tensors_real.RData, 80 exporters x 25 partners x 21 HS sections
#      x 8 channels;
#   X: processed/preprocessing_sensitivity_inputs/v3_scope_sqrt_raw_center_X_raw.csv
#      (6 features per exporter);
#   folds: outputs/tables/official_rank3_fold_assignments.csv (5 folds of 16).
# In each fold, X is standardized with the training means and standard
# deviations and Y is centered with the training mean of each channel.
#
# Methods: CP-ALS, RGD and RGN from the same CPCA start, RRR, Tucker-RGD,
# Tucker-RGN and TT-IHT, each with rank 3 in its own format and run for 300
# iterations, and HOPLS with 1, 2 and 3 latent blocks. RGD uses the fixed step
# 1 / lambda_max(X'X), i.e. alpha * sigma_X^2 = 1 with sigma_X^2 the largest
# eigenvalue of X'X / n (the six covariates are correlated). RGN and RGD use
# SegreCP with ridge 0.
#
# Output (out_dir):
#   accuracy.csv  one row per fold, method and checkpoint (iteration 3, 5, 10,
#                 30, 300; HOPLS: blocks 1, 2, 3): test RMSE, test R^2 and
#                 train RMSE. Table 1 uses iteration 30.
#   session.txt   sessionInfo().

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(SegreCP)
code <- file.path("reproducibility", "code")
source(file.path(code, "regression_adapters.R"))
source(file.path(code, "tucker_rgn.R"))

RANK <- if (length(args) >= 3) as.integer(args[[3]]) else 3
ITERATIONS <- 30                    # Table 1
LONG_RUN <- 300                     # fitting paths
CHECKPOINTS <- c(3, 5, 10, 30, 300)
INPUT_DIMS <- 6
OUTPUT_DIMS <- c(25, 21, 8)
CP_METHODS <- c("CP-ALS", "RGD", "RGN")
OTHER_METHODS <- c("Tucker-RGD", "TT-IHT")
if (requireNamespace("MultiwayRegression", quietly = TRUE)) {
  OTHER_METHODS <- c("RRR", OTHER_METHODS)
}

# ------------------------------------------------------------------- data
env <- new.env()
load(file.path(data_dir, "processed", "tensors_real.RData"), envir = env)
Y_full <- env$Y                                   # 80 x 25 x 21 x 8
exporters <- env$exporters$exp_iso3
x_table <- read.csv(file.path(data_dir, "processed",
  "preprocessing_sensitivity_inputs", "v3_scope_sqrt_raw_center_X_raw.csv"))
x_table <- x_table[match(exporters, x_table$exp_iso3), ]
stopifnot(identical(x_table$exp_iso3, exporters))
X_raw <- as.matrix(x_table[, -1])                 # 80 x 6
folds <- read.csv(file.path(data_dir, "outputs", "tables",
                            "official_rank3_fold_assignments.csv"))
fold_id <- folds$fold[match(exporters, folds$exporter)]
stopifnot(!anyNA(fold_id), all(table(fold_id) == 16))

# Training/test split and preprocessing of one fold.
prepare_fold <- function(fold) {
  test <- which(fold_id == fold)
  train <- setdiff(seq_along(fold_id), test)
  x_mu <- colMeans(X_raw[train, ])
  x_sd <- pmax(apply(X_raw[train, ], 2, sd), 1e-12)
  y_mu <- apply(Y_full[train, , , ], 4, mean)     # one mean per channel
  center <- function(rows) {
    Y <- sweep(Y_full[rows, , , , drop = FALSE], 4, y_mu, "-")
    matrix(as.numeric(Y), nrow = length(rows))
  }
  Y_train <- center(train)
  list(X_train = scale(X_raw[train, ], x_mu, x_sd),
       X_test = scale(X_raw[test, ], x_mu, x_sd),
       Y_train = Y_train, Y_test = center(test),
       Y_train_array = array(Y_train, c(length(train), OUTPUT_DIMS)))
}

# Coefficient arrays use the input-then-output convention (6 x 25 x 21 x 8).
scores <- function(coefficient, dat) {
  B <- matrix(as.numeric(coefficient), nrow = ncol(dat$X_train))
  test_error <- dat$Y_test - dat$X_test %*% B
  train_error <- dat$Y_train - dat$X_train %*% B
  c(test_rmse = sqrt(mean(test_error^2)),
    test_r2 = 1 - sum(test_error^2) / sum(dat$Y_test^2),
    train_rmse = sqrt(mean(train_error^2)))
}

# CPCA start: CPCA of the moment X'Y / n (sigma_X = 1, since X is
# standardized), with the balanced split of the coefficient modes
# {feature, partner} (6 x 25 = 150, closest to sqrt(25200)).
cpca_start <- function(dat) {
  make_shared_cpca_initialization(dat$X_train, dat$Y_train_array, INPUT_DIMS,
    OUTPUT_DIMS, RANK, split = 1:2, design_variance = nrow(dat$X_train),
    normalize_regression = "none")
}

method_options <- function(diagnostics) {
  options <- regression_trajectory_options()  # ridge 0; RGD fixed step (below)
  options$protocol <- "real_trade"  # TT-IHT step 1/lambda_max(X'X); Tucker-RGD
                                    # starts from OLS + T-HOSVD (as before)
  options$tucker_ranks <- rep(RANK, 4)
  options$tt_ranks <- rep(RANK, 3)
  options$diagnostics <- diagnostics
  options$collect_coefficients <- diagnostics
  options
}

fit_method <- function(method, dat, start, fold, diagnostics,
                       iterations = ITERATIONS, blocks = RANK) {
  if (method == "Tucker-RGN") {
    return(tucker_rgn_tot(dat$X_train, dat$Y_train, INPUT_DIMS, OUTPUT_DIMS,
                          rep(RANK, 4), iterations))
  }
  options <- method_options(diagnostics)
  options$hopls_blocks <- blocks
  if (method == "RGD") {
    options$step_size <- 1 / max(eigen(crossprod(dat$X_train))$values)
  }
  prepared <- prepare_regression_kernel(method, dat$X_train, dat$Y_train_array,
    INPUT_DIMS, OUTPUT_DIMS, rank = RANK,
    shared_init = if (method %in% CP_METHODS) start else NULL,
    iterations = iterations, seed = 2549 + fold, ridge = 0, options = options)
  fit <- fit_regression_kernel(prepared)
  if (diagnostics) fit <- finalize_regression_kernel(fit, prepared)
  fit
}

coefficient_path <- function(fit) {
  if (!is.null(fit$coefficient_history)) return(fit$coefficient_history)
  lapply(fit$component_history, SegreCP::reconstruct_cp)
}

# --------------------------------------------------------------- accuracy
accuracy <- list()
for (f in 1:5) {
  dat <- prepare_fold(f)
  start <- cpca_start(dat)
  for (method in c(CP_METHODS, OTHER_METHODS, "Tucker-RGN")) {
    path <- coefficient_path(fit_method(method, dat, start, f, TRUE,
                                        iterations = LONG_RUN))
    stopifnot(length(path) == LONG_RUN + 1)
    for (k in CHECKPOINTS) {
      accuracy[[length(accuracy) + 1]] <- data.frame(fold = f,
        method = method, checkpoint = k, t(scores(path[[k + 1]], dat)))
    }
  }
  for (blocks in 1:RANK) {
    fit <- fit_method("HOPLS", dat, start, f, TRUE, blocks = blocks)
    accuracy[[length(accuracy) + 1]] <- data.frame(fold = f,
      method = "HOPLS", checkpoint = blocks, t(scores(fit$coefficient, dat)))
  }
}
write.csv(do.call(rbind, accuracy), file.path(out_dir, "accuracy.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------- session
cpu <- if (file.exists("/proc/cpuinfo")) {
  grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]
} else NA
writeLines(c(paste("cpu:", cpu),
             capture.output(sessionInfo())),
           file.path(out_dir, "session.txt"))
