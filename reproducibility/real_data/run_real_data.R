# Trade-data analysis: Table 1 and the fitting-path table of Appendix D.
#
# Run from the repository root:
#   Rscript reproducibility/real_data/run_real_data.R <prepared_dir> <out_dir> [rank]
# <prepared_dir> holds response_Y.csv, design_X.csv and folds.csv, written by
# 3_build_data.R or taken from the supplementary material of the paper.
# rank defaults to 3.
#
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
OTHER_METHODS <- c("RRR", "Tucker-RGD", "TT-IHT")

# ------------------------------------------------------------------- data
# Y: 80 exporters x 25 partners x 21 HS sections x 8 channels.
y <- read.csv(file.path(data_dir, "response_Y.csv"))
exporters <- unique(y$exporter)
cell <- cbind(match(y$exporter, exporters), match(y$partner, unique(y$partner)),
              match(y$section, unique(y$section)))
channels <- names(y)[-(1:3)]
Y_full <- array(NA_real_, c(length(exporters), OUTPUT_DIMS))
for (k in seq_along(channels)) Y_full[cbind(cell, k)] <- y[[channels[k]]]
# X: 6 features per exporter.
x_table <- read.csv(file.path(data_dir, "design_X.csv"))
X_raw <- as.matrix(x_table[match(exporters, x_table$exporter), -1])   # 80 x 6
folds <- read.csv(file.path(data_dir, "folds.csv"))
fold_id <- folds$fold[match(exporters, folds$exporter)]

# Training/test split and preprocessing of one fold.
prepare_fold <- function(fold) {
  test <- which(fold_id == fold)
  train <- setdiff(seq_along(fold_id), test)
  x_mu <- colMeans(X_raw[train, ])
  x_sd <- apply(X_raw[train, ], 2, sd)
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
  regression_start(dat$X_train, dat$Y_train_array, INPUT_DIMS, OUTPUT_DIMS,
                   RANK, split = 1:2, design_variance = nrow(dat$X_train))
}

# TT-IHT uses the step 1 / lambda_max(X'X), and Tucker-RGD starts from the
# T-HOSVD of the least-squares estimate (real_data = TRUE).
fit_method <- function(method, dat, start, fold, iterations = ITERATIONS,
                       blocks = RANK) {
  if (method == "Tucker-RGN") {
    return(tucker_rgn_tot(dat$X_train, dat$Y_train, INPUT_DIMS, OUTPUT_DIMS,
                          rep(RANK, 4), iterations))
  }
  fit_regression_method(method, dat$X_train, dat$Y_train_array, INPUT_DIMS,
    OUTPUT_DIMS, rank = RANK, start = start, iterations = iterations,
    seed = 2549 + fold,
    step_size = 1 / max(eigen(crossprod(dat$X_train))$values),
    blocks = blocks, real_data = TRUE)
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
    path <- coefficient_path(fit_method(method, dat, start, f,
                                        iterations = LONG_RUN))
    for (k in CHECKPOINTS) {
      accuracy[[length(accuracy) + 1]] <- data.frame(fold = f,
        method = method, checkpoint = k, t(scores(path[[k + 1]], dat)))
    }
  }
  for (blocks in 1:RANK) {
    fit <- fit_method("HOPLS", dat, start, f, blocks = blocks)
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
