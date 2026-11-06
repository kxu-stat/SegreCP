# Figures and simulation tables of the paper.
#
# Run from the repository root:
#   Rscript reproducibility/simulations/plot_figures.R <results_dir> <figure_dir> <table_dir> paper
# results_dir contains the results_<study>.csv files written by
# run_simulations.R (or their .csv.gz copies in output/simulations).
# figure_dir receives the four figures (PDF):
#   regression_cp.pdf               Figure 1
#   regression_other_formats.pdf    Figure 2
#   decomposition.pdf               Figure 3
#   decomposition_gauss_newton.pdf  Figure 4
# table_dir receives the LaTeX tables of Appendix C (table_errors_regression.tex
# = Table 3, table_errors_decomposition.tex = Table 4, timing_summary_*.csv =
# Table 5, table_time_to_level.tex = Table 6, table_initialization.tex =
# Table 7) and the summary CSV files behind them.
# With "paper", figures are made at the size they are printed: 5.5 inches wide
# with 8-point text. Otherwise 16 inches wide with 12-point text.
#
# Curves: root mean square over replications, sqrt(mean(error^2)), of the
# relative Frobenius error; values below 1e-15 are drawn at 1e-15 and values
# above 1 are not drawn.
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
results_dir <- args[[1]]
# results_<study>.csv, or its compressed copy results_<study>.csv.gz
result_file <- function(study) {
  file <- file.path(results_dir, paste0("results_", study, ".csv"))
  if (file.exists(file)) file else paste0(file, ".gz")
}
figure_dir <- args[[2]]
table_dir <- args[[3]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
PAPER <- length(args) >= 4 && args[[4]] == "paper"
WIDTH <- if (PAPER) 5.5 else 16         # inches
BASE <- if (PAPER) 8 else 12            # font size
LW <- if (PAPER) 0.6 else 1             # line widths and point sizes

# One colour and one point shape per method.
method_cols <- c(
  "RGD" = "#F8766D", "RGN" = "#00BFC4", "CP-ALS" = "#228B22",
  "Orth-ALS" = "#3ABEFF", "ICO" = "#DAA520",
  "NLS-GN" = "#7B3294", "GN-CG" = "#E7298A", "CP-ALS+ELS" = "#8C510A",
  "Hybrid AMDM" = "#737373",
  "RRR" = "#08306B", "Tucker-RGD" = "#FDB462", "Tucker-RGN" = "#D95F02",
  "TT-IHT" = "#000000", "HOPLS" = "#B15928")
method_shapes <- c(
  "RGD" = 16, "RGN" = 17, "CP-ALS" = 15, "Orth-ALS" = 3, "ICO" = 4,
  "NLS-GN" = 18, "GN-CG" = 8, "CP-ALS+ELS" = 0, "Hybrid AMDM" = 1,
  "RRR" = 2, "Tucker-RGD" = 5, "Tucker-RGN" = 6, "TT-IHT" = 7, "HOPLS" = NA)
proposed <- c("RGD", "RGN")

snr_lab <- function(s) factor(ifelse(is.infinite(s), "Noiseless",
                                     paste0("SNR: ", s)),
  levels = c("Noiseless", "SNR: 30", "SNR: 10", "SNR: 3"))
coh_lab <- function(x) factor(paste0("Coherence: ", x),
  levels = paste0("Coherence: ", c(0, 0.25, 0.5, 0.75)))
# One strip label per panel, e.g. "SNR 30, rho = 0.25" (plotmath, for
# label_parsed), ordered by noise level and then coherence.
panel_lab <- function(s, coherence) {
  noise <- ifelse(is.infinite(s), '"Noiseless"', paste0("SNR~", s))
  label <- paste0(noise, '*","~rho==', coherence)
  grid <- expand.grid(coherence = c(0, 0.25, 0.5, 0.75), s = c(Inf, 30, 10, 3))
  levels <- paste0(ifelse(is.infinite(grid$s), '"Noiseless"', paste0("SNR~", grid$s)),
                   '*","~rho==', grid$coherence)
  factor(label, levels = levels[levels %in% label])
}

old_theme <- theme_minimal(base_size = BASE) + theme(
  legend.position = "top",
  plot.title = element_text(hjust = 0.5, face = "bold"),
  strip.background = element_rect(fill = "grey95", color = "grey50"),
  panel.grid.major = element_line(color = "grey85"))

# height: inches in paper mode; the screen version keeps the same shape.
save_figure <- function(p, file, height = 7.8) {
  if (!PAPER) height <- height * WIDTH / 5.5
  ggsave(file.path(figure_dir, paste0(file, ".pdf")), p, width = WIDTH,
         height = height)
}

rmse_curves <- function(d) {
  d <- subset(d, !is.na(iteration))
  aggregate(error ~ method + coherence + snr + iteration, d,
            function(e) sqrt(mean(e^2)))
}

# Panels may differ in length (decomposition: 100 iterations at coherence 0.75,
# 20 elsewhere); each panel then has its own x axis, with a point every
# iteration in 20-iteration panels and every 10 iterations in longer ones.
convergence_plot <- function(d, methods, title, file, hline_method = NULL,
                             height = 7.8) {
  m <- rmse_curves(subset(d, method %in% methods))
  if (!is.null(hline_method)) {        # a non-iterative method: flat dashed line
    h <- subset(d, method == hline_method)
    h <- aggregate(error ~ method + coherence + snr, h, function(e) sqrt(mean(e^2)))
    h <- merge(h, data.frame(iteration = 0:max(m$iteration)))
    m <- rbind(m, h[names(m)])
    methods <- c(methods, hline_method)
  }
  # Errors above 1 (worse than the zero estimate; only the random start of RRR)
  # are not drawn, so that they do not stretch the axes.
  m <- subset(m, error <= 1)
  m$log_rmse <- log10(pmax(m$error, 1e-15))
  m$panel <- panel_lab(m$snr, m$coherence)
  m$method <- factor(m$method, levels = methods)

  # Each panel spans at least 0.6 decades, so flat panels are not blown up.
  span <- aggregate(log_rmse ~ panel, m, function(y) mean(range(y)))
  span <- rbind(transform(span, log_rmse = log_rmse - 0.3),
                transform(span, log_rmse = log_rmse + 0.3))

  linetypes <- setNames(ifelse(methods == "HOPLS", "dashed", "solid"), methods)
  last <- ave(m$iteration, m$panel, FUN = max)
  pts <- subset(m, iteration %% ifelse(last > 20, 10, 1) == 0)
  same_length <- length(unique(last)) == 1
  p <- ggplot(m, aes(iteration, log_rmse, colour = method, shape = method,
                     linetype = method, group = method)) +
    geom_blank(data = span, aes(x = 0, y = log_rmse), inherit.aes = FALSE) +
    geom_line(data = subset(m, !method %in% proposed), linewidth = 0.6 * LW) +
    geom_point(data = subset(pts, !method %in% proposed), size = 1.4 * LW,
               na.rm = TRUE) +
    geom_line(data = subset(m, method %in% proposed), linewidth = 1.1 * LW) +
    geom_point(data = subset(pts, method %in% proposed), size = 1.9 * LW) +
    facet_wrap(vars(panel), ncol = 4,
               scales = if (same_length) "free_y" else "free",
               labeller = label_parsed) +
    scale_colour_manual(values = method_cols, breaks = methods) +
    scale_shape_manual(values = method_shapes, breaks = methods) +
    scale_linetype_manual(values = linetypes, breaks = methods) +
    scale_x_continuous(breaks = function(lim) if (lim[2] > 30) seq(0, 100, 25)
                                              else seq(0, 20, 5)) +
    labs(title = title, x = "Iteration", y = expression(log[10](RMSE)),
         colour = "Method", shape = "Method", linetype = "Method") +
    guides(colour = guide_legend(nrow = if (length(methods) > 5) 2 else 1)) +
    old_theme
  save_figure(p, file, height)
  invisible(m)
}

# A method that stops before the last iteration (RRR when its updates stop
# changing the fit) keeps its last estimate: the path is extended with its
# last error, so that every curve and table entry averages all replications.
fill_paths <- function(d) {
  it <- subset(d, !is.na(iteration))
  last <- max(it$iteration)
  key <- c("method", "snr", "coherence", "replication")
  ends <- it[order(-it$iteration), ]
  ends <- ends[!duplicated(ends[key]) & ends$iteration < last, ]
  if (!nrow(ends)) return(d)
  extra <- do.call(rbind, lapply(seq_len(nrow(ends)), function(k) {
    e <- ends[rep(k, last - ends$iteration[k]), ]
    e$iteration <- (ends$iteration[k] + 1):last
    e
  }))
  rbind(d, extra)
}

# Table entries: mean (standard deviation) over replications of the relative
# error times 10^3, with two decimals below 10, one below 100 and none above.
mean_sd <- function(e) c(mean = 1000 * mean(e), sd = 1000 * sd(e))
decimals <- function(m) ifelse(m < 10, 2, ifelse(m < 100, 1, 0))
cell_text <- function(m, sd, bold = rep(FALSE, length(m))) {
  k <- decimals(m)
  value <- sprintf("%.*f", k, m)
  value[bold] <- paste0("\\textbf{", value[bold], "}")
  paste0(value, " (", sprintf("%.*f", k, sd), ")")
}
table_head <- function(caption, label, first) {
  c("\\begin{table}[!htbp]", "\\centering", paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"), "{\\small",
    paste0("\\begin{tabular}{", strrep("l", length(first)), "cccc}"), "\\toprule",
    paste0(paste(first, collapse = " & "),
           " & $\\rho=0$ & $\\rho=0.25$ & $\\rho=0.5$ & $\\rho=0.75$ \\\\"))
}
table_tail <- c("\\bottomrule", "\\end{tabular}}", "\\end{table}")

# LaTeX table of the errors at iteration 20 in the settings with noise (without
# noise they are at machine precision): one row per noise level and method, one
# column per coherence; bold: lowest mean in each setting at the shown precision.
error_table <- function(d, methods, label, caption, file) {
  e <- subset(d, is.finite(snr) & (iteration == 20 | is.na(iteration)) &
                 method %in% methods)
  a <- do.call(data.frame, aggregate(error ~ method + snr + coherence, e, mean_sd))
  a$shown <- as.numeric(sprintf("%.*f", decimals(a$error.mean), a$error.mean))
  a$bold <- a$shown == ave(a$shown, a$snr, a$coherence, FUN = min)
  lines <- table_head(caption, label, c("$s$", "Method"))
  for (snr_value in c(30, 10, 3)) {
    lines <- c(lines, "\\midrule")
    for (k in seq_along(methods)) {
      x <- a[a$snr == snr_value & a$method == methods[k], ]
      x <- x[order(x$coherence), ]
      lines <- c(lines, paste0(if (k == 1) snr_value else "", " & ", methods[k], " & ",
        paste(cell_text(x$error.mean, x$error.sd, x$bold), collapse = " & "),
        " \\\\"))
    }
  }
  writeLines(c(lines, table_tail), file.path(table_dir, file))
}

# ------------------------------------------------------------------ regression
reg <- fill_paths(read.csv(result_file("regression")))
cp_format <- intersect(c("CP-ALS", "RRR", "RGD", "RGN"), reg$method)
convergence_plot(reg, cp_format,
  "Convergence of Methods for Tensor Regression - Rank 3",
  "regression_cp", height = 4.3)
convergence_plot(reg, c("Tucker-RGD", "Tucker-RGN", "TT-IHT", "RGN"),
  "Convergence of Methods for Tensor Regression - Rank 3",
  "regression_other_formats", hline_method = "HOPLS")

# --------------------------------------------------------------- decomposition
dec <- fill_paths(read.csv(result_file("decomposition")))
# At coherence 0.75 the figures show the 100-iteration run (study
# decomposition_long: same data sets, so its first 20 iterations equal those
# of the main study).
dec_plot <- dec
long_file <- result_file("decomposition_long")
if (file.exists(long_file)) {
  long <- fill_paths(read.csv(long_file))
  dec_plot <- rbind(subset(dec, coherence != 0.75), long[names(dec)])
  write.csv(rmse_curves(long),
            file.path(table_dir, "rmse_curves_decomposition_long.csv"),
            row.names = FALSE)
}
convergence_plot(dec_plot, c("RGD", "RGN", "CP-ALS", "Orth-ALS", "ICO"),
  "Convergence of Methods for Tensor Decomposition - Rank 3",
  "decomposition")
convergence_plot(dec_plot, c("NLS-GN", "GN-CG", "CP-ALS+ELS", "Hybrid AMDM",
                             "CP-ALS", "RGN"),
  "Convergence of Methods for Tensor Decomposition - Rank 3",
  "decomposition_gauss_newton")

# ------------------------------------------------------------------- tables
error_caption <- paste("mean (standard deviation) over 20 replications of the",
  "relative Frobenius error at iteration 20, multiplied by $10^3$, in the",
  "settings with noise. Bold: lowest mean in each setting.")
error_table(reg, c("RGD", "RGN", "CP-ALS", "RRR", "Tucker-RGD", "Tucker-RGN",
                   "TT-IHT", "HOPLS"),
  "tab:regression_errors",
  paste("Tensor-on-tensor regression:", error_caption, "HOPLS is not iterative."),
  "table_errors_regression.tex")
error_table(dec, c("RGD", "RGN", "CP-ALS", "Orth-ALS", "ICO", "NLS-GN", "GN-CG",
                   "CP-ALS+ELS", "Hybrid AMDM"),
  "tab:decomposition_errors", paste("CP decomposition:", error_caption),
  "table_errors_decomposition.tex")
write.csv(rmse_curves(reg), file.path(table_dir, "rmse_curves_regression.csv"),
          row.names = FALSE)
write.csv(rmse_curves(dec), file.path(table_dir, "rmse_curves_decomposition.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------- timing
# Median wall-clock seconds over all settings from the timing studies
# (diagnostics off, one thread, 20 iterations per fit); seconds_per_iteration =
# seconds / 20 (not defined for HOPLS, which is not iterative, or the CPCA start).
# Time to the error level, per setting with noise. The level is the lowest RMSE
# at the last iteration shown in the figures (20, or 100 at coherence 0.75 in
# decomposition) among CP-ALS, RGD and RGN; iterations = first iteration with
# RMSE within 1% of the level (NA if never); milliseconds = iterations x the
# median time per iteration of the method in that setting (timing study).
timing_table <- function(task, d) {
  file <- result_file(paste0("timing_", task))
  if (!file.exists(file)) return(NULL)
  t <- read.csv(file)
  t$seconds_per_iteration <- t$seconds / 20
  t$seconds_per_iteration[t$method %in% c("HOPLS", "CPCA start")] <- NA
  overall <- aggregate(cbind(seconds, seconds_per_iteration) ~ method, t,
                       median, na.action = na.pass)
  write.csv(overall, file.path(table_dir, paste0("timing_summary_", task, ".csv")),
            row.names = FALSE)
  m <- subset(rmse_curves(d), is.finite(snr))
  m$last <- ave(m$iteration, m$snr, m$coherence, FUN = max)
  level <- aggregate(error ~ snr + coherence,
    subset(m, iteration == last & method %in% c("CP-ALS", "RGD", "RGN")), min)
  m <- merge(m, setNames(level, c("snr", "coherence", "level")))
  first <- aggregate(iteration ~ method + snr + coherence,
                     subset(m, error <= 1.01 * level), min)
  cell <- aggregate(seconds_per_iteration ~ method + snr + coherence, t, median)
  out <- merge(first, cell)
  out$milliseconds <- 1000 * out$iteration * out$seconds_per_iteration
  out <- out[order(out$snr, out$coherence, out$method), ]
  write.csv(out, file.path(table_dir, paste0("time_to_level_", task, ".csv")),
            row.names = FALSE)
  out
}
time_dec <- timing_table("decomposition", dec_plot)
time_reg <- timing_table("regression", reg)

# LaTeX table of the times to the error level: one row per problem, noise
# level and method (only methods that reach the level), one column per
# coherence; "--" = not reached; bold = fastest in each setting.
if (!is.null(time_dec) && !is.null(time_reg)) {
  lines <- table_head(paste("Wall-clock time in milliseconds to come within",
      "1\\% of the error level, the lowest root mean square error at the last",
      "iteration among CP-ALS, RGD and RGN (iteration 20, or 100 at $\\rho=0.75$",
      "in decomposition), in the settings with noise. --: not reached. Bold:",
      "fastest in each setting. RRR and the regression methods in other tensor",
      "formats reach the level in no setting, and Orth-ALS and ICO only at",
      "$\\rho=0$."),
    "tab:time_to_level", c("$s$", "Method"))
  blocks <- list(list(name = "Tensor-on-tensor regression", d = time_reg,
                      methods = c("RGD", "RGN", "CP-ALS")),
                 list(name = "CP decomposition", d = time_dec,
                      methods = c("RGD", "RGN", "CP-ALS", "NLS-GN", "GN-CG",
                                  "CP-ALS+ELS", "Hybrid AMDM")))
  for (b in blocks) {
    lines <- c(lines, "\\midrule", paste0("\\multicolumn{6}{l}{\\emph{", b$name,
                                          "}} \\\\"))
    for (snr_value in c(30, 10, 3)) {
      if (snr_value != 30) lines <- c(lines, "\\addlinespace")
      for (k in seq_along(b$methods)) {
        cells <- sapply(c(0, 0.25, 0.5, 0.75), function(rho) {
          x <- b$d[b$d$snr == snr_value & b$d$coherence == rho, ]
          ms <- x$milliseconds[x$method == b$methods[k]]
          if (!length(ms)) return("--")
          value <- sprintf("%.0f", ms)
          best <- min(round(x$milliseconds[x$method %in% b$methods]))
          if (round(ms) == best) paste0("\\textbf{", value, "}") else value
        })
        lines <- c(lines, paste0(if (k == 1) snr_value else "", " & ",
          b$methods[k], " & ", paste(cells, collapse = " & "), " \\\\"))
      }
    }
  }
  writeLines(c(lines, table_tail), file.path(table_dir, "table_time_to_level.tex"))
}

# ------------------------------------------------ initialization sensitivity
# Study initialization (settings with noise): CP-ALS, RGD and RGN started from
# CPCA, T-HOSVD and random factors. One LaTeX table in the format of the error
# tables (one row per problem, noise level, start and method); the CPCA start
# is already in the error tables.
init_file <- result_file("initialization")
if (file.exists(init_file)) {
  ini <- subset(read.csv(init_file), iteration == 20)
  a <- do.call(data.frame, aggregate(error ~ problem + start + method + snr +
                                       coherence, ini, mean_sd))
  write.csv(a, file.path(table_dir, "errors_initialization.csv"), row.names = FALSE)
  lines <- table_head(paste("Sensitivity to the initialization: mean (standard",
      "deviation) over 20 replications of the relative Frobenius error at",
      "iteration 20, multiplied by $10^3$, from the T-HOSVD start and a random",
      "start. The CPCA start is in Tables~\\ref{tab:regression_errors}",
      "and~\\ref{tab:decomposition_errors}."),
    "tab:initialization", c("$s$", "Start", "Method"))
  for (problem in c("regression", "decomposition")) {
    lines <- c(lines, "\\midrule", paste0("\\multicolumn{7}{l}{\\emph{",
      if (problem == "regression") "Tensor-on-tensor regression" else
        "CP decomposition", "}} \\\\"))
    for (snr_value in c(30, 10, 3)) {
      if (snr_value != 30) lines <- c(lines, "\\addlinespace")
      for (start in c("T-HOSVD", "Random")) {
        for (k in 1:3) {
          method <- c("RGD", "RGN", "CP-ALS")[k]
          x <- a[a$problem == problem & a$snr == snr_value & a$start == start &
                 a$method == method, ]
          x <- x[order(x$coherence), ]
          lines <- c(lines, paste0(
            if (start == "T-HOSVD" && k == 1) snr_value else "", " & ",
            if (k == 1) start else "", " & ", method, " & ",
            paste(cell_text(x$error.mean, x$error.sd), collapse = " & "), " \\\\"))
        }
      }
    }
  }
  writeLines(c(lines, table_tail), file.path(table_dir, "table_initialization.tex"))
}
