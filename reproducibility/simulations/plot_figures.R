# Figures and simulation tables of the paper.
#
# Run from the repository root:
#   Rscript reproducibility/simulations/plot_figures.R <results_dir> <figure_dir> <table_dir> paper
# results_dir contains the results_<study>.csv files written by
# run_simulations.R (or their .csv.gz copies in output/simulations).
# figure_dir receives the figures (PDF):
#   regression_cp.pdf               Figure 1 (main text)
#   regression_other_formats.pdf    Figure 2 (Appendix C)
#   decomposition.pdf               Figure 3
#   decomposition_gauss_newton.pdf  Figure 4
#   profile_regression.pdf, profile_decomposition.pdf   Figure 5
# table_dir receives the LaTeX tables of Appendix C (timing_summary_*.csv =
# Table 3, table_time_regression.tex = Table 4, table_time_decomposition.tex =
# Table 5, table_errors_regression.tex = Table 6, table_errors_decomposition.tex
# = Table 7, table_initialization_regression.tex = Table 8,
# table_initialization_decomposition.tex = Table 9) and the summary CSV files
# behind them (time_to_level_*.csv, profile_*.csv, rmse_curves_*.csv,
# errors_initialization.csv).
# With "paper", figures are made at the size they are printed: 5.5 inches wide
# (the NeurIPS text width) with 8-point text, so LaTeX does not shrink them;
# appendix figures use most of a page. Otherwise 16 inches wide, 12-point text.
#
# Style: each panel has its own y-axis and one strip label ("SNR 30, rho =
# 0.25"); each method has a colour and a point shape, all lines have the same
# width, and RGD and RGN are drawn last. The points are staggered so that
# coinciding curves stay distinguishable: with M iterative methods in a figure,
# the k-th method marks the iterations k-1, k-1+M, k-1+2M, ... (twice these in
# 100-iteration panels).
# Curves: root mean square over the replications, sqrt(mean(error^2)), of the
# relative Frobenius error; values below 1e-15 are drawn at 1e-15 and values
# above 1 are not drawn.
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
results_dir <- args[[1]]
figure_dir <- args[[2]]
table_dir <- args[[3]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
PAPER <- length(args) >= 4 && args[[4]] == "paper"
# results_<study>.csv, or its compressed copy results_<study>.csv.gz
result_file <- function(study) {
  file <- file.path(results_dir, paste0("results_", study, ".csv"))
  if (file.exists(file)) file else paste0(file, ".gz")
}
WIDTH <- if (PAPER) 5.5 else 16         # inches
BASE <- if (PAPER) 8 else 12            # font size
LW <- if (PAPER) 0.6 else 1             # line widths and point sizes

# Colours and point shapes: within a figure all colours differ (checked for
# colour-vision deficiencies, no black or grey) and so do the shapes.
method_cols <- c(
  "RGD" = "#F8766D", "RGN" = "#00BFC4", "CP-ALS" = "#228B22",
  "Orth-ALS" = "#3ABEFF", "ICO" = "#DAA520",
  "RRR" = "#332288",                                             # Figure 1
  "Tucker-RGD" = "#E69F00", "Tucker-RGN" = "#E7298A",            # Figure 2
  "TT-IHT" = "#332288", "HOPLS" = "#882255",
  "NLS-GN" = "#332288", "GN-CG" = "#E7298A",                     # Figure 4
  "CP-ALS+ELS" = "#E69F00", "Hybrid AMDM" = "#882255")
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

figure_theme <- theme_minimal(base_size = BASE) + theme(
  legend.position = "top",
  strip.background = element_rect(fill = "grey95", color = "grey50"),
  panel.grid.major = element_line(color = "grey85"))

# height: inches in paper mode; the screen version keeps the same shape. The
# fonts are embedded (as the proceedings require) when Ghostscript is found.
save_figure <- function(p, file, height = 7.8) {
  if (!PAPER) height <- height * WIDTH / 5.5
  pdf_file <- file.path(figure_dir, paste0(file, ".pdf"))
  ggsave(pdf_file, p, width = WIDTH, height = height)
  if (nzchar(Sys.which("gs"))) embedFonts(pdf_file, options = "-dPDFSETTINGS=/prepress")
}

rmse_curves <- function(d) {
  d <- subset(d, !is.na(iteration))
  aggregate(error ~ method + coherence + snr + iteration, d,
            function(e) sqrt(mean(e^2)))
}

# Panels may differ in length (decomposition: 100 iterations at coherence 0.75,
# 20 elsewhere); each panel then has its own x axis, with a point every
# iteration in 20-iteration panels and every 10 iterations in longer ones.
convergence_plot <- function(d, methods, file, hline_method = NULL,
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
  # staggered points: method k (0-based, among the iterative methods) at the
  # iterations gap * k + gap * M * j, j = 0, 1, ...
  iterative <- setdiff(methods, hline_method)
  k <- match(as.character(m$method), iterative) - 1
  gap <- ifelse(last > 20, 2, 1)
  pts <- subset(m, !is.na(k) & (iteration - gap * k) %% (gap * length(iterative)) == 0)
  same_length <- length(unique(last)) == 1
  p <- ggplot(m, aes(iteration, log_rmse, colour = method, shape = method,
                     linetype = method, group = method)) +
    geom_blank(data = span, aes(x = 0, y = log_rmse), inherit.aes = FALSE) +
    geom_line(data = subset(m, !method %in% proposed), linewidth = 0.8 * LW) +
    geom_point(data = subset(pts, !method %in% proposed), size = 1.4 * LW,
               na.rm = TRUE) +
    geom_line(data = subset(m, method %in% proposed), linewidth = 0.8 * LW) +
    geom_point(data = subset(pts, method %in% proposed), size = 1.9 * LW) +
    facet_wrap(vars(panel), ncol = 4,
               scales = if (same_length) "free_y" else "free",
               labeller = label_parsed) +
    scale_colour_manual(values = method_cols, breaks = methods) +
    scale_shape_manual(values = method_shapes, breaks = methods) +
    scale_linetype_manual(values = linetypes, breaks = methods) +
    scale_x_continuous(breaks = function(lim) if (lim[2] > 30) seq(0, 100, 25)
                                              else seq(0, 20, 5)) +
    labs(x = "Iteration", y = expression(log[10](RMSE)),
         colour = NULL, shape = NULL, linetype = NULL) +
    guides(colour = guide_legend(nrow = if (length(methods) > 5) 2 else 1)) +
    figure_theme
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
# error in units of 10^-3, with two decimals below 10, one below 100 and none above.
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
convergence_plot(reg, cp_format, "regression_cp", height = 4.3)
convergence_plot(reg, c("Tucker-RGD", "Tucker-RGN", "TT-IHT", "RGN"),
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
  "decomposition")
convergence_plot(dec_plot, c("NLS-GN", "GN-CG", "CP-ALS+ELS", "Hybrid AMDM",
                             "CP-ALS", "RGN"),
  "decomposition_gauss_newton")

# ------------------------------------------------------------------- tables
error_caption <- paste("mean (standard deviation) over 20 replications of the",
  "relative Frobenius error at iteration 20, in units of $10^{-3}$, in the",
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
# Table 3 of the paper reports 1000 x seconds_per_iteration.
for (task in c("regression", "decomposition")) {
  file <- result_file(paste0("timing_", task))
  if (!file.exists(file)) next
  t <- read.csv(file)
  t$seconds_per_iteration <- t$seconds / 20
  t$seconds_per_iteration[t$method %in% c("HOPLS", "CPCA start")] <- NA
  overall <- aggregate(cbind(seconds, seconds_per_iteration) ~ method, t,
                       median, na.action = na.pass)
  write.csv(overall, file.path(table_dir, paste0("timing_summary_", task, ".csv")),
            row.names = FALSE)
}

# ------------------------------------------------ time to the error level
# Per setting with noise: the error level is the lowest RMSE at the last
# iteration (20, or 100 at coherence 0.75 in decomposition) among CP-ALS, RGD
# and RGN; iterations = first iteration with RMSE within 1% of the level (a
# method that never comes within 1% has no row); milliseconds = iterations x the
# method's median time per iteration in the setting (timing study). Written to
# time_to_level_<task>.csv.
time_to_level <- function(d, task) {
  file <- result_file(paste0("timing_", task))
  if (!file.exists(file)) return(NULL)
  t <- read.csv(file)
  t$ms_per_iteration <- 1000 * t$seconds / 20
  cell <- aggregate(ms_per_iteration ~ method + snr + coherence, t, median)
  m <- subset(rmse_curves(d), is.finite(snr))
  last <- ave(m$iteration, m$snr, m$coherence, FUN = max)
  level <- aggregate(error ~ snr + coherence,
    subset(m, iteration == last & method %in% c("CP-ALS", "RGD", "RGN")), min)
  names(level)[3] <- "level"
  m <- merge(m, level)
  first <- aggregate(iteration ~ method + snr + coherence,
                     subset(m, error <= 1.01 * level), min)
  out <- merge(first, cell)
  out$milliseconds <- out$iteration * out$ms_per_iteration
  out <- out[order(out$snr, out$coherence, out$method), ]
  write.csv(out, file.path(table_dir, paste0("time_to_level_", task, ".csv")),
            row.names = FALSE)
  out
}

# LaTeX table of the times to the error level: one row per noise level and
# method, one column per coherence; "--": not reached; bold: fastest in each
# setting (at the shown precision).
time_table <- function(out, methods, label, caption, file) {
  lines <- table_head(caption, label, c("$s$", "Method"))
  for (snr_value in c(30, 10, 3)) {
    lines <- c(lines, "\\midrule")
    for (k in seq_along(methods)) {
      cells <- sapply(c(0, 0.25, 0.5, 0.75), function(rho) {
        x <- out[out$snr == snr_value & out$coherence == rho, ]
        ms <- round(x$milliseconds[x$method == methods[k]])
        if (!length(ms)) return("--")
        best <- min(round(x$milliseconds[x$method %in% methods]))
        if (ms == best) paste0("\\textbf{", ms, "}") else as.character(ms)
      })
      lines <- c(lines, paste0(if (k == 1) snr_value else "", " & ", methods[k],
                               " & ", paste(cells, collapse = " & "), " \\\\"))
    }
  }
  writeLines(c(lines, table_tail), file.path(table_dir, file))
}

# Performance profiles (Dolan and More, 2002): for each method, the fraction
# of the 12 settings with noise in which its value (time to the error level,
# or RMSE at the last iteration) is at most tau times the best value in the
# setting. A method without a value in a setting (it does not reach the level)
# is not within any factor there. Where the best time is 0 (every method starts
# at the level), the methods with time 0 are within factor 1.
profile_rows <- function(x, value, methods, taus) {
  best <- aggregate(x[[value]] ~ snr + coherence, x, min)
  names(best)[3] <- "best"
  x <- merge(x, best)
  x$ratio <- ifelse(x$best == 0, ifelse(x[[value]] == 0, 1, Inf), x[[value]] / x$best)
  rows <- NULL
  for (m in methods) {
    r <- x$ratio[x$method == m]
    rows <- rbind(rows, data.frame(method = m, tau = taus,
      fraction = sapply(taus, function(tau) sum(r <= tau) / 12)))
  }
  rows
}
profile_figure <- function(d, out, methods, file) {
  out$ms <- round(out$milliseconds)
  time <- profile_rows(out, "ms", methods, 2^seq(0, 4, by = 0.05))
  e <- subset(d, is.finite(snr) & method %in% methods)
  last <- ave(e$iteration, e$snr, e$coherence, FUN = function(i) max(i, na.rm = TRUE))
  e <- subset(e, is.na(iteration) | iteration == last)    # HOPLS has no iteration
  e <- aggregate(error ~ method + snr + coherence, e, function(v) sqrt(mean(v^2)))
  err <- profile_rows(e, "error", methods, seq(1, 2, by = 0.005))
  rows <- rbind(data.frame(measure = "Time to the error level", time),
                data.frame(measure = "Error at the last iteration", err))
  rows$measure <- factor(rows$measure, levels = unique(rows$measure))
  rows$method <- factor(rows$method, levels = methods)
  pts <- subset(rows, (seq_len(nrow(rows)) %% 8) == 0)
  shapes <- method_shapes
  shapes["HOPLS"] <- 9                                    # a point shape for HOPLS
  p <- ggplot(rows, aes(tau, fraction, colour = method, shape = method)) +
    geom_step(linewidth = 0.8 * LW) +
    geom_point(data = pts, size = 1.8 * LW, na.rm = TRUE) +
    facet_wrap(vars(measure), scales = "free_x") +
    scale_x_continuous(trans = "log2", breaks = function(lim)
      if (lim[2] > 3) c(1, 2, 4, 8, 16) else c(1, 1.25, 1.5, 2)) +
    scale_colour_manual(values = method_cols, breaks = methods) +
    scale_shape_manual(values = shapes, breaks = methods) +
    labs(x = expression(tau), y = "Fraction of settings",
         colour = NULL, shape = NULL) +
    guides(colour = guide_legend(nrow = if (length(methods) > 5) 2 else 1)) +
    figure_theme
  save_figure(p, file, height = 2.7)
  write.csv(rows, file.path(table_dir, paste0(file, ".csv")), row.names = FALSE)
}

time_reg <- time_to_level(reg, "regression")
time_dec <- time_to_level(dec_plot, "decomposition")
time_caption <- paste("wall-clock time in milliseconds to come within 1\\% of the",
  "error level of the setting, the lowest root mean square error at the last",
  "iteration among CP-ALS, RGD and RGN, in the settings with noise. The time is",
  "the number of iterations needed times the median time per iteration of the",
  "method in the setting. --: not reached. Bold: fastest in each setting.")
if (!is.null(time_reg)) {
  time_table(time_reg, c("RGD", "RGN", "CP-ALS"), "tab:time_regression",
    paste("Tensor-on-tensor regression:", time_caption, "The last iteration is 20.",
          "RRR, Tucker-RGD, Tucker-RGN, TT-IHT and HOPLS do not reach the level",
          "in any setting."),
    "table_time_regression.tex")
  profile_figure(reg, time_reg, c("RGD", "RGN", "CP-ALS", "RRR", "Tucker-RGD",
                                  "Tucker-RGN", "TT-IHT", "HOPLS"), "profile_regression")
}
if (!is.null(time_dec)) {
  time_table(time_dec, c("RGD", "RGN", "CP-ALS", "Orth-ALS", "ICO", "NLS-GN", "GN-CG",
                         "CP-ALS+ELS", "Hybrid AMDM"), "tab:time_decomposition",
    paste("CP decomposition:", time_caption, "The last iteration is 20, or 100 at",
          "$\\rho=0.75$. At $s=30$ and $\\rho=0$, every method starts at the level."),
    "table_time_decomposition.tex")
  profile_figure(dec_plot, time_dec, c("RGD", "RGN", "CP-ALS", "Orth-ALS", "ICO", "NLS-GN",
                                       "GN-CG", "CP-ALS+ELS", "Hybrid AMDM"),
                 "profile_decomposition")
}

# ------------------------------------------------ initialization sensitivity
# Study initialization (settings with noise): CP-ALS, RGD and RGN started from
# CPCA, T-HOSVD and random factors. One LaTeX table per problem: one row per
# noise level, coherence and method; for each start, the mean error at
# iterations 0, 3, 5, 10 and 20 (iteration 0 is the error of the start).
# Means and standard deviations at these iterations: errors_initialization.csv.
init_file <- result_file("initialization")
if (file.exists(init_file)) {
  shown <- c(0, 3, 5, 10, 20)
  starts <- c("CPCA", "T-HOSVD", "Random")
  ini <- subset(read.csv(init_file), iteration %in% shown)
  a <- do.call(data.frame, aggregate(error ~ problem + start + method + snr +
                                       coherence + iteration, ini, mean_sd))
  write.csv(a, file.path(table_dir, "errors_initialization.csv"), row.names = FALSE)
  captions <- c(
    regression = paste("Sensitivity to the initialization in tensor-on-tensor",
      "regression: mean over 20 replications of the relative Frobenius error,",
      "in units of $10^{-3}$, at iterations $t=0,3,5,10,20$ from the CPCA,",
      "T-HOSVD and random starts, in the settings with noise. At $t=0$ the",
      "entry is the error of the start, the same for all three methods."),
    decomposition = paste("Sensitivity to the initialization in CP",
      "decomposition, in the format of Table~\\ref{tab:initialization_regression}."))
  for (problem in names(captions)) {
    lines <- c("\\begin{table}[!htbp]", "\\centering",
      paste0("\\caption{", captions[[problem]], "}"),
      paste0("\\label{tab:initialization_", problem, "}"),
      "{\\small\\setlength{\\tabcolsep}{2.5pt}",
      "\\begin{tabular}{lll*{15}{r}}", "\\toprule",
      paste0(" & & & ", paste0("\\multicolumn{5}{c}{", starts, "}",
                               collapse = " & "), " \\\\"),
      "\\cmidrule(lr){4-8}\\cmidrule(lr){9-13}\\cmidrule(l){14-18}",
      paste0("$s$ & $\\rho$ & Method & ",
             paste(rep(shown, 3), collapse = " & "), " \\\\"))
    for (snr_value in c(30, 10, 3)) {
      lines <- c(lines, "\\midrule")
      for (rho in c(0, 0.25, 0.5, 0.75)) {
        if (rho > 0) lines <- c(lines, "\\addlinespace")
        for (k in 1:3) {
          method <- c("RGD", "RGN", "CP-ALS")[k]
          cells <- sapply(starts, function(start) {
            x <- a[a$problem == problem & a$snr == snr_value &
                   a$coherence == rho & a$start == start & a$method == method, ]
            m <- x$error.mean[order(x$iteration)]
            sprintf("%.*f", decimals(m), m)
          })
          lines <- c(lines, paste0(
            if (rho == 0 && k == 1) snr_value else "", " & ",
            if (k == 1) rho else "", " & ", method, " & ",
            paste(cells, collapse = " & "), " \\\\"))
        }
      }
    }
    writeLines(c(lines, table_tail),
               file.path(table_dir, paste0("table_initialization_", problem, ".tex")))
  }
}
