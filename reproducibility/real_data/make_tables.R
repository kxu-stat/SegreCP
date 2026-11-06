# LaTeX rows of Table 1 and of the fitting-path table (Table 8), from the
# output of run_real_data.R.
# Run from the repository root:
#   Rscript reproducibility/real_data/make_tables.R <out_dir>
# Writes to <out_dir>:
#   table1_rows.tex  rows of Table 1 (test RMSE and R^2 after 30 iterations),
#                    CP-format methods first, then the other formats, each group
#                    under its own heading; the lowest mean test RMSE of each
#                    group is in bold;
#   paths_rows.tex   rows of Table 8 (iterations 3, 5, 10, 30, 300) and the
#                    HOPLS row (1, 2, 3 latent blocks).
# A method that was not run (RRR without MultiwayRegression) gets dashes.
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1]]
accuracy <- read.csv(file.path(out_dir, "accuracy.csv"))

se <- function(x) sd(x) / sqrt(length(x))
pm <- function(x) sprintf("%.4f\\pm%.4f", mean(x), se(x))

cp_methods <- c("CP-ALS", "RRR", "RGD", "RGN")
other_methods <- c("Tucker-RGD", "Tucker-RGN", "TT-IHT", "HOPLS")
model <- c("CP-ALS" = "CP rank $=3$", "RRR" = "CP rank $=3$",
           "RGD" = "CP rank $=3$", "RGN" = "CP rank $=3$",
           "Tucker-RGD" = "Tucker rank $=(3,3,3,3)$",
           "Tucker-RGN" = "Tucker rank $=(3,3,3,3)$",
           "TT-IHT" = "TT rank $=(3,3,3)$",
           "HOPLS" = "3 latent blocks")
# Fit after 30 iterations (HOPLS: 3 blocks).
terminal <- subset(accuracy, checkpoint == ifelse(method == "HOPLS", 3, 30))

group_rows <- function(methods, heading) {
  run <- methods[methods %in% terminal$method]
  means <- sapply(run, function(m) mean(subset(terminal, method == m)$test_rmse))
  best <- names(means)[which.min(means)]
  lines <- sprintf("\\multicolumn{4}{l}{\\textit{%s}} \\\\", heading)
  for (m in methods) {
    x <- subset(terminal, method == m)
    if (!nrow(x)) {
      lines <- c(lines, sprintf("%-10s & %s & -- & -- \\\\ %% not run",
                                m, model[[m]]))
      next
    }
    rmse <- pm(x$test_rmse)
    if (m == best) rmse <- paste0("\\mathbf{", rmse, "}")
    lines <- c(lines, sprintf("%-10s & %s & $%s$ & $%s$ \\\\", m,
                              model[[m]], rmse, pm(x$test_r2)))
  }
  lines
}
table1 <- c(group_rows(cp_methods, "CP format"), "\\midrule",
            group_rows(other_methods, "Other tensor formats"))
writeLines(table1, file.path(out_dir, "table1_rows.tex"))

path_row <- function(m, checkpoints) {
  x <- subset(accuracy, method == m)
  if (!nrow(x)) {
    return(paste0(sprintf("%-10s", m),
                  strrep(" & --", length(checkpoints)), " \\\\ % not run"))
  }
  v <- sapply(checkpoints, function(k) mean(subset(x, checkpoint == k)$test_rmse))
  paste0(sprintf("%-10s", m), paste(sprintf(" & $%.4f$", v), collapse = ""), " \\\\")
}
iters <- c(3, 5, 10, 30, 300)
paths <- c(sapply(cp_methods, path_row, checkpoints = iters), "\\midrule",
           sapply(setdiff(other_methods, "HOPLS"), path_row, checkpoints = iters),
           "\\midrule",
           "& \\multicolumn{3}{c}{Panel B: HOPLS latent blocks} \\\\",
           "\\cmidrule(lr){2-4}",
           "Method & 1 & 2 & 3 \\\\",
           "\\midrule",
           path_row("HOPLS", 1:3))
writeLines(paths, file.path(out_dir, "paths_rows.tex"))
cat(table1, "", paths, sep = "\n")
