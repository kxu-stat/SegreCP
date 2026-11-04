#!/usr/bin/env Rscript
# Fresh signal4 decomposition study through the installed public API.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
source(file.path(dirname(script), "helpers.R"))
run_signal4_simulation("decomposition", script)
