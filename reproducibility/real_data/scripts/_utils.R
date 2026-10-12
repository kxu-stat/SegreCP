# Shared constants and output helpers for the official-data preprocessing path.
# This file is sourced by the tensor builder; it does not load model code.

data_dir <- Sys.getenv("SEGRECP_REAL_DATA_DIR", unset = "")
if (!nzchar(data_dir)) stop("SEGRECP_REAL_DATA_DIR must name a private data directory")
data_dir <- normalizePath(data_dir, mustWork = TRUE)
analysis_dir <- file.path(data_dir, "analysis_official")

paths <- list(
  raw = file.path(analysis_dir, "raw"),
  processed = file.path(analysis_dir, "processed"),
  logs = file.path(analysis_dir, "outputs", "logs"),
  qc = file.path(analysis_dir, "outputs", "qc"),
  tables = file.path(analysis_dir, "outputs", "tables")
)
for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(seed = 2026L, n_hs1_sections = 21L)
INPUT_CHANNELS <- c(
  "tariff", "ntb_proxy", "log_gdp_partner", "log_distance",
  "rta_indicator", "world_demand"
)
OUTPUT_CHANNELS <- c(
  "log_export_value", "log_export_quantity", "log_unit_value",
  "exporter_share", "partner_market_share", "rca_balassa",
  "extensive_margin", "intensive_margin"
)
HS_SECTIONS <- data.frame(
  hs1 = sprintf("HS%02d", 1:21),
  roman = c(
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
    "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII",
    "XIX", "XX", "XXI"
  ),
  hs2_lo = c(
    1, 6, 15, 16, 25, 28, 39, 41, 44, 47, 50,
    64, 68, 71, 72, 84, 86, 90, 93, 94, 97
  ),
  hs2_hi = c(
    5, 14, 15, 24, 27, 38, 40, 43, 46, 49, 63,
    67, 70, 71, 83, 85, 89, 92, 93, 96, 97
  ),
  label = c(
    "Live animals & animal products", "Vegetable products",
    "Fats & oils", "Prepared foodstuffs", "Mineral products",
    "Chemicals", "Plastics & rubber", "Hides, skins, leather",
    "Wood & cork", "Pulp & paper", "Textiles & apparel",
    "Footwear & headgear", "Stone, ceramics, glass",
    "Pearls & precious metals", "Base metals",
    "Machinery & electrical", "Vehicles & transport",
    "Optical & medical", "Arms & ammunition",
    "Misc manufactured", "Works of art"
  ),
  stringsAsFactors = FALSE
)
stopifnot(length(INPUT_CHANNELS) == 6L,
          length(OUTPUT_CHANNELS) == 8L,
          nrow(HS_SECTIONS) == PARAMS$n_hs1_sections)
set.seed(PARAMS$seed)

log_msg <- function(msg, stage = "pipeline", level = "INFO") {
  line <- sprintf("[%s] [%-5s] %s",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "",
      file = file.path(paths$logs, paste0(stage, ".log")), append = TRUE)
}

.write_meta <- function(target_path, extra = list()) {
  meta <- c(
    list(
      saved_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      r_version = getRversion() |> as.character(),
      pipeline_seed = PARAMS$seed
    ),
    extra
  )
  jsonlite::write_json(
    meta, paste0(target_path, "_meta.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
}

save_csv <- function(df, target_path, extra_meta = list()) {
  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, target_path, row.names = FALSE)
  .write_meta(target_path,
              c(extra_meta, list(rows = nrow(df), cols = ncol(df))))
  invisible(df)
}

save_rdata <- function(target_path, ..., extra_meta = list()) {
  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
  expr <- as.list(match.call(expand.dots = FALSE)$`...`)
  names_to_save <- vapply(expr, function(x) {
    if (!is.symbol(x)) stop("save_rdata accepts bare object names only")
    as.character(x)
  }, character(1))
  save(list = names_to_save, file = target_path, envir = parent.frame())
  .write_meta(target_path, extra_meta)
  invisible(target_path)
}

assert_shape <- function(x, expected, name = deparse(substitute(x))) {
  if (!identical(as.integer(dim(x)), as.integer(expected))) {
    stop(name, " has dimensions ", paste(dim(x), collapse = "x"),
         "; expected ", paste(expected, collapse = "x"))
  }
  invisible(TRUE)
}
