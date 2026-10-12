#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(xml2)
})

data_dir <- Sys.getenv("SEGRECP_REAL_DATA_DIR", unset = "")
if (!nzchar(data_dir)) stop("SEGRECP_REAL_DATA_DIR must name a private data directory")
root <- normalizePath(data_dir, mustWork = TRUE)
raw_dir <- file.path(root, "raw_official")
meta_dir <- file.path(root, "metadata")

candidates <- fread(
  file.path(meta_dir, "wits_reporters_2019.csv"),
  colClasses = list(character = c(
    "baci_country_code", "imp_iso3", "wits_reporter_code", "wits_source_iso3"
  ))
)
candidates[, c(
  "wits_reporter_code", "wits_country_name", "wits_nomenclature",
  "wits_aveestimated_available", "wits_last_updated"
) := NULL]

doc <- read_xml(file.path(raw_dir, "wits_trains_dataavailability_2018.xml"))
nodes <- xml_find_all(doc, ".//wits:reporter", xml_ns(doc))
availability <- data.table(
  wits_reporter_code = xml_attr(nodes, "countrycode"),
  wits_source_iso3 = xml_attr(nodes, "iso3Code"),
  wits_country_name = xml_text(xml_find_first(nodes, "wits:name", xml_ns(doc))),
  wits_nomenclature = xml_attr(
    xml_find_first(nodes, "wits:reporternernomenclature", xml_ns(doc)),
    "reporternernomenclaturecode"
  ),
  wits_aveestimated_available = xml_text(
    xml_find_first(
      nodes,
      "wits:isspecificdutyexpressionestimatedavailable",
      xml_ns(doc)
    )
  ),
  wits_last_updated = xml_text(
    xml_find_first(nodes, "wits:lastupdateddate", xml_ns(doc))
  )
)

candidates_2018 <- merge(
  candidates,
  availability,
  by = "wits_source_iso3",
  all.x = TRUE
)
if (candidates_2018[, anyNA(wits_reporter_code)]) {
  stop("At least one selected importer lacks a WITS 2018 schedule")
}
setorder(candidates_2018, query_order)
fwrite(candidates_2018, file.path(meta_dir, "wits_reporters_2018.csv"))

unique_queries <- unique(candidates_2018[, .(
  wits_reporter_code,
  wits_source_iso3,
  wits_country_name,
  wits_nomenclature,
  wits_aveestimated_available,
  wits_last_updated
)])
setorder(unique_queries, wits_source_iso3)
unique_queries[, wits_query_order := seq_len(.N)]
setcolorder(unique_queries, "wits_query_order")
fwrite(unique_queries, file.path(meta_dir, "wits_unique_queries_2018.csv"))

message(
  "Prepared ", nrow(unique_queries), " WITS 2018 schedules for ",
  nrow(candidates_2018), " selected importers."
)
message(
  "HS2017 schedules: ", unique_queries[wits_nomenclature == "H5", .N],
  "; other schedules to harmonize at HS section level: ",
  unique_queries[wits_nomenclature != "H5", .N]
)
print(unique_queries)
