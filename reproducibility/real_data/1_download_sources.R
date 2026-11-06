# Step 1 of the trade data preparation: download the public sources.
#
#   Rscript reproducibility/real_data/1_download_sources.R <data_dir>
#
# Writes to <data_dir>/sources:
#   BACI_HS17_Y2018_V202601.csv, BACI_HS17_Y2019_V202601.csv,
#   country_codes_V202601.csv   CEPII BACI trade flows (HS 2017, version 202601)
#   Gravity_V202211.rds         CEPII Gravity (version 202211)
#   wdi_*.json, lpi_2018.json   World Bank WDI and LPI indicators
#   wits_availability_*.xml     WITS/UNCTAD TRAINS tariff schedules per country
# The archives (about 1 GB) are kept in <data_dir>/archives.

data_dir <- commandArgs(trailingOnly = TRUE)[1]
source_dir <- file.path(data_dir, "sources")
archive_dir <- file.path(data_dir, "archives")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)  # the BACI archive is large

download <- function(url, file) {
  if (!file.exists(file)) download.file(url, file, mode = "wb")
}

# CEPII BACI: bilateral trade flows by HS6 product.
baci_zip <- file.path(archive_dir, "BACI_HS17_V202601.zip")
download("https://www.cepii.fr/DATA_DOWNLOAD/baci/data/BACI_HS17_V202601.zip",
         baci_zip)
unzip(baci_zip, exdir = source_dir, junkpaths = TRUE,
      files = c("BACI_HS17_Y2018_V202601.csv", "BACI_HS17_Y2019_V202601.csv",
                "country_codes_V202601.csv"))

# CEPII Gravity: distance and trade agreements of country pairs.
gravity_zip <- file.path(archive_dir, "Gravity_rds_V202211.zip")
download("https://www.cepii.fr/DATA_DOWNLOAD/gravity/data/Gravity_rds_V202211.zip",
         gravity_zip)
unzip(gravity_zip, files = "Gravity_V202211.rds", exdir = source_dir,
      junkpaths = TRUE)

# World Bank indicators: GDP, GDP per capita and population (2018-2019), and
# the Logistics Performance Index (2018). The World Bank revises these
# indicators; the files used in the paper are in the supplementary material,
# and files already in <data_dir>/sources are not downloaded again.
world_bank <- "https://api.worldbank.org/v2/country/all/indicator/"
indicators <- c(
  wdi_gdp_2018_2019 = "NY.GDP.MKTP.CD?date=2018:2019",
  wdi_gdp_per_capita_2018_2019 = "NY.GDP.PCAP.CD?date=2018:2019",
  wdi_population_2018_2019 = "SP.POP.TOTL?date=2018:2019",
  lpi_2018 = "LP.LPI.OVRL.XQ?date=2018"
)
for (name in names(indicators)) {
  download(paste0(world_bank, indicators[[name]], "&format=json&per_page=20000"),
           file.path(source_dir, paste0(name, ".json")))
}

# WITS: which countries have a tariff schedule in 2018 and 2019.
for (year in c(2018, 2019)) {
  download(paste0("https://wits.worldbank.org/API/V1/wits/datasource/trn/",
                  "dataavailability/country/all/year/", year),
           file.path(source_dir, paste0("wits_availability_", year, ".xml")))
}
