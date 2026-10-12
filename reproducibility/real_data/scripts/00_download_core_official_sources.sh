#!/usr/bin/env bash

# Download the versioned CEPII archives and the World Bank API payloads used by
# the official 2018-to-2019 trade analysis.  Downloads are written to a
# `.partial` file, validated, and moved atomically.  Existing valid files are
# retained, so the script is safe to rerun.

set -euo pipefail
export LC_ALL=C
export LANG=C

data_root="${SEGRECP_REAL_DATA_DIR:?Set SEGRECP_REAL_DATA_DIR to a private data directory}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
verify_script="${script_dir}/verify_source.R"
r_bin="${R_BIN:-$(command -v Rscript || true)}"
[[ -n "${r_bin}" ]] || { echo "Rscript is required" >&2; exit 1; }
archive_dir="${data_root}/source_archives"
raw_dir="${data_root}/raw_official"
metadata_dir="${data_root}/metadata"
log_path="${metadata_dir}/official_download_log.tsv"
force_download="${FORCE_DOWNLOAD:-0}"

mkdir -p "${archive_dir}" "${raw_dir}" "${metadata_dir}"
if [[ ! -f "${log_path}" ]]; then
  printf 'checked_at_utc\taction\tsource_url\tlocal_path\tbytes\tsha256\n' > "${log_path}"
fi

sha256_file() {
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    echo "A SHA-256 utility (shasum or sha256sum) is required." >&2
    return 127
  fi
}

record_file() {
  local action="$1"
  local url="$2"
  local path="$3"
  local rel_path="${path#${data_root}/}"
  local bytes
  bytes="$(wc -c < "${path}" | tr -d ' ')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "${action}" "${url}" \
    "${rel_path}" "${bytes}" "$(sha256_file "${path}")" >> "${log_path}"
}

download_payload() {
  local url="$1"
  local destination="$2"
  local marker="$3"
  local partial="${destination}.partial"

  if [[ "${force_download}" != "1" && -s "${destination}" ]]; then
    grep -q "${marker}" "${destination}" ||
      { echo "Invalid existing response: ${destination}" >&2; return 1; }
    "${r_bin}" "${verify_script}" "${destination}" "${destination#${data_root}/}"
    record_file "verified-existing" "${url}" "${destination}"
    return 0
  fi

  curl --fail --location --silent --show-error --retry 5 --retry-delay 5 \
    --connect-timeout 30 --max-time 600 --output "${partial}" "${url}"
  if [[ ! -s "${partial}" ]] || ! grep -q "${marker}" "${partial}"; then
    echo "Invalid response from ${url}" >&2
    exit 1
  fi
  "${r_bin}" "${verify_script}" "${partial}" "${destination#${data_root}/}"
  mv "${partial}" "${destination}"
  record_file "downloaded" "${url}" "${destination}"
}

download_archive() {
  local url="$1"
  local destination="$2"
  local expected_sha="$3"
  local partial="${destination}.partial"

  if [[ "${force_download}" == "1" || ! -s "${destination}" ]]; then
    curl --fail --location --silent --show-error --retry 5 --retry-delay 5 \
      --connect-timeout 30 --max-time 7200 --output "${partial}" "${url}"
    unzip -tq "${partial}" >/dev/null
    actual_sha="$(sha256_file "${partial}")"
    if [[ "${actual_sha}" != "${expected_sha}" ]]; then
      echo "SHA-256 mismatch for staged archive ${partial}: ${actual_sha}" >&2
      exit 1
    fi
    mv "${partial}" "${destination}"
    action="downloaded"
  else
    unzip -tq "${destination}" >/dev/null
    action="verified-existing"
    actual_sha="$(sha256_file "${destination}")"
    if [[ "${actual_sha}" != "${expected_sha}" ]]; then
      echo "SHA-256 mismatch for existing archive ${destination}: ${actual_sha}" >&2
      exit 1
    fi
  fi
  record_file "${action}" "${url}" "${destination}"
}

baci_url="https://www.cepii.fr/DATA_DOWNLOAD/baci/data/BACI_HS17_V202601.zip"
baci_zip="${archive_dir}/BACI_HS17_V202601.zip"
download_archive "${baci_url}" "${baci_zip}" \
  "16aaa22d3b000cc87744c2c28ac97ce44394c3adf315ec259451d5bb4fcae609"
for member in \
  BACI_HS17_Y2018_V202601.csv \
  BACI_HS17_Y2019_V202601.csv \
  country_codes_V202601.csv \
  product_codes_HS17_V202601.csv; do
  # Re-extract from the hash-verified archive on every run so a stale or
  # partially corrupted extracted file cannot be silently accepted.
  unzip -jo "${baci_zip}" "${member}" -d "${raw_dir}" >/dev/null
  [[ -s "${raw_dir}/${member}" ]] || { echo "Missing extracted ${member}" >&2; exit 1; }
  record_file "verified-extracted" "${baci_url}#${member}" "${raw_dir}/${member}"
done

gravity_url="https://www.cepii.fr/DATA_DOWNLOAD/gravity/data/Gravity_rds_V202211.zip"
gravity_zip="${archive_dir}/Gravity_rds_V202211.zip"
download_archive "${gravity_url}" "${gravity_zip}" \
  "aa5f2e512622c6f429280c6c98c9c0c99fd56bc6efc005158c1fe0d9f11fd60e"
gravity_members=(
  Gravity_V202211.rds Countries_V202211.csv
  Label_gdp_source_d_V202211.csv Label_gdp_source_o_V202211.csv
  Label_legal_new_d_V202211.csv Label_legal_new_o_V202211.csv
  Label_legal_old_d_V202211.csv Label_legal_old_o_V202211.csv
  Label_pop_source_d_V202211.csv Label_pop_source_o_V202211.csv
  Label_rta_coverage_V202211.csv Label_rta_type_V202211.csv
)
for member in "${gravity_members[@]}"; do
  unzip -jo "${gravity_zip}" "${member}" -d "${raw_dir}" >/dev/null
  [[ -s "${raw_dir}/${member}" ]] || { echo "Missing extracted ${member}" >&2; exit 1; }
  record_file "verified-extracted" "${gravity_url}#${member}" "${raw_dir}/${member}"
done

download_payload \
  "https://api.worldbank.org/v2/country/all/indicator/NY.GDP.MKTP.CD?date=2018:2019&format=json&per_page=20000" \
  "${raw_dir}/wdi_gdp_current_usd_2018_2019.json" 'NY.GDP.MKTP.CD'
download_payload \
  "https://api.worldbank.org/v2/country/all/indicator/NY.GDP.PCAP.CD?date=2018:2019&format=json&per_page=20000" \
  "${raw_dir}/wdi_gdp_pc_current_usd_2018_2019.json" 'NY.GDP.PCAP.CD'
download_payload \
  "https://api.worldbank.org/v2/country/all/indicator/SP.POP.TOTL?date=2018:2019&format=json&per_page=20000" \
  "${raw_dir}/wdi_population_2018_2019.json" 'SP.POP.TOTL'
download_payload \
  "https://api.worldbank.org/v2/country/all/indicator/LP.LPI.OVRL.XQ?date=2018&format=json&per_page=20000" \
  "${raw_dir}/world_bank_lpi_overall_2018.json" 'LP.LPI.OVRL.XQ'
download_payload \
  "https://wits.worldbank.org/API/V1/wits/datasource/trn/country/ALL" \
  "${raw_dir}/wits_trains_country_metadata.xml" 'wits:countries'
for year in 2018 2019; do
  download_payload \
    "https://wits.worldbank.org/API/V1/wits/datasource/trn/dataavailability/country/all/year/${year}" \
    "${raw_dir}/wits_trains_dataavailability_${year}.xml" 'wits:dataavailability'
done

echo "Core official sources verified. Download log: ${log_path}"
