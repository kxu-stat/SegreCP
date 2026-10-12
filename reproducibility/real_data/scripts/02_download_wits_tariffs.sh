#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C
export LANG=C

data_root="${SEGRECP_REAL_DATA_DIR:?Set SEGRECP_REAL_DATA_DIR to a private data directory}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
verify_script="${script_dir}/verify_source.R"
r_bin="${R_BIN:-$(command -v Rscript || true)}"
[[ -n "${r_bin}" ]] || { echo "Rscript is required" >&2; exit 1; }
year="${WITS_YEAR:-2019}"
if [[ ! "${year}" =~ ^[0-9]{4}$ ]]; then
  echo "WITS_YEAR must be a four-digit year" >&2
  exit 2
fi
reporter_csv="${data_root}/metadata/wits_unique_queries_${year}.csv"
log_path="${data_root}/metadata/official_download_log.tsv"
if [[ ! -f "${log_path}" ]]; then
  printf 'checked_at_utc\taction\tsource_url\tlocal_path\tbytes\tsha256\n' > "${log_path}"
fi
datatype="${WITS_DATATYPE:-reported}"
if [[ "${datatype}" != "reported" && "${datatype}" != "aveestimated" ]]; then
  echo "WITS_DATATYPE must be reported or aveestimated" >&2
  exit 2
fi
if [[ "${datatype}" == "reported" ]]; then
  output_dir="${data_root}/raw_official/wits_tariffs_${year}_by_reporter"
  filename_suffix="${year}"
else
  output_dir="${data_root}/raw_official/wits_tariffs_${year}_aveestimated_by_reporter"
  filename_suffix="aveestimated_${year}"
fi
mkdir -p "${output_dir}"

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

download_one() {
  local reporter_numeric="$1"
  local imp_iso3="$2"
  local destination="${output_dir}/${imp_iso3}_${reporter_numeric}_${filename_suffix}.json"
  local partial="${destination}.partial"
  local url
  url="https://wits.worldbank.org/API/V1/SDMX/V21/datasource/TRN/reporter/${reporter_numeric}/partner/000/product/all/year/${year}/datatype/${datatype}?format=JSON"

  if [[ -s "${destination}" ]]; then
    grep -q 'DF_WITS_Tariff_TRAINS' "${destination}" ||
      { echo "invalid existing WITS response: ${destination}" >&2; return 1; }
    "${r_bin}" "${verify_script}" "${destination}" "${destination#${data_root}/}"
    record_file "verified-existing" "${url}" "${destination}"
    echo "verified existing ${imp_iso3} (${reporter_numeric})"
    return 0
  fi

  curl --fail --location --silent --show-error --retry 4 --retry-delay 5 \
    --connect-timeout 30 --max-time 300 \
    --output "${partial}" "${url}"

  if [[ ! -s "${partial}" ]] || ! grep -q 'DF_WITS_Tariff_TRAINS' "${partial}"; then
    echo "invalid WITS response for ${imp_iso3} (${reporter_numeric})" >&2
    return 1
  fi

  "${r_bin}" "${verify_script}" "${partial}" "${destination#${data_root}/}"
  mv "${partial}" "${destination}"
  record_file "downloaded" "${url}" "${destination}"
  echo "downloaded ${imp_iso3} (${reporter_numeric})"
}

export output_dir
export datatype
export filename_suffix
export year
export data_root
export log_path
export r_bin
export verify_script
export -f sha256_file
export -f record_file
export -f download_one

awk -F, -v datatype="${datatype}" \
  'NR > 1 && (datatype == "reported" || $6 == "Yes") { print $2, $3 }' \
  "${reporter_csv}" |
  xargs -n 2 -P 4 bash -c 'download_one "$1" "$2"' _

expected_count="$(awk -F, -v datatype="${datatype}" \
  'NR > 1 && (datatype == "reported" || $6 == "Yes") { count++ } END { print count + 0 }' \
  "${reporter_csv}")"
actual_count="$(find "${output_dir}" -type f -name "*_${filename_suffix}.json" | wc -l | tr -d ' ')"

if [[ "${actual_count}" != "${expected_count}" ]]; then
  echo "WITS download incomplete: expected ${expected_count}, found ${actual_count}" >&2
  exit 1
fi

echo "WITS ${datatype} download complete: ${actual_count} official reporter files"
