#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
raw_dir="${project_root}/data/raw/GSE285983"
h5_dir="${raw_dir}/h5"

mkdir -p "${raw_dir}" "${h5_dir}"

curl -L --fail --retry 3 -C - \
  -o "${raw_dir}/GSE285983_metadata_all.csv.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285983/suppl/GSE285983_metadata_all.csv.gz"

curl -L --fail --retry 3 -C - \
  -o "${raw_dir}/GSE285983_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285983/suppl/GSE285983_RAW.tar"

tar -xf "${raw_dir}/GSE285983_RAW.tar" -C "${h5_dir}"

