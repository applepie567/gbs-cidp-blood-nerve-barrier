#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
output_root="$project_root/data/raw/PRJNA1293757"
mkdir -p "$output_root"

# The public Mendeley deposits contain complete Cell Ranger output archives
# (including large BAM files).  The filtered_feature_bc_matrix.h5 file is the
# first archive member.  A bounded HTTP range therefore retrieves only the
# analysis-ready matrix and avoids downloading irrelevant BAM data.
download_matrix() {
  local sample="$1"
  local url="$2"
  local sample_dir="$output_root/$sample"
  local temp_dir
  temp_dir="$(mktemp -d)"
  mkdir -p "$sample_dir"

  set +e
  curl -L --fail --retry 3 --silent --show-error --range 0-70000000 "$url" \
    | tar -xzf - -C "$temp_dir" outs/filtered_feature_bc_matrix.h5 2>/dev/null
  local pipeline_status=("${PIPESTATUS[@]}")
  set -e

  # A truncated-archive error is expected after the requested first member;
  # integrity is checked by the HDF5 reader in the analysis script.
  if [[ ! -s "$temp_dir/outs/filtered_feature_bc_matrix.h5" ]]; then
    printf 'Failed to extract %s (curl=%s, tar=%s)\n' \
      "$sample" "${pipeline_status[0]}" "${pipeline_status[1]}" >&2
    return 1
  fi
  mv "$temp_dir/outs/filtered_feature_bc_matrix.h5" \
    "$sample_dir/filtered_feature_bc_matrix.h5"
  rmdir "$temp_dir/outs" "$temp_dir"
  printf '%s\t%s bytes\n' "$sample" \
    "$(stat -c '%s' "$sample_dir/filtered_feature_bc_matrix.h5")"
}

download_matrix HPBMC1 \
  'https://data.mendeley.com/public-files/datasets/f2523nkt3m/files/76364a09-c551-4648-80b2-2ae4d5e868b9/file_downloaded' &
download_matrix HPBMC2 \
  'https://data.mendeley.com/public-files/datasets/f2523nkt3m/files/e6f0e310-82ee-481b-b33e-a3af41cd5e5d/file_downloaded' &
download_matrix PBMC1 \
  'https://data.mendeley.com/public-files/datasets/f2523nkt3m/files/99b434b4-dc24-4e20-8d78-bb4c12da9231/file_downloaded' &
download_matrix PBMC2 \
  'https://data.mendeley.com/public-files/datasets/f2523nkt3m/files/8dedb442-be70-43da-8de6-ba12edc82114/file_downloaded' &
download_matrix PBMC3 \
  'https://data.mendeley.com/public-files/datasets/f2523nkt3m/files/bd84c8f0-640c-4d43-b410-2bc5e979ec20/file_downloaded' &
wait

