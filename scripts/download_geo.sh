#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 --geo_id <GEO_ID> --out_dir <OUTDIR>

  --geo_id   GEO series accession (e.g. GSE141044)
  --out_dir  Directory to save downloaded files
  -h         Show this help message
EOF
  exit 1
}

# Parse command-line options
GEO_ID=""
OUTDIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --geo_id)
      GEO_ID="$2"; shift 2 ;;
    --out_dir)
      OUTDIR="$2"; shift 2 ;;
    -h|--help)
      usage ;;
    *)
      echo "Unknown option: $1" >&2
      usage ;;
  esac
done

# Validate parameters
if [[ -z "$GEO_ID" || -z "$OUTDIR" ]]; then
  echo "Error: --geo_id and --out_dir are required." >&2
  usage
fi

# Prepare output directory
mkdir -p "${OUTDIR}/${GEO_ID}"

# GEO supplementary FTP base
FTP_BASE="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${GEO_ID:0:6}nnn/${GEO_ID}/suppl"

echo "Downloading 10× files for ${GEO_ID} into ${OUTDIR}/${GEO_ID}…"

# List of files to fetch
files=(
  "${GEO_ID}_matrix.mtx.gz"
  "${GEO_ID}_Genes.csv.gz"
  "${GEO_ID}_barcode.csv.gz"
  "${GEO_ID}_barcode_sequence.txt.gz"
)

for file in "${files[@]}"; do
  wget -q \
       -O "${OUTDIR}/${GEO_ID}/${file}" \
       "${FTP_BASE}/${file}"
done

echo "Download complete."
