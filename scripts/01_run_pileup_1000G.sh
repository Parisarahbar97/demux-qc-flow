#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 01_run_pileup_1000G.sh --sam BAM --barcodes FILE --out PREFIX [options]

Run popscle dsc-pileup using a population site list (e.g., 1000G common SNPs).
This keeps your donor VCF untouched for demuxlet but gives VAR allele
frequencies stable population AFs as recommended by the Demuxafy docs.

Required:
  --sam PATH        Coordinate-sorted BAM/CRAM with CB/UB tags
  --barcodes FILE   Barcode whitelist (one per line)
  --out PREFIX      Output prefix (produces PREFIX.{plp,var,umi,cel}.gz)

Options:
  --vcf PATH        Population VCF/BCF (bgzip + index)
                    [default: /home/pr422/RDS/live/Users/Parisa/alex_output/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf]
  --tag-group STR   SAM tag for barcodes (default: CB)
  --tag-umi STR     SAM tag for UMIs (default: UB)
  --min-mq INT      Minimum mapping quality (default: 20)
  --min-bq INT      Minimum base quality (default: 13)
  --cap-bq INT      Base quality cap (default: 40)
  --host-root PATH  Host path to bind for Docker (default: /home/pr422)
  --image NAME      popscle Docker image (default: parisa/demux:2.1)
  -h, --help        Show this message
USAGE
}

SAM=""
BARCODES=""
OUT_PREFIX=""
VCF=${POPUL_VCF:-/home/pr422/RDS/live/Users/Parisa/alex_output/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf}
TAG_GROUP=CB
TAG_UMI=UB
MIN_MQ=20
MIN_BQ=13
CAP_BQ=40
HOST_ROOT=${HOST_ROOT:-/home/pr422}
IMAGE=${POPSCLE_IMAGE:-parisa/demux:2.1}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sam) SAM="$2"; shift 2 ;;
    --barcodes) BARCODES="$2"; shift 2 ;;
    --out) OUT_PREFIX="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --tag-group) TAG_GROUP="$2"; shift 2 ;;
    --tag-umi) TAG_UMI="$2"; shift 2 ;;
    --min-mq) MIN_MQ="$2"; shift 2 ;;
    --min-bq) MIN_BQ="$2"; shift 2 ;;
    --cap-bq) CAP_BQ="$2"; shift 2 ;;
    --host-root) HOST_ROOT="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
done

for var in SAM BARCODES OUT_PREFIX VCF; do
  if [[ -z "${!var}" ]]; then
    echo "[ERROR] --${var,,} is required" >&2
    show_help >&2
    exit 1
  fi
done

for f in "$SAM" "$SAM.bai" "$BARCODES" "$VCF"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Missing file: $f" >&2
    exit 1
  fi
done

if [[ ! -f "$VCF.tbi" && ! -f "$VCF.csi" ]]; then
  echo "[ERROR] VCF index (.tbi or .csi) not found for $VCF" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT_PREFIX")"

echo "[INFO] Running dsc-pileup with population VCF: $VCF"
# shellcheck disable=SC2086
docker run --rm -u "$(id -u)":"$(id -g)" \
  -v "$HOST_ROOT":"$HOST_ROOT" \
  --entrypoint /opt/conda/bin/popscle "$IMAGE" dsc-pileup \
    --sam "$SAM" \
    --vcf "$VCF" \
    --group-list "$BARCODES" \
    --tag-group "$TAG_GROUP" \
    --tag-UMI "$TAG_UMI" \
    --min-MQ "$MIN_MQ" \
    --min-BQ "$MIN_BQ" \
    --cap-BQ "$CAP_BQ" \
    --out "$OUT_PREFIX"

echo "[INFO] pileup complete: ${OUT_PREFIX}.{plp,var,umi,cel}.gz"
