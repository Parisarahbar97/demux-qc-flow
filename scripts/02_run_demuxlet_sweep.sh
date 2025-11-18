#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 02_run_demuxlet_sweep.sh --plp PREFIX --vcf VCF --barcodes FILE --outdir DIR [options]

Run popscle demuxlet (GT mode) on a population-VCF pileup and sweep a set of
fixed genotype error offsets (default: 0.01,0.02,0.03,0.05). This follows the
Demuxafy recommendation to combine a population site list (for pileup VAR AF)
with slightly larger genotype error tolerance for GT-only donor VCFs.

Required:
  --plp PREFIX       Output prefix from dsc-pileup (without .plp.gz suffix)
  --vcf PATH         Pool-specific donor VCF/BCF (bgzip + index, BAM-ordered)
  --barcodes FILE    Whitelist used for pileup
  --outdir DIR       Directory for demuxlet runs (one subprefix per error)

Options:
  --errs CSV         Comma-separated geno-error offsets (default: 0.01,0.02,0.03,0.05)
  --doublet-prior F  Prior doublet rate (default: 0.05)
  --alpha-list CSV   Alpha grid (default: 0,0.1,0.2,0.3,0.4,0.5)
  --host-root PATH   Host path to bind for Docker (default: /home/pr422)
  --image NAME       popscle Docker image (default: parisa/demux:2.1)
  -h, --help         Show this message
USAGE
}

PLP=""
VCF=""
BARCODES=""
OUTDIR=""
ERRS=${DEMUX_ERRS:-0.15}
DOUBLETP=0.05
ALPHAS=${DEMUX_ALPHAS:-0,0.5}
HOST_ROOT=${HOST_ROOT:-/home/pr422}
IMAGE=${POPSCLE_IMAGE:-parisa/demux:2.1}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --plp) PLP="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --barcodes) BARCODES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --errs) ERRS="$2"; shift 2 ;;
    --doublet-prior) DOUBLETP="$2"; shift 2 ;;
    --alpha-list) ALPHAS="$2"; shift 2 ;;
    --host-root) HOST_ROOT="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
done

for var in PLP VCF BARCODES OUTDIR; do
  if [[ -z "${!var}" ]]; then
    echo "[ERROR] --${var,,} is required" >&2
    show_help >&2
    exit 1
  fi
done

for f in "$PLP.plp.gz" "$PLP.var.gz" "$PLP.umi.gz" "$PLP.cel.gz" "$VCF" "$BARCODES"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Missing file: $f" >&2
    exit 1
  fi
done

if [[ ! -f "$VCF.tbi" && ! -f "$VCF.csi" ]]; then
  echo "[ERROR] VCF index (.tbi/.csi) missing for $VCF" >&2
  exit 1
fi

IFS=',' read -r -a ERR_ARR <<< "$ERRS"
IFS=',' read -r -a ALPHA_ARR <<< "$ALPHAS"
ALPHA_ARGS=()
for a in "${ALPHA_ARR[@]}"; do
  ALPHA_ARGS+=(--alpha "$a")
done

mkdir -p "$OUTDIR"

echo "[INFO] Running demuxlet sweep over geno-error-offsets: ${ERRS}"
for err in "${ERR_ARR[@]}"; do
  prefix="$OUTDIR/demuxlet_err${err}"
  echo "[INFO]  -> geno-error-offset=$err (output prefix: $prefix)"
  # shellcheck disable=SC2086
  docker run --rm -u "$(id -u)":"$(id -g)" \
    -v "$HOST_ROOT":"$HOST_ROOT" \
    --entrypoint /opt/conda/bin/popscle "$IMAGE" demuxlet \
      --plp "$PLP" \
      --vcf "$VCF" \
      --field GT \
      --group-list "$BARCODES" \
      --doublet-prior "$DOUBLETP" \
      --geno-error-coeff 0.0 \
      --geno-error-offset "$err" \
      --out "$prefix" \
      "${ALPHA_ARGS[@]}"
done

echo "[INFO] demuxlet sweep complete (results under $OUTDIR)"
