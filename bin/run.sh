#!/usr/bin/env bash
set -euo pipefail

# Simple helper to run the pipeline using the provided profile
SAMPLES=${1:-samples.csv}
PROFILE=${2:-standard}

nextflow run main.nf --samples ${SAMPLES} -profile ${PROFILE}
