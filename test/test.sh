#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-HG002}"

mkdir -p "$OUTDIR"

/usr/bin/time -v fufihla --fa HG002.fa.gz --out HG002_out --refdir ../share/fufihla/dps-dat/ --debug 1>$OUTDIR.out 2>$OUTDIR.err

echo "[FuFiHLA test] Finished. Results in $OUTDIR/, $OUTDIR.out, logs: / $OUTDIR.err"

