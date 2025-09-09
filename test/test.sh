#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-HG002_out}"

mkdir -p "$OUTDIR"

/usr/bin/time -v fufihla HG002.fa.gz "$OUTDIR" \
  1>"$OUTDIR.out" 2>"$OUTDIR.err" || true

echo "[FuFiHLA test] Finished. Results in $OUTDIR/, $OUTDIR.out, logs: / $OUTDIR.err"

