#!/usr/bin/env bash
set -Eeuo pipefail

install -d "$PREFIX/bin"
install -m 0755 "bin/fufihla" "$PREFIX/bin/"

install -d "$PREFIX/share/fufihla"
cp -r "share/fufihla/"* "$PREFIX/share/fufihla/"

chmod 0755 "$PREFIX/share/fufihla/scripts/FuFiHLA.sh" || true

