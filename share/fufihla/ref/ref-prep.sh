#!/usr/bin/env bash

# Resolve paths
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
echo "${SCRIPT_DIR}"
EASYIPD="${SCRIPT_DIR}/easyIPD.py"
GENELIST="${SCRIPT_DIR}/gene.list"

mkdir -p ref_data
cd ref_data

# Download IMGT/HLA database
wget -O hla.dat.zip "https://github.com/ANHIG/IMGTHLA/raw/refs/heads/Latest/hla.dat.zip"

# Extract and build
unzip -p hla.dat.zip > hla.dat
python3 "$EASYIPD" --elements < hla.dat > gene.csv
grep -f "$GENELIST" gene.csv | grep -v partial > gene.info
cut -f4 gene.info | cut -d'=' -f2 > sel.allele.list
python3 "$EASYIPD" --geneSeq --alleleF sel.allele.list < hla.dat | gzip -c > ref.gene.fa.gz

echo "Done: ref.gene.fa.gz (and gene.list) in $(pwd)"

