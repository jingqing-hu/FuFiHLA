# FuFiHLA: HLA Typing Pipeline for Long Reads

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11+-blue.svg)]()
[![Conda](https://img.shields.io/badge/install%20with-conda-green.svg)](https://docs.conda.io/)

FuFiHLA is a modular pipeline for accurate HLA gene typing and consensus sequence construction from long-read sequencing data.  
It currently supports PacBio HiFi data on six clinically important transplant genes: **HLA-A, -B, -C, -DQA1, -DQB1, -DRB1**.

## Highlights
- **Reference-free**: does *not* depend on GRCh38 or other whole-genome assemblies  
- **Improved consensus accuracy** compared to [StarPhase](#) (link TBD)  

**Citation:** TBD

---

## Installation

Install from **Bioconda** (recommended):

```
conda install -c bioconda -c conda-forge fufihla
```

## Usage
```
fufihla <input_reads.fa.gz> <output_dir>
```
Arguments
- `<input_reads.fa.gz>` : raw PacBio HiFi reads (FASTA/FASTQ, gzipped)
- `<output_dir>` : directory for pipeline outputs

## Quick Test
A small toy dataset is included in `test/` for installation checks:
```
cd test
bash test.sh
```
This runs FuFiHLA on `HG002.fa.gz` and produces:
- `HG002_out/` → pipeline logs
- `HG002_out.out` → result output
- `HG002_out.err` → stderr log

## Outputs
A typical run produces:
```
<outdir>/consensus/*_asm*.fa        → consensus allele FASTA sequences
```
Allele calls are printed to standard output in **PAF-like format** with minimap2 tags.
Example:
```
HLA-A*01:01:01:01  HLA-A*01_01_01_01  ...  cs:Z::3503
HLA-A*26:01:01:01  HLA-A*26_01_01_01  ...  cs:Z::3517
```

- **Column 1** → the **allele name** called by FuFiHLA  
- **Last column** (`cs:Z`) → minimap2 cs tag encoding base-level matches/mismatches:
  - **Known Alleles**: `cs:Z::3503` → perfect match over 3503 bp 
  - **Novel Alleles** → cs:Z contains substitutions (*), insertions (+), or deletions (-)


