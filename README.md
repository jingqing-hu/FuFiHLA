# FuFiHLA: HLA Typing Pipeline for Long Reads

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11+-blue.svg)]()

FuFiHLA is a pipeline for accurate HLA gene full field typing and consensus sequence construction from long-read sequencing data.  
It currently supports PacBio HiFi data on six clinically important transplant genes: **HLA-A, -B, -C, -DQA1, -DQB1, -DRB1**.

## Highlights
- **Reference-free**: does *not* depend on a specific version of reference genome such GRCh38 or CHM13  
- **Improved consensus accuracy** compared to [StarPhase](#) (link TBD)  

**Citation:** TBD

---

## Installation

Install from **Bioconda** (recommended):

```
conda install -c bioconda -c conda-forge fufihla
```

## Usage
To download the latest reference allele sequences `ref.gene.fa.gz` from IMGT:
```
fufihla-ref-prep
```
To run the pipeline:
```
fufihla --fa <input_reads.fa.gz> --out <output_dir>
fufihla --fa <input_reads.fa.gz> --refdir <reference data directory> --out <output_dir> --debug
```
Arguments
- `<input_reads.fa.gz>` : raw PacBio HiFi reads (FASTA/FASTQ, gzipped)
- `<output_dir>` : directory for pipeline outputs
- `--refdir <reference_data_directory>`(optional): path to reference allele dataset; if omitted, uses the default bundled set
- `--debug`(optional): keep all intermediate files; otherwise only consensus results are kept

## Quick Test
A small toy dataset is included in `test/` for installation checks:
```
<path to installation>/test/HG002.fa.gz
```
You can run with: 
```
fufihla --fa <path to installation>/test/HG002.fa.gz --out HG002
```
The output includes:
- `HG002/` → pipeline logs
- `HG002.out` → result output
- `HG002.err` → stderr log

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

## Running tips

1) extract reads from exist bam files can also generate similar result as using WGS reads.

```bash
echo "
chr6	29942254	29945755
chr6	31268254	31272571
chr6	31353362	31357442
chr6	32578769	32589848
chr6	32636717	32643200
chr6	32660031	32667132" > sel.bed

samtools view -bh ${bam} --region-file sel.bed | samtools fasta | gzip -c > out.fa.gz

```

