# FuFiHLA: HLA Typing Pipeline for Long Reads

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11+-blue.svg)]()
[![Conda](https://img.shields.io/badge/install%20with-conda-green.svg)](https://docs.conda.io/)

FuFiHLA is a modular pipeline for accurate HLA gene typing and consensus sequence construction from long-read sequencing data.
Currently, it only supports PacBio HiFi data on six graft transplant genes (HLA-A, -B, -C, -DQA1, -DQB1, -DRB1).

Highlights:
  * No dependent on reference genome (eg. GRCh38)
  * More accurate than starphase[add link here] in consensus sequence construction.

Citation: TBD

---

## Installation

Clone the repository and use conda for installation

```bash
git clone git@github.com:jingqing-hu/FuFiHLA.git
cd FuFiHLA
conda env create -f env/environment.yml
conda activate fufihla
```


## Usage

```bash
bash FuFiHLA.sh <input_reads.fa.gz> <output_dir> [gene_list.txt]
```
Arguments:
<input_reads.fq.gz>: raw PacBio HiFi reads (FASTA/FASTQ, gzipped)
<output_dir>: directory for outputs
[gene_list.txt] (optional): file listing HLA genes to analyze, one gene each line (currently only supporting six genes). Example:
```bash
HLA-A
HLA-B
HLA-C
HLA-DRB1
HLA-DQB1
HLA-DQA1
```
If no gene list is provided, the default (dps-dat/6genes.txt) is used.

We provide a tiny toy dataset in test/ for installation checks:
```bash
cd test/
bash test.sh
```

## Outputs
A typical run produces:
``` bash
<outdir>/consensus/   → consensus allele sequences (FASTA)
<outdir>/new_allele.fa → novel consensus alleles
<outdir>/known_allele.paf → mappings to known alleles
```
Allele calls are printed in standard output with alignment information
