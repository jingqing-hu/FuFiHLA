# FuFiHLA: Full Field HLA allele typing for Long Reads

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11+-blue.svg)]()

FuFiHLA is a pipeline for full-field HLA allele typing and consensus sequence construction from long-read sequencing data.
It performs complete allele typing and can additionally report partial (CDS-supported) alleles when mutations are detected within coding exons.

It currently supports PacBio HiFi data on six clinically important transplant genes: **HLA-A, -B, -C, -DQA1, -DQB1, -DRB1**.

## Highlights

* **Reference-free**: does *not* depend on a specific version of reference genome such GRCh38 or CHM13
* **Improved consensus accuracy** compared to StarPhase

**Citation:** Jingqing Hu, Qian Qin, Heng Li, Ying Zhou, FuFiHLA: A tool for Full-Field HLA typing from long-read data, Bioinformatics, 2026;, btag231, https://doi.org/10.1093/bioinformatics/btag231

---

## Installation

Install from **Bioconda** (recommended):

```
conda install -c bioconda -c conda-forge fufihla
```

---

## Quick Test

With `test.fa.gz` under the folder `test`, run:

```
fufihla --fa test.fa.gz --out test_dir
```

The output includes:

* `test_dir/` → pipeline logs
* `test_dir.out` → result output
* `test_dir.err` → stderr log


---

## Usage

To use the latest reference allele sequences from IMGT, type:

```
fufihla-ref-prep
```

This will create a directory called `ref_data`, which would contain the reference allele sequence `ref.gene.fa.gz`.

Run the pipeline:

```
# default bundled reference (IPD-IMGT/HLA)
fufihla --fa <input_reads.fa.gz> --out <output_dir>

# specify reference directory and options
fufihla --fa <input_reads.fa.gz> --out <output_dir> --refdir <reference data directory> --hifi/--ont --hom-threshold <float> --debug
```

### Arguments

* `<input_reads.fa.gz>` : raw long reads (.fa/.fa.gz/.fq/.fq.gz)
* `<output_dir>` : directory for pipeline outputs
* `--refdir <reference_data_directory>` (optional): path to reference allele dataset; if omitted, uses the default bundled set
* `--hifi/--ont` (optional): choose HiFi or Nanopore input reads, default is `--hifi`
* `--hom-threshold <float>` (optional): threshold for homozygous detection.
  When two haplotypes are highly similar, FuFiHLA suppresses a second allele call to avoid false heterozygous typing.
* `--debug` (optional): keep all intermediate files; otherwise only consensus results are retained

---

## Outputs

A typical run produces:

```
<outdir>/consensus/*_asm*.fa
```

Consensus allele FASTA sequences for each gene haplotype.

Allele calls are printed to `<output_dir>.out` in **PAF-like format** with minimap2 tags.

Example:

```
HLA-A*01:01:01:01  cons_HLA-A*01_01_01_01  ...  cs:Z::3503
HLA-A*26:01:01:01  cons_HLA-A*26_01_01_01  ...  cs:Z::3517
```

* **Column 1** → final allele call
* **Column 2** → consensus sequence identifier
* **Last column (`cs:Z`)** → minimap2 cs tag describing base-level alignment

Interpretation:

* **Known alleles**: `cs:Z::3503` → perfect match across the gene
* **Novel alleles**: substitutions `*`, insertions `+`, or deletions `-` appear in `cs:Z`

When partial alleles are reported, the allele originates from CDS alignment rather than full gene alignment.
This occurs only when exon mutations prevent confident full-length allele assignment.

---

## Running Tips

Extract reads from existing BAM files can also generate similar results:

```
echo "
chr6	29942254	29945755
chr6	31268254	31272571
chr6	31353362	31357442
chr6	32578769	32589848
chr6	32636717	32643200
chr6	32660031	32667132" > sel.bed

samtools view -bh ${bam} --region-file sel.bed | samtools fasta | gzip -c > out.fa.gz
```
