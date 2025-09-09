#!/usr/bin/env python3
import os
import sys
import argparse
import gzip
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("raw_reads_id")
    parser.add_argument("s2_inp_paf")
    parser.add_argument("target_genes_file")
    parser.add_argument("out_prefix")
    args = parser.parse_args()    
    with open(args.target_genes_file, 'r') as f:
        target_genes = set(line.strip() for line in f if line.strip())

    with open(args.raw_reads_id, 'r') as f:
        read_ids_of_interest = set(line.strip() for line in f if line.strip())

    alignments = defaultdict(list)
    # Parse PAF
    with gzip.open(args.s2_inp_paf, 'rt') as f:
        for line in f:
            fields = line.strip().split("\t")
            read_id = fields[0]
            ref_len = int(fields[6])
            ref_start = int(fields[7])
            ref_end = int(fields[8])
            match_len = ref_end - ref_start
            coverage = match_len / ref_len
            allele = fields[5]
            gene = allele.split('*')[0]
            if coverage < 0.4 and ref_start > 200 and ref_end < ref_len - 200:
                continue
            nm = int(fields[12].split(":")[2])
            # Score = mismatches/matching length
            score = nm / match_len
            alignments[read_id].append({
                'gene': gene,
                'ref_length': ref_len,
                'match_length': match_len,
                'nm': nm,
                'score': score,
                'line': line.strip()
            })

    # Prepare a per-gene dictionary for read-IDs
    gene_read_ids = defaultdict(set)
    gene_list_dir = f"{args.out_prefix}_list"
    os.makedirs(gene_list_dir, exist_ok=True)

    # For each read ID, select the best alignment (lowest score)
    # Then, if any alignment's score is within 10x of the best score, print the read ID and write the best alignment line.
    for read_id in read_ids_of_interest:
        if read_id not in alignments:
            continue
        best_score = min(aln['score'] for aln in alignments[read_id])
        valid_alignments = [aln for aln in alignments[read_id] if aln['score'] <= 10 * best_score and aln['gene'] in target_genes]
        if valid_alignments:
            print(read_id)  # Print read ID to s2_out.read-id.txt
            for aln in valid_alignments:
                # Record this read_id under that gene
                gene_read_ids[aln['gene']].add(read_id)

    for gene, ids in gene_read_ids.items():
        fn = os.path.join(gene_list_dir, f"{gene}.list")
        with open(fn, 'w') as gf:
            for rid in sorted(ids):
                gf.write(rid + "\n")

if __name__ == "__main__":
    main()
