#!/usr/bin/env python3
import sys
import re

###########################################################
# Call Closest Allele
###########################################################

def parse_gene_annot(file_path):
    """
    Parses the gene annotation file and returns a dict mapping allele
    names (e.g. "HLA-C*03:04:18") to a list of CDS intervals (1-based inclusive).
    """
    gene_data = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            info = {k: v for k, v in (fld.split('=',1) for fld in fields if '=' in fld)}
            if 'allele' in info and 'CDS' in info:
                allele = info['allele']
                cds_intervals = []
                for token in info['CDS'].split(','):
                    if '..' in token:
                        try:
                            s,e = token.split('..')
                            cds_intervals.append((int(s), int(e)))
                        except ValueError:
                            pass
                gene_data[allele] = cds_intervals
    return gene_data

def check_exon_coverage(cds_intervals, ref_start, ref_end, threshold=1.0):
    for s,e in cds_intervals:
        length = e - s + 1
        # Compute overlap between [interval_start,interval_end] and [ref_start,ref_end]
        ov = max(0, min(e, ref_end) - max(s, ref_start) + 1)
        if ov / length < threshold:
            return False
    return True

# Parse the cs tag
def parse_cs(cs, ref_offset=1):
    token_pattern = re.compile(r'(:\d+)|(\*[A-Za-z]{2})|([\+\-][A-Za-z]+)')
    pos = ref_offset
    muts = []
    for m in token_pattern.finditer(cs):
        tok = m.group(0)
        if tok.startswith(':'):
            pos += int(tok[1:])
        elif tok.startswith('*'):
            muts.append(('sub', pos, tok))
            pos += 1
        elif tok.startswith('+'):
            muts.append(('ins', pos, tok))
        elif tok.startswith('-'):
            length = len(tok) - 1
            muts.append(('del', pos, tok))
            pos += length
    return muts

# Penalty & PID calculation
PENALTY_CDS     = 10
PENALTY_NONCDS  = 1

def mutation_penalty(mutations, cds_intervals):
    pen = 0
    for _, pos, _ in mutations:
        in_cds = any(s <= pos <= e for s,e in cds_intervals)
        pen += (PENALTY_CDS if in_cds else PENALTY_NONCDS)
    return pen

def calculate_pid(cs):
    matches = total = 0
    token_pattern = re.compile(r'(:\d+)|(\*[A-Za-z]{2})|([\+\-][A-Za-z]+)')
    for m in token_pattern.finditer(cs):
        tok = m.group(0)
        if tok.startswith(':'):
            l = int(tok[1:])
            matches += l
            total   += l
        elif tok.startswith('*'):
            total   += 1
        elif tok.startswith('-'):
            total   += (len(tok) - 1)
        # Insertions are ignored in PID
    pid = (matches / total * 100) if total else 0.0
    return matches, total, pid

def find_csstr(line):
    m = re.search(r'cs:Z:(\S+)', line)
    return m.group(1) if m else ''

def main():
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: final_call_commonprefix.py <gene_annot.info>\n")
        sys.exit(1)
    gene_data = parse_gene_annot(sys.argv[1])
    da = {}

    for line in sys.stdin:
        line = line.rstrip()
        if not line: continue
        ll = line.split('\t')
        if len(ll) < 10: continue

        allele     = ll[5]
        rec_seq    = ll[0]
        try:
            ref_len  = float(ll[6])
            ref_s    = int(ll[7])
            ref_e    = int(ll[8])
            aln_len  = ref_e - ref_s


        except (ValueError, IndexError):
            continue


        # Skip short alignments
        if aln_len / ref_len < 0.5:
            continue

        # If we have CDS intervals, enforce exon coverage
        if allele in gene_data:
            cds = gene_data[allele]
            # Convert ref_start to 1-based
            if not check_exon_coverage(cds, ref_s + 1, ref_e, threshold=0.95):
                continue

        # Account for PAF’s 0-based target start:
        ref_s = int(ll[7])  
        ref_start_1based = ref_s + 1
        csstr    = find_csstr(line)
        mutations = parse_cs(csstr, ref_offset=ref_start_1based)

        matches, total_bases, pid = calculate_pid(csstr)
        # Compute penalty
        if allele in gene_data:
            pen = mutation_penalty(mutations, gene_data[allele])
        else:
            pen = csstr.count('*') + csstr.count('+') + csstr.count('-')

        da.setdefault(rec_seq, []).append([allele, aln_len, pen, pid, line])

    # Pick best per reconstructed seq
    for rec in da:
        # Sort by penalty ↑, then aln_len ↓
        best = sorted(da[rec], key=lambda x: (x[2], -x[1]))[0]
        allele, aln_len, pen, pid, raw = best
        print(f"{allele}\t{aln_len}\t{pen}\t{pid:.2f}\t@{raw}")

if __name__ == '__main__':
    main()

