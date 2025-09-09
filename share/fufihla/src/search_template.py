import sys
import math
import re
import random
import pprint as pp
from itertools import combinations_with_replacement

########################################################################
# Read Input Lines and Populate Dictionary 'dat'
########################################################################

dat = {}

ds_regex = re.compile(r'ds:Z:(.*?)\n')
err_regex = re.compile(r'([+\-])([ACGTNacgtn\[\]]+)')

# Parse PAF input
for line in sys.stdin:
    if not line:
        break
    llst = line.split('\t')
    reads = llst[5]
    allele = llst[0]
    if not allele[-1].isdigit(): continue
    gene = allele.split('*')[0]
    allele_len = int(llst[1])
    ref_start = llst[7]
    ref_end   = llst[8]
    
    match_len = int(ref_end) - int(ref_start)
    if match_len <= 0:
        continue

    # Mismatchness (NM)
    nm_field = llst[12].split(':')
    if len(nm_field) < 3:
        continue
    nm = int(nm_field[2])



    ds_ret = ds_regex.search(line)
    if not ds_ret:
        continue

    ds_str = ds_ret.group(1)

    # Count raw mismatch symbols (* + -)
    mismatch = ds_str.count('*') + ds_str.count('+') + ds_str.count('-')

    # Ignore mismatch for bracketed events <=2bp -- possible sequencing error
    err_ret = err_regex.findall(ds_str)
    for sign, chunk in err_ret:
        stripped = re.sub(r'[\[\]]', '', chunk)
        if len(stripped) <= 2:
            mismatch -= 1

    mismatch_rate = round(mismatch * 1.0 / match_len, 7)
    coverage = match_len * 1.0 / allele_len

    # Only keep alignments with mismatch_rate < 0.03
    if mismatch_rate < 0.03:
        if reads not in dat:
            dat[reads] = []
        dat[reads].append([
            gene,          # [0]
            allele,        # [1]
            match_len,     # [2]
            mismatch,      # [3]
            nm,            # [4]
            coverage,      # [5]
            mismatch_rate, # [6]
            ref_start,     # [7]
            ref_end,       # [8]
            ds_str         # [9]
        ])

########################################################################
# Regex & Helper Functions
########################################################################

# Regex: we capture bracket characters in the insertion/deletion chunk,
# then strip them out for later.
ds_pattern = re.compile(r"""
    :(\d+)                                 # (1) match length
  | \*([ACGTNacgtn])([ACGTNacgtn])         # (2,3) single-base mismatch
  | \+([ACGTNacgtn\[\]]+)                  # (4) insertion chunk (includes brackets)
  | -([ACGTNacgtn\[\]]+)                  # (5) deletion chunk (includes brackets)
""", re.VERBOSE)

def parse_ds_string_full(ds_str, ref_start=0):
    """
    Stage-1 parse: build events for match/mismatch/indel.
    For insertions/deletions, store the entire chunk in 'all_bases'
    (including bracket chars). Will strip those out later.
    """
    events = []
    ref_pos = ref_start

    tokens = ds_pattern.findall(ds_str)
    # Each match => 5-tuple: (match_len, snp_from, snp_to, ins_chunk, del_chunk)

    for (match_len, snp_from, snp_to, ins_chunk, del_chunk) in tokens:

        # A) Match block
        if match_len:
            length = int(match_len)
            events.append({
                'type': 'match',
                'ref_start': ref_pos,
                'ref_end':   ref_pos + length
            })
            ref_pos += length
            continue

        # B) SNP => single-base mismatch
        if snp_from and snp_to:
            events.append({
                'type': 'mismatch',
                'ref_start': ref_pos,
                'ref_end':   ref_pos + 1
            })
            ref_pos += 1
            continue

        # C) Insertion => +(...) capturing group
        if ins_chunk:
            events.append({
                'type': 'insertion',
                'ref_start': ref_pos,
                'ref_end':   ref_pos,    # insertion doesn't consume index
                'all_bases': ins_chunk
            })
            continue

        # D) Deletion => -(...) capturing group
        if del_chunk:
            length = len(del_chunk)
            events.append({
                'type': 'deletion',
                'ref_start': ref_pos,
                'ref_end':   ref_pos + length,
                'all_bases': del_chunk
            })
            ref_pos += length
            continue

    return events

def subtract_mismatches_from_start(ds_str, length):
    """
    Subtract mismatches from the first `length` bases [0..length).
    Return (mismatch_sub, nm_sub).

    Logic:
      - For SNP => +1 mismatch each
      - For insertion/deletion => remove bracket chars from 'all_bases'. 
        Let stripped_len = len(stripped_bases). 
         * If stripped_len > 2 => +1 mismatch
         * Else => ignore
    """
    events = parse_ds_string_full(ds_str, ref_start=0)
    mismatch_sub = 0
    nm_sub = 0

    region_start = 0
    region_end   = length

    for e in events:
        e_start = e['ref_start']
        e_end   = e['ref_end']

        # Skip if event is fully before or fully beyond [region_start..region_end)
        if e_end <= region_start:
            continue
        if e_start >= region_end and e['type'] != 'insertion':
            break

        overlap_start = max(e_start, region_start)
        overlap_end   = min(e_end, region_end)
        overlap_len   = overlap_end - overlap_start

        if e['type'] == 'match':
            continue  # no mismatch

        elif e['type'] == 'mismatch':
            # Single-base => +1 if it overlaps
            mismatch_sub += 1
            nm_sub       += 1

        elif e['type'] == 'insertion':
            # Insertion => zero-length on reference
            # If e_start < region_end => it belongs
            if e_start <= region_end:
                # Strip bracket chars
                stripped = re.sub(r'[\[\]]', '', e['all_bases'])
                if len(stripped) > 2:
                    mismatch_sub += 1
                nm_sub       += len(stripped)
                # else => ignore

        elif e['type'] == 'deletion':
            # Partial overlap => fraction
            if overlap_len <= 0:
                continue
            event_span = e_end - e_start
            fraction   = overlap_len / float(event_span)

            # Remove brackets from 'all_bases'
            stripped = re.sub(r'[\[\]]', '', e['all_bases'])
            stripped_len = len(stripped)

            # If stripped_len >2 => +1 mismatch event
            if stripped_len > 2:
                mismatch_sub += 1
                # Increment nm_sub by overlap portion
                overlap_bases = int(round(stripped_len * fraction))
            nm_sub       += len(stripped)
            # else => ignore

    return mismatch_sub, nm_sub

def subtract_mismatches_from_end(ds_str, length):
    """
    Subtract mismatches from the last `length` bases of the alignment.
    Return (mismatch_sub, nm_sub).

    Same logic for partial overlap as 'subtract_mismatches_from_start'.
    """
    events = parse_ds_string_full(ds_str, ref_start=0)
    if not events:
        return 0, 0

    mismatch_sub = 0
    nm_sub = 0

    ref_alignment_length = max(e['ref_end'] for e in events)
    region_end   = ref_alignment_length
    region_start = max(0, ref_alignment_length - length)

    for e in events:
        e_start = e['ref_start']
        e_end   = e['ref_end']

        if e_end <= region_start:
            continue
        if e_start >= region_end and e['type'] != 'insertion':
            continue

        overlap_start = max(e_start, region_start)
        overlap_end   = min(e_end, region_end)
        overlap_len   = overlap_end - overlap_start

        if e['type'] == 'match':
            continue

        elif e['type'] == 'mismatch':
            mismatch_sub += 1
            nm_sub       += 1

        elif e['type'] == 'insertion':
            if e_start <= region_end:
                # Strip bracket chars
                stripped = re.sub(r'[\[\]]', '', e['all_bases'])
                if len(stripped) > 2:
                    mismatch_sub += 1
                nm_sub       += len(stripped)

        elif e['type'] == 'deletion':
            if overlap_len <= 0:
                continue
            event_span   = e_end - e_start
            fraction     = overlap_len / float(event_span)
            stripped     = re.sub(r'[\[\]]', '', e['all_bases'])
            stripped_len = len(stripped)

            if stripped_len > 2:
                mismatch_sub += 1
                overlap_bases = int(round(stripped_len * fraction))
            nm_sub       += len(stripped)

    return mismatch_sub, nm_sub

########################################################################
# Select Best Allele by Mismatch in the Overlapped Mapping Interval
########################################################################

weights = {}
best_for_read = {}

for read in dat:
    cand = dat[read]
    if not cand:
        continue

    allele_best = cand[0]
    best_set = {tuple(allele_best)}

    for allele in cand:
        if allele[1] == allele_best[1]:
            continue

        mismatch0   = allele_best[3]
        nm0         = allele_best[4]
        start_best  = int(allele_best[7])
        end_best    = int(allele_best[8])
        ds_str0     = allele_best[9]

        mismatch1   = allele[3]
        nm1         = allele[4]
        start_curr  = int(allele[7])
        end_curr    = int(allele[8])
        ds_str1     = allele[9]

        # Process unshared START
        if start_curr > start_best:
            change_start = start_curr - start_best
            sub_m0, sub_nm0 = subtract_mismatches_from_start(ds_str0, change_start)
            mismatch0 -= sub_m0
            nm0       -= sub_nm0
        elif start_curr < start_best:
            change_start = start_best - start_curr
            sub_m1, sub_nm1 = subtract_mismatches_from_start(ds_str1, change_start)
            mismatch1 -= sub_m1
            nm1       -= sub_nm1

        # Process unshared END
        if end_curr > end_best:
            change_end = end_curr - end_best
            sub_m1, sub_nm1 = subtract_mismatches_from_end(ds_str1, change_end)
            mismatch1 -= sub_m1
            nm1       -= sub_nm1
        elif end_curr < end_best:
            change_end = end_best - end_curr
            sub_m0, sub_nm0 = subtract_mismatches_from_end(ds_str0, change_end)
            mismatch0 -= sub_m0
            nm0       -= sub_nm0

        # Decide which allele is better
        if mismatch1 < mismatch0:
            allele_best = allele
            best_set = {tuple(allele_best)}
        elif mismatch1 == mismatch0:
            if nm1 < nm0:
                allele_best = allele
                best_set = {tuple(allele_best)}
            elif nm1 == nm0:
                best_set.add(tuple(allele))
            else:
                pass
        else:
            pass

    # Update weights
    if len(best_set) == 1:
        only_allele = list(best_set)[0]
        gene_best   = only_allele[0]
        allele_name = only_allele[1]
        coverage = only_allele[5]
        weights.setdefault(gene_best, {}).setdefault(allele_name, 0)
        weights[gene_best][allele_name] += coverage
    else:
        for tie_allele in best_set:
            gene_best   = tie_allele[0]
            allele_name = tie_allele[1]
            coverage = tie_allele[5]
            weights.setdefault(gene_best, {}).setdefault(allele_name, 0)
            weights[gene_best][allele_name] += coverage/len(best_set)

    best_for_read[read] = [(al[1], al[5], al[3], al[4]) for al in best_set]

########################################################################
# Write read->allele mapping
########################################################################

with open(sys.argv[1] + "_mapping.csv", "wt") as fp:
    for read, best_alleles in best_for_read.items():
        allele_names = [allele for allele, _, _, _ in best_alleles]
        allele_str = ",".join(allele_names) 
        fp.write(f"{read}\t{allele_str}\n")



########################################################################
# Pick top two alleles per gene -> "<prefix>_template_allelelist.txt"
########################################################################

from itertools import combinations_with_replacement

genes = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1']
template_alleles = []

def short_naming(allele_name):
    # Keep only the first 4 ':' fields, e.g., "HLA-A*02:01:01:01"
    return ":".join(allele_name.split(":")[:4])

with open(sys.argv[1] + "_template_allelelist.txt", "wt") as fp:
    for gene in genes:
        if gene not in weights:
            continue

        # Gather (allele_name, total_coverage), ignoring zero coverage
        arr = [(allele_name, vote_count) for allele_name, vote_count in weights[gene].items() if vote_count > 0]

        # Sort alleles by coverage descending
        arr.sort(key=lambda x: x[1], reverse=True)
        #print(arr, file=sys.stderr)
        if not arr:
            continue

        # Select top 15 alleles (or fewer if less than 15)
        top_alleles = [allele for allele, _ in arr[:15]]

        # Generate all possible allele pairs (a,a), (a,b), (b,b), but not (b,a)
        allele_pairs = list(combinations_with_replacement(top_alleles, 2))
        
        # Store coverage + mismatch + nm information for each allele pair
        pair_scores = {pair: {"coverage": 0, "mismatch": 0, "nm": 0, "count": 0} for pair in allele_pairs}

        # Assign read support to allele pairs based on best_for_read
        for read, supported_alleles in best_for_read.items():
            read_allele_coverage = {allele: coverage for allele, coverage, _, _ in supported_alleles}

            for pair in allele_pairs:
                a, b = pair
                if a in read_allele_coverage and b in read_allele_coverage:
                    # If the read supports both alleles, add the **average** of their coverages
                    avg_coverage = (read_allele_coverage[a] + read_allele_coverage[b]) / 2
                    pair_scores[pair]["coverage"] += avg_coverage
                    pair_scores[pair]["mismatch"] += sum(m for allele, _, m, _ in supported_alleles if allele in pair) / 2
                    pair_scores[pair]["nm"] += sum(n for allele, _, _, n in supported_alleles if allele in pair) / 2
                    pair_scores[pair]["count"] += 1
                elif a in read_allele_coverage:
                    pair_scores[pair]["coverage"] += read_allele_coverage[a]
                    pair_scores[pair]["mismatch"] += sum(m for allele, _, m, _ in supported_alleles if allele == a)
                    pair_scores[pair]["nm"] += sum(n for allele, _, _, n in supported_alleles if allele == a)
                    pair_scores[pair]["count"] += 1
                elif b in read_allele_coverage:
                    pair_scores[pair]["coverage"] += read_allele_coverage[b]
                    pair_scores[pair]["mismatch"] += sum(m for allele, _, m, _ in supported_alleles if allele == b)
                    pair_scores[pair]["nm"] += sum(n for allele, _, _, n in supported_alleles if allele == b)
                    pair_scores[pair]["count"] += 1

        # Compute a weighted score for each allele pair
        sorted_pairs = sorted(
            pair_scores.items(),
            key=lambda x: (
                x[1]["count"],
                x[1]["coverage"],          
                -x[1]["mismatch"],         
                -x[1]["nm"]     
            ),
            reverse=True
        )

        #print(sorted_pairs[:10], file=sys.stderr)  # Print top 10 pairs for debugging

        # Get the best pair based on the highest composite score
        best_pair = sorted_pairs[0][0]  # Extract allele pair from sorted list
        allele0, allele1 = best_pair

        # Append the best pair to the final selection
        template_alleles.append(allele0)
        if allele0 != allele1:  # Avoid adding duplicates if the best pair is (a, a)
            template_alleles.append(allele1)

    # Write them all to the template file
    for allele_name in template_alleles:
        fp.write(allele_name + "\n")

