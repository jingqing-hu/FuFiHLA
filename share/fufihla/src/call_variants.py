#!/usr/bin/env python3
import sys, os, re, gzip
import glob, subprocess
from collections import defaultdict

if len(sys.argv) != 6:
    sys.stderr.write(
        "Usage: call_variants_bcf.py "
        "<template_list> <filtered_paf.gz> "
        "<reads_fa.gz> <all_allele_seq.fa.gz> <outdir>\n"
    )
    sys.exit(1)

template_file    = sys.argv[1]
filtered_paf_gz  = sys.argv[2]
reads_fa         = sys.argv[3]
allele_refs_fa   = sys.argv[4]
OUTDIR           = sys.argv[5]

SAMPLE = os.path.basename(OUTDIR)
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
KEEP_MAJOR = os.path.join(SCRIPT_DIR, "keep_major.py")




ds_regex = re.compile(r'ds:Z:(.*)')
err_regex = re.compile(r'([+\-])([ACGTNacgtn\[\]]+)')
###########################################################
# Load Template Alleles from sys.argv[1]
###########################################################
def load_templates(template_file):
    template_alleles = {}
    with open(template_file, 'rt') as fp:
        for line in fp:
            allele = line.strip()
            if not allele:
                continue
            gene = allele.split('*')[0]
            template_alleles.setdefault(gene, []).append(allele)
    return template_alleles

###########################################################
# Candidate Alignment Processing from s4_filtered.paf.gz
###########################################################
# Read candidate alignments from sys.argv[2] (a gzip file).
# Extract alignment information from template allele 

dat = {}
template_alleles = load_templates(sys.argv[1])
with gzip.open(sys.argv[2], 'rt') as fp:
    for line in fp:
        line = line.strip()
        if not line:
            continue
        llst = line.split('\t')
        # In s4_filtered.paf.gz:
        # Field[0]: candidate allele (template allele)
        # Field[1]: allele length.
        # Field[5]: read identifier.
        # Field[7]: ref_start.
        # Field[8]: ref_end.
        read = llst[5]
        allele = llst[0]
        gene = allele.split('*')[0]
        # Only keep candidate alignments if the allele is in the template set.
        if allele not in template_alleles.get(gene, []):
            continue
        allele_len = int(llst[1])
        ref_start = llst[7]
        ref_end = llst[8]
        match_len = int(ref_end) - int(ref_start)
        if match_len <= 0:
            continue
        nm_field = llst[12].split(':')
        if len(nm_field) < 3:
            continue
        nm = int(nm_field[2])
        m = ds_regex.search(line)
        if not m:
            continue
        ds_str = m.group(1)
        # Count raw mismatch symbols.
        mismatch = ds_str.count('*') + ds_str.count('+') + ds_str.count('-')
        # Ignore mismatch for bracketed events less than or equal to 2 bp -- possible sequencing error
        for sign, chunk in err_regex.findall(ds_str):
            stripped = re.sub(r'[\[\]]', '', chunk)
            if len(stripped) <= 2:
                mismatch -= 1
        mismatch_rate = round(mismatch * 1.0 / match_len, 7)
        coverage = match_len * 1.0 / allele_len
        # Only keep candidate alignments with mismatch_rate < 0.03.
        if mismatch_rate < 0.03:
            dat.setdefault(read, []).append([
                gene,          # [0]
                allele,        # [1]
                match_len,     # [2]
                mismatch,      # [3]
                nm,            # [4]
                coverage,      # [5]
                mismatch_rate, # [6]
                ref_start,     # [7] (string)
                ref_end,       # [8] (string)
                ds_str         # [9]
            ])

###########################################################
# Read Assignment by Mismatch in the Overlapped Mapping Interval
###########################################################
def subtract_mismatches_from_start(ds_str, length):
    events = parse_ds_string_full(ds_str, ref_start=0)
    mismatch_sub = 0
    nm_sub = 0
    region_start = 0
    region_end = length
    for e in events:
        e_start = e['ref_start']
        e_end = e['ref_end']
        if e_end <= region_start:
            continue
        if e_start >= region_end and e['type'] != 'insertion':
            break
        overlap_len = min(e_end, region_end) - max(e_start, region_start)
        if e['type'] == 'mismatch':
            mismatch_sub += 1
            nm_sub += 1
        elif e['type'] == 'insertion':
            if e_start < region_end:
                stripped = re.sub(r'[\[\]]', '', e['all_bases'])
                if len(stripped) > 0:
                    mismatch_sub += 1
                nm_sub += len(stripped)
        elif e['type'] == 'deletion':
            if overlap_len > 0:
                stripped = re.sub(r'[\[\]]', '', e['all_bases'])
                if len(stripped) > 0:
                    mismatch_sub += 1
                nm_sub += len(stripped)
    return mismatch_sub, nm_sub

def subtract_mismatches_from_end(ds_str, length):
    events = parse_ds_string_full(ds_str, ref_start=0)
    if not events:
        return 0, 0
    mismatch_sub = 0
    nm_sub = 0
    ref_alignment_length = max(e['ref_end'] for e in events)
    region_end = ref_alignment_length
    region_start = max(0, ref_alignment_length - length)
    for e in events:
        e_start = e['ref_start']
        e_end = e['ref_end']
        if e_end <= region_start:
            continue
        if e_start >= region_end and e['type'] != 'insertion':
            continue
        overlap_len = min(e_end, region_end) - max(e_start, region_start)
        if e['type'] == 'mismatch':
            mismatch_sub += 1
            nm_sub += 1
        elif e['type'] == 'insertion':
            if e_start < region_end:
                stripped = re.sub(r'[\[\]]', '', e['all_bases'])
                if len(stripped) > 0:
                    mismatch_sub += 1
                nm_sub += len(stripped)
        elif e['type'] == 'deletion':
            if overlap_len > 0:
                stripped = re.sub(r'[\[\]]', '', e['all_bases'])
                if len(stripped) > 0:
                    mismatch_sub += 1
                nm_sub += len(stripped)
    return mismatch_sub, nm_sub

ds_pattern = re.compile(r"""
    :(\d+)
  | \*([ACGTNacgtn])([ACGTNacgtn])
  | \+([ACGTNacgtn\[\]]+)
  | -([ACGTNacgtn\[\]]+)
""", re.VERBOSE)

def parse_ds_string_full(ds_str, ref_start=0):
    events = []
    ref_pos = ref_start
    tokens = ds_pattern.findall(ds_str)
    for (match_len, snp_from, snp_to, ins_chunk, del_chunk) in tokens:
        if match_len:
            length = int(match_len)
            events.append({
                'type': 'match',
                'ref_start': ref_pos,
                'ref_end': ref_pos + length
            })
            ref_pos += length
            continue
        if snp_from and snp_to:
            events.append({
                'type': 'mismatch',
                'ref_start': ref_pos,
                'ref_end': ref_pos + 1
            })
            ref_pos += 1
            continue
        if ins_chunk:
            events.append({
                'type': 'insertion',
                'ref_start': ref_pos,
                'ref_end': ref_pos,
                'all_bases': ins_chunk
            })
            continue
        if del_chunk:
            length = len(del_chunk)
            events.append({
                'type': 'deletion',
                'ref_start': ref_pos,
                'ref_end': ref_pos + length,
                'all_bases': del_chunk
            })
            ref_pos += length
            continue
    return events

def csstr_to_var(csstr, fro, to) :
    cslst = re.split(r'(:|\*|\+|\-|\~)', csstr)[1:]
    cslst_L = int(len(cslst)/2)
    cur_pos = fro
    out = [fro, to]
    for i in range(cslst_L) :
        tag = cslst[i*2]
        signal = cslst[i*2+1]
        if tag == ':' :
            cur_pos += int(signal)
        elif tag == '*' :
            out.append(str(cur_pos)+':'+'*'+ signal)
            cur_pos += 1
        elif tag == '+' :
            out.append(str(cur_pos)+':'+'+'+ signal)
        elif tag == '-' :
            out.append(str(cur_pos)+':'+'-'+ signal)
            cur_pos += len(signal)
    return(out)


# Now assign read to allele based on mismatch in the overlapped mapping interval.
best_for_read = {}
for read, cand_list in dat.items():
    if not cand_list:
        continue
    best_candidate = cand_list[0]
    best_set = {tuple(best_candidate)}
    for candidate in cand_list[1:]:
        # Extract values for best_candidate.
        mismatch0 = best_candidate[3]
        nm0 = best_candidate[4]
        start0 = int(best_candidate[7])
        end0 = int(best_candidate[8])
        ds_str0 = best_candidate[9]
        # Extract candidate's values.
        mismatch1 = candidate[3]
        nm1 = candidate[4]
        start1 = int(candidate[7])
        end1 = int(candidate[8])
        ds_str1 = candidate[9]
        
        # Adjust differences at the start.
        if start1 > start0:
            diff = start1 - start0
            sub0, subnm0 = subtract_mismatches_from_start(ds_str0, diff)
            mismatch0_adj = mismatch0 - sub0
            nm0_adj = nm0 - subnm0
            mismatch1_adj = mismatch1
            nm1_adj = nm1
        elif start1 < start0:
            diff = start0 - start1
            sub1, subnm1 = subtract_mismatches_from_start(ds_str1, diff)
            mismatch1_adj = mismatch1 - sub1
            nm1_adj = nm1 - subnm1
            mismatch0_adj = mismatch0
            nm0_adj = nm0
        else:
            mismatch0_adj = mismatch0
            nm0_adj = nm0
            mismatch1_adj = mismatch1
            nm1_adj = nm1
        
        # Adjust differences at the end.
        if end1 > end0:
            diff = end1 - end0
            sub1, subnm1 = subtract_mismatches_from_end(ds_str1, diff)
            mismatch1_adj -= sub1
            nm1_adj -= subnm1
        elif end1 < end0:
            diff = end0 - end1
            sub0, subnm0 = subtract_mismatches_from_end(ds_str0, diff)
            mismatch0_adj -= sub0
            nm0_adj -= subnm0

        # Select candidate with lower adjusted mismatches.
        if mismatch1_adj < mismatch0_adj:
            best_candidate = candidate
        elif mismatch1_adj == mismatch0_adj:
            if nm1_adj < nm0_adj:
                best_candidate = candidate
            # Else, if equal, keep current.
    best_for_read[read] = best_candidate

with open(os.path.join(OUTDIR,"read_assignment.txt"),"w") as f:    
    for key, value in best_for_read.items():
        f.write(f"{key}\t{value[1]}\n")

allele_to_reads = defaultdict(list)

for read_id, best_candidate in best_for_read.items():
    allele = best_candidate[1]
    allele_to_reads[allele].append(read_id)

with open(os.path.join(OUTDIR,"allele_to_reads.txt"), "w") as out:
    for allele, reads in allele_to_reads.items():
        out.write(f"{allele}\t{','.join(reads)}\n")


###########################################################
# Consensus Sequence Reconstruction based on Templates
###########################################################
INPUT_READS_FA = reads_fa
REF_ALLELES_FA = allele_refs_fa
THREADS          = 4
MM2OPTS          = ["-ax","asm10","--ds",f"-t{THREADS}"]
GENES            = ["HLA-A","HLA-B","HLA-C","HLA-DQA1","HLA-DQB1","HLA-DRB1"]

subdirs = ["lists","fas","bam","vcf","consensus","paf","vcf_log"]
for sub in subdirs:
    os.makedirs(os.path.join(OUTDIR, sub), exist_ok=True)


# Detect and force homozygosity when one allele has ≥4× more reads
gene_alleles = {}
for gene in GENES:
    # Find all called alleles for this gene
    alleles = [a for a in allele_to_reads if a.startswith(gene + "*")]
    if len(alleles) == 2:
        r0 = allele_to_reads[alleles[0]]
        r1 = allele_to_reads[alleles[1]]
        if len(r0) >= 4 * len(r1):
            # Allele[0] is dominant → homozygous
            gene_alleles[gene] = [alleles[0], alleles[0]]
        elif len(r1) >= 4 * len(r0):
            # Allele[1] is dominant → homozygous
            gene_alleles[gene] = [alleles[1], alleles[1]]
        else:
            # Balanced → keep as-is
            gene_alleles[gene] = alleles
    else:
        # Unexpected count (e.g. 0 or >2): just pass through
        gene_alleles[gene] = alleles


new_allele_list = []
known_allele_list = []

for gene in GENES:
    alleles = gene_alleles.get(gene, [])
    if len(alleles) != 2:
        print(f"Warning: expected 2 alleles for {gene}, found {len(alleles)} → {alleles}", file=sys.stderr)

    for idx, allele in enumerate(alleles, start=1):
        allele_safe = allele.replace(":", "_")
        list_file = os.path.join(OUTDIR,    "lists",      f"{allele_safe}.list")
        reads_fa_out = os.path.join(OUTDIR, "fas", f"{allele_safe}_reads.fa")
        ref_normal   = os.path.join(OUTDIR, "fas", f"{allele}.fa")
        ref_fa    = os.path.join(OUTDIR, "fas", f"{allele_safe}.fa")
        bam       = os.path.join(OUTDIR,    "bam",        f"{gene}_asm{idx}.bam")
        vcf       = os.path.join(OUTDIR,    "vcf",        f"{gene}_asm{idx}.vcf.gz")
        cons_fa   = os.path.join(OUTDIR,    "consensus",  f"{gene}_asm{idx}.fa")

        # 1) Write the list of reads
        with open(list_file, "w") as lf:
            lf.write("\n".join(allele_to_reads[allele]) + "\n")
        
        # 2) Extract reads FASTA
        reads_fa_out = os.path.join(OUTDIR, "fas", f"{gene}_{allele_safe}_reads.fa")
        with open(reads_fa_out, "w") as rf:
            subprocess.run(
                ["seqtk", "subseq", INPUT_READS_FA, list_file],
                stdout=rf, check=True
            )

        # 3) Extract reference allele sequence via seqtk:
        # 3a) Write a one‑line “ref_list” of the sanitized (":" to "_") allele name
        ref_list = os.path.join(OUTDIR, "lists", f"{gene}_{allele_safe}_ref.list")
        with open(ref_list, "w") as rl:
            rl.write(allele_safe + "\n")


        # 3b) Decompress the big allele FASTA, sanitize headers, and subseq
        with open(ref_fa, "w") as rf2:
            cmd = (
                f"gzip -dc {allele_refs_fa} | "
                f"sed 's/:/_/g' | "
                f"seqtk subseq - {ref_list}"
            )
            subprocess.run(cmd, shell=True, stdout=rf2, check=True)


        with open(ref_normal, "w") as rfn:
            subprocess.run(
                f"sed '1s/_/:/g' {ref_fa}",
                shell=True,
                check=True,
                stdout=rfn
            )


        # 3c) Index the per‑allele FASTA for downstream tools
        subprocess.run(["samtools", "faidx", ref_fa], check=True)       

        # 4) Map reads → allele, sort & index
        p1 = subprocess.Popen(
            ["minimap2"] + MM2OPTS + [ref_fa, reads_fa_out],
            stdout=subprocess.PIPE
        )
        p2 = subprocess.Popen(
            ["samtools","sort", "-o", bam, f"-@{THREADS}", "-m2g"],
            stdin=p1.stdout
        )
        p1.stdout.close()
        p2.communicate()
        subprocess.run(["samtools","index", bam], check=True)

        # 5) Call variants & index using longcallD by Yan Gao 
        # "https://github.com/yangao07/longcallD"
        
        cl = subprocess.Popen(
            ["longcallD", "call", ref_fa, bam],
            stdout=subprocess.PIPE
        )

        # write raw VCF for log purposes
        vcf_log = os.path.join(OUTDIR, "vcf_log", f"{gene}_asm{idx}.vcf")
        tee = subprocess.Popen(
            ["tee", vcf_log],
            stdin=cl.stdout,
            stdout=subprocess.PIPE
        )
        cl.stdout.close()

        km = subprocess.Popen(
            ["python3", KEEP_MAJOR],
            stdin=tee.stdout,
            stdout=subprocess.PIPE
        )
        tee.stdout.close()

        # compress & write final VCF
        with open(vcf, "wb") as vf:
            subprocess.run(
                ["bgzip","-c"],
                stdin=km.stdout,
                stdout=vf,
                check=True
            )
        km.stdout.close()

        # index it
        subprocess.run(["bcftools","index","-t", vcf], check=True)

        # 6a) Check whether the vcf is empty; if so, skip consensus/mapping
        with gzip.open(vcf, "rt") as vaf:
            has_variant = any(line.strip() and not line.startswith('#') for line in vaf)
        if has_variant:
            new_allele_list.append(cons_fa)
        else:
            known_allele_list.append(cons_fa) 

        # 6b) Build consensus
        with open(cons_fa, "w") as cf:
            subprocess.run(
                ["bcftools","consensus","-f", ref_fa, vcf],
                stdout=cf,
                check=True
            )

          
with open(os.path.join(OUTDIR, "new_allele.fa"), "w") as out:
    for fname in new_allele_list:
        with open(fname, "r") as infile:
            out.write(infile.read())


def normalize_header(header_line: str) -> str:
    # Take the first token (">HLA-A*02_06_01_01"), change back to normal nomenclature
    token = header_line.strip().split()[0]      
    allele = token.lstrip('>')                    
    return allele.replace('_', ':') # → "HLA-A*02:06:01:01"

with open(os.path.join(OUTDIR, "known_allele.paf"), "w") as out_paf:
    for cons_fa in known_allele_list:
        # Read the header from the consensus file to get the allele name
        with open(cons_fa) as fh:
            header = next(l for l in fh if l.startswith(">"))
        allele_name = normalize_header(header)
        allele_fa = os.path.join(OUTDIR, "fas", f"{allele_name}.fa")

        # Map assembled consensus (query) to the known allele (reference)
        subprocess.run(
            ["minimap2", "-c", "--cs", allele_fa, cons_fa],
            stdout=out_paf,
            check=True
        )
