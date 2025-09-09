#!/usr/bin/env python3
import sys

def filter_paf(file_input):
    """
    Filter PAF lines: keep those with (NM / alignment_length) <= 0.05 and NM <= 200.
    """
    for line in file_input:
        line = line.rstrip("\n")
        fields = line.split("\t")
        aln_start = int(fields[7])
        aln_end = int(fields[8])
        length = aln_end - aln_start
        nm_field = next(f for f in fields if f.startswith("NM:i:"))
        nm = int(nm_field.split(":", 2)[2])
        if nm <= 200 and (nm / length) <= 0.05:
            sys.stdout.write(line + "\n")

def main():
    if len(sys.argv) == 1:
        filter_paf(sys.stdin)
    elif len(sys.argv) == 2:
        with open(sys.argv[1], 'r') as f:
            filter_paf(f)
    else:
        prog = sys.argv[0]
        print(f"Usage: {prog} [input.paf]  (or pipe PAF into it)", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()

