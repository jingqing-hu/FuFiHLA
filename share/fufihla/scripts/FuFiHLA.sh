#!/usr/bin/env bash
set -euo pipefail

# Directory of this script
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && cd .. && pwd -P)"

# Usage check
if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <input_reads.fq.gz> <output_dir> [<gene_list_file|CSV>]"
  echo "  - If <gene_list> is omitted, default is: \${SCRIPT_DIR}/dps-dat/6genes.txt"
  echo "  - <gene_list> must be:"
  echo "      * a file with one gene per line (e.g. HLA-A), OR"
  echo "      * a CSV string (e.g. HLA-A,HLA-B,HLA-C)"
  exit 1
fi

INP_READS="$1"
OUT_DIR="$2"
USER_GENE_ARG="${3:-}"   # optional
mkdir -p "${OUT_DIR}"

DEFAULT_GENE_LIST="${SCRIPT_DIR}/dps-dat/6genes.txt"
GENE_LIST=""

if [[ -z "${USER_GENE_ARG}" ]]; then
  # Use default
  GENE_LIST="${DEFAULT_GENE_LIST}"
elif [[ -f "${USER_GENE_ARG}" ]]; then
  # Use user-provided file
  GENE_LIST="$(cd "$(dirname "${USER_GENE_ARG}")" && pwd -P)/$(basename "${USER_GENE_ARG}")"
else
  # User gave CSV → write file
  GENE_LIST="${OUT_DIR}/genes.list.txt"
  IFS=',' read -r -a _genes_csv <<< "${USER_GENE_ARG}"
  printf "%s\n" "${_genes_csv[@]}" \
    | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' \
    | awk 'NF' \
    | awk '!seen[$0]++' > "${GENE_LIST}"
fi

if ! awk 'BEGIN{ok=1} {if ($0 !~ /^HLA-[A-Za-z0-9]+$/) {ok=0; exit}} END{exit ok?0:1}' "${GENE_LIST}"; then
  echo "ERROR: Gene list must have one gene per line, format 'HLA-<NAME>' (e.g., HLA-A)" >&2
  echo "Invalid content found in: ${GENE_LIST}" >&2
  exit 1
fi

echo "[info] Using gene list: ${GENE_LIST}"
echo "${GENE_LIST}" | sed 's/^/[info]   /'

export GENE_LIST


# Define which steps to run (true/false)
declare -A STEP=(
  [0]=true
  [1]=true
  [2]=true
  [3]=true
  [4]=true
  [5]=true
  [6]=true
  [7]=true
)

log_time() {
  local STEP_NUM="$1"
  >&2 echo "Step ${STEP_NUM} completed at: $(date '+%Y-%m-%d %H:%M:%S')"
}

step0() {
  echo "Step 0: Verify input reads file"
  [[ -e "${INP_READS}" ]] || { echo "Error: ${INP_READS} not found"; exit 1; }
  [[ -r "${INP_READS}" ]] || { echo "Error: ${INP_READS} not readable"; exit 1; }
}

step1() {
  echo "Step 1: Raw alignment"
  local OUT_PAF="${OUT_DIR}/s1_raw_reads.paf.gz"
  local THREADS=20
  local TEMPLATE="${SCRIPT_DIR}/dps-dat/136ref.fa.gz"
  local OPTS="-t ${THREADS} -e10000 -k19 -w19 -g10k -A1 -B4 -O6,26 -E2,1 -s350 -c --end-bonus=10"
  minimap2 ${OPTS} "${TEMPLATE}" "${INP_READS}" | gzip -c > "${OUT_PAF}"
}

step2() {
  echo "Step 2: Gain informative reads from step-1 alignments"
  local s2_inp_reads="${INP_READS}"
  local s1_out_paf="${OUT_DIR}/s1_raw_reads.paf.gz"
  local s2_out_pref="${OUT_DIR}/s2_out"
  local target_genes_file="${GENE_LIST}"
  local raw_reads_id="${s2_out_pref}.raw_read-id.txt"
  local output_reads_id="${s2_out_pref}.read-id.txt"
  local output_reads_fa="${s2_out_pref}.read-fa.gz"
  local s2_filter="${SCRIPT_DIR}/src/raw_reads_filter.py"

  mkdir -p "$(dirname "${raw_reads_id}")"
  zgrep -f "${target_genes_file}" "${s1_out_paf}" \
    | cut -f1 | sort -u > "${raw_reads_id}"

  python3 "${s2_filter}" \
    "${raw_reads_id}" \
    "${s1_out_paf}" \
    "${target_genes_file}" \
    "${s2_out_pref}" \
    > "${output_reads_id}"

  seqtk subseq "${s2_inp_reads}" "${output_reads_id}" \
    | seqtk seq -A | gzip -c > "${output_reads_fa}"
}

step3() {
  echo "Step 3: Map filtered reads to per-gene allele sets"
  local s2_out_pref="${OUT_DIR}/s2_out"
  local s3_inp_reads="${s2_out_pref}.read-fa.gz"
  local list_dir="${s2_out_pref}_list"
  local s3_out_paf="${OUT_DIR}/s3_out.paf.gz"
  local mm2thread=20
  local template_gene_seq="${SCRIPT_DIR}/dps-dat/6gene.fa.gz"
  local target_genes_file="${GENE_LIST}"
  local mm2opt="-t ${mm2thread} -e10000 -k19 -w19 -g10k \
    -A1 -B4 -O6,26 -E2,1 -s350 -I10k -c --ds --end-bonus=10"

  mkdir -p "${list_dir}"
  > "${s3_out_paf}"

  # for each gene, build a little allele FASTA and read FASTA, then align
  while read -r GENE; do
    # gather allele IDs
    grep "^>${GENE}\\*" "${SCRIPT_DIR}/dps-dat/6gene.list" \
      | sed 's/^>//' \
      > "${list_dir}/${GENE}.alleles.ids"

    # if no alleles, skip
    if [[ ! -s "${list_dir}/${GENE}.alleles.ids" ]]; then
      echo "no alleles for ${GENE}, skipping" >&2
      continue
    fi

    # build allele reference FASTA
    seqtk subseq "${template_gene_seq}" "${list_dir}/${GENE}.alleles.ids" \
      | gzip -c > "${list_dir}/${GENE}.alleles.fa.gz"

    # build read FASTA for that gene
    seqtk subseq "${s3_inp_reads}" "${list_dir}/${GENE}.list" \
      | gzip -c > "${list_dir}/${GENE}.reads.fa.gz"

    # map reads → alleles, append to one big PAF
    minimap2 ${mm2opt} \
      "${list_dir}/${GENE}.reads.fa.gz" \
      "${list_dir}/${GENE}.alleles.fa.gz" \
      | gzip -c >> "${s3_out_paf}"
  done < "${target_genes_file}"

  >&2 echo "Step 3: all‐genes alignment finished"
}

step4() {
  echo "Step 4: Vote best alleles as templates"
  local IN_PAF="${OUT_DIR}/s3_out.paf.gz"
  local FILTER_SCRIPT="${SCRIPT_DIR}/src/filter_paf.py"
  local SEARCH_SCRIPT="${SCRIPT_DIR}/src/search_template.py"
  local PREF="${OUT_DIR}/s4_out"

  zcat "${IN_PAF}" | python3 "${FILTER_SCRIPT}" | gzip -c > "${PREF}_filtered.paf.gz"
  zcat "${PREF}_filtered.paf.gz" \
    | python3 "${SEARCH_SCRIPT}" "${PREF}"
}


step5() {
  s5_inp_mapping="${OUT_DIR}/s4_out_filtered.paf.gz"
  s5_inp_template_allele="${OUT_DIR}/s4_out_template_allelelist.txt"
  s5_inp_reads="${OUT_DIR}/s2_out.read-fa.gz"

  all_allele_seq="${SCRIPT_DIR}/dps-dat/6gene.fa.gz"
  call_variants="${SCRIPT_DIR}/src/call_variants.py"
  time python3 ${call_variants} \
  ${s5_inp_template_allele} \
  ${s5_inp_mapping} \
  ${s5_inp_reads} \
  ${all_allele_seq} \
  ${OUT_DIR}
  >&2 echo "s5 time:"

}


step6() {
  echo "Step 6: Realign reconstructed sequences"
  local QUERY="${SCRIPT_DIR}/dps-dat/6gene.fa.gz"
  local CONS="${OUT_DIR}/new_allele.fa"
  local OUT_PAF="${OUT_DIR}/s6_out.paf.gz"
  local THREADS=20
  local OPTS="-t ${THREADS} -s350 -c --cs --end-bonus=10 -I10"
  minimap2 ${OPTS} "${QUERY}" "${CONS}" | gzip -c > "${OUT_PAF}"
}

step7() {
  echo "Step 7: Final allele calls"
  local IN_PAF="${OUT_DIR}/s6_out.paf.gz"
  local REF_INFO="${SCRIPT_DIR}/dps-dat/gene_annot.info"
  local OUT_PAF="${OUT_DIR}/new_allele.paf"
  local FINAL_SCRIPT="${SCRIPT_DIR}/src/final_call.py"

  time zcat "${IN_PAF}" | python3 "${FINAL_SCRIPT}" "${REF_INFO}" > "${OUT_PAF}"
  awk '{ print $6, $0 }' "${OUT_DIR}/known_allele.paf" \
    > "${OUT_DIR}/known_allele.output.paf"
  cat "${OUT_PAF}" "${OUT_DIR}/known_allele.output.paf" \
    | sort -t$'\t' -k1,1
  >&2 echo "Final call complete"
}

main() {
  for i in {0..7}; do
    if [[ "${STEP[$i]:-false}" == true ]]; then
      step${i}
      log_time "$i"
    fi
  done
}

main "$@"
