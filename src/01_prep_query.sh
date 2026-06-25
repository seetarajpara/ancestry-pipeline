#!/usr/bin/env bash
# =============================================================================
# 01_prep_query.sh  <query.vcf.gz>  <chrom>
# Cleans one query VCF for one chromosome so it is FLARE/Beagle-ready:
#   - subsets to the chromosome
#   - left-aligns/normalizes, splits multiallelics
#   - (optionally) restricts to biallelic SNPs
#   - harmonizes chromosome naming to $CHR_PREFIX
#   - drops sites with no ALT / all-missing genotypes
# Output: ${OUT_DIR}/<sample_tag>/<chrom>.query.prepped.vcf.gz
# =============================================================================
set -euo pipefail
source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../config.sh"

QVCF="$1"; CHR="$2"
TAG="$(basename "${QVCF}" | sed -E 's/\.vcf\.gz$//; s/\.vcf$//')"
WORK="${OUT_DIR}/${TAG}"; mkdir -p "${WORK}"
OUT="${WORK}/${CHR}.query.prepped.vcf.gz"

# Build a chromosome-rename map both directions so we can match the query to the
# reference's naming regardless of how the query is labelled.
RENAME="${WORK}/.chr_rename.txt"
if [[ "${CHR_PREFIX}" == "chr" ]]; then
  # map  1->chr1 ... and leave chrN untouched
  for n in $(seq 1 22) X Y MT M; do echo "$n chr$n"; done > "${RENAME}"
else
  for n in $(seq 1 22) X Y MT M; do echo "chr$n $n"; done > "${RENAME}"
fi

# Region string in the query's own naming is ambiguous, so we annotate-rename
# FIRST, then subset by the (now harmonized) chromosome.
TMP1="${WORK}/.${CHR}.tmp1.vcf.gz"

${BCFTOOLS} annotate --rename-chrs "${RENAME}" "${QVCF}" -Oz -o "${TMP1}"
${BCFTOOLS} index -f -t "${TMP1}"

FILTER_SNPS=""
if [[ "${KEEP_ONLY_BIALLELIC_SNPS}" == "1" ]]; then
  FILTER_SNPS="-m2 -M2 -v snps"
fi

# Normalize (needs the same reference FASTA the data was called against for a
# perfect left-align; without it we still split multiallelics safely).
${BCFTOOLS} view "${TMP1}" -r "${CHR}" ${FILTER_SNPS} -Ou \
  | ${BCFTOOLS} norm -m -any -Ou \
  | ${BCFTOOLS} view -e 'ALT="." || F_MISSING==1' -Oz -o "${OUT}"
${BCFTOOLS} index -f -t "${OUT}"
rm -f "${TMP1}" "${TMP1}.tbi"

N=$(${BCFTOOLS} index -n "${OUT}")
echo "[$TAG / $CHR] prepped sites: ${N}  ->  ${OUT}"
[[ "${N}" -eq 0 ]] && echo "  WARNING: 0 sites — check build match and CHR_PREFIX."
