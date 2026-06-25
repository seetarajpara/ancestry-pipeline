#!/usr/bin/env bash
# =============================================================================
# 02_phase.sh  <sample_tag>  <chrom>
# Phases the prepped query against the 1000G reference panel with Beagle.
# Reference-based phasing matters a LOT here: with only ~5 query samples,
# statistical phasing on the query alone would be poor. Phasing against 1000G
# gives haplotypes consistent with the same panel FLARE uses downstream.
# Output: ${OUT_DIR}/<tag>/<chrom>.query.phased.vcf.gz
# =============================================================================
set -euo pipefail
source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../config.sh"

TAG="$1"; CHR="$2"
WORK="${OUT_DIR}/${TAG}"
IN="${WORK}/${CHR}.query.prepped.vcf.gz"
REF_VCF="${REF_VCF_TMPL/\{CHR\}/$CHR}"
MAP="${MAP_TMPL/\{CHR\}/$CHR}"
OUTPFX="${WORK}/${CHR}.query.phased"

for f in "${IN}" "${REF_VCF}" "${MAP}"; do
  [[ -s "$f" ]] || { echo "MISSING input: $f"; exit 1; }
done

${JAVA} -Xmx${JAVA_MEM} -jar "${BEAGLE_JAR}" \
  gt="${IN}" \
  ref="${REF_VCF}" \
  map="${MAP}" \
  out="${OUTPFX}" \
  impute=false \
  nthreads=${NTHREADS}

${BCFTOOLS} index -f -t "${OUTPFX}.vcf.gz"
echo "[$TAG / $CHR] phased -> ${OUTPFX}.vcf.gz"
