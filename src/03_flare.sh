#!/usr/bin/env bash
# =============================================================================
# 03_flare.sh  <sample_tag>  <chrom>
# Runs FLARE: phased reference + phased query + genetic map + panel map.
# Produces, per chromosome:
#   <chrom>.flare.anc.vcf.gz   local ancestry (AN1/AN2 per marker per sample)
#   <chrom>.flare.global.anc.gz  global ancestry proportions per sample
#   <chrom>.flare.model / .log
# =============================================================================
set -euo pipefail
source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../config.sh"

TAG="$1"; CHR="$2"
WORK="${OUT_DIR}/${TAG}"
GT="${WORK}/${CHR}.query.phased.vcf.gz"
REF_VCF="${REF_VCF_TMPL/\{CHR\}/$CHR}"
MAP="${MAP_TMPL/\{CHR\}/$CHR}"
OUTPFX="${WORK}/${CHR}.flare"

for f in "${GT}" "${REF_VCF}" "${MAP}" "${REFPANEL_FILE}"; do
  [[ -s "$f" ]] || { echo "MISSING input: $f"; exit 1; }
done

${JAVA} -Xmx${JAVA_MEM} -jar "${FLARE_JAR}" \
  ref="${REF_VCF}" \
  ref-panel="${REFPANEL_FILE}" \
  gt="${GT}" \
  map="${MAP}" \
  out="${OUTPFX}" \
  probs=true \
  nthreads=${NTHREADS}

${BCFTOOLS} index -f -t "${OUTPFX}.anc.vcf.gz" 2>/dev/null || true
echo "[$TAG / $CHR] FLARE done -> ${OUTPFX}.anc.vcf.gz (+ .global.anc.gz)"
echo "  ancestry code map (remember for interpretation):"
gzip -dc "${OUTPFX}.anc.vcf.gz" | grep -m1 '##ANCESTRY=' || true
