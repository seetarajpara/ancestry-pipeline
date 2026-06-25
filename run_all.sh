#!/usr/bin/env bash
# =============================================================================
# run_all.sh  —  end-to-end driver (resumable, loud on failure)
# Loops over every query VCF in $QUERY_DIR and every chromosome in $CHROMS:
#   01 prep -> 02 phase -> 03 FLARE  (per sample-file, per chromosome)
# then 04 arm aggregation once per sample-file across all chromosomes.
#
# Design notes:
#   - NO global `set -e`: a failure in one chromosome must not silently kill the
#     whole run. Each step's exit code is checked explicitly and reported.
#   - RESUMABLE: a chromosome whose FLARE output already exists is skipped, so
#     re-running after a failure picks up where it left off. Delete
#     results/<tag>/<chrom>.flare.* to force a redo of that chromosome.
# =============================================================================
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${HERE}/config.sh"
mkdir -p "${OUT_DIR}"

shopt -s nullglob
QUERIES=( "${QUERY_DIR}"/${QUERY_GLOB} )
if [[ ${#QUERIES[@]} -eq 0 ]]; then
  echo "No query VCFs found in ${QUERY_DIR}/${QUERY_GLOB}"; exit 1
fi

run_step () {  # run_step "label" cmd args...   -> echoes status, returns rc
  local label="$1"; shift
  "$@"
  local rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "    !! FAILED: ${label} (exit ${rc})"
  fi
  return $rc
}

echo "Queries: ${#QUERIES[@]}   Chromosomes: ${CHROMS}"
overall_fail=0
for QVCF in "${QUERIES[@]}"; do
  TAG="$(basename "${QVCF}" | sed -E 's/\.vcf\.gz$//; s/\.vcf$//')"
  echo "######################################################"
  echo "# QUERY: ${TAG}"
  echo "######################################################"
  for CHR in ${CHROMS}; do
    FLARE_OUT="${OUT_DIR}/${TAG}/${CHR}.flare.anc.vcf.gz"
    if [[ -s "${FLARE_OUT}" ]]; then
      echo "----- ${TAG} : ${CHR}  (already done, skipping) -----"
      continue
    fi
    echo "----- ${TAG} : ${CHR} -----"
    run_step "01 prep ${CHR}"  bash "${HERE}/src/01_prep_query.sh" "${QVCF}" "${CHR}" || { overall_fail=1; continue; }
    run_step "02 phase ${CHR}" bash "${HERE}/src/02_phase.sh"      "${TAG}"  "${CHR}" || { overall_fail=1; continue; }
    run_step "03 FLARE ${CHR}" bash "${HERE}/src/03_flare.sh"      "${TAG}"  "${CHR}" || { overall_fail=1; continue; }
  done

  # arm aggregation across whatever FLARE outputs exist for this tag
  ANC=( "${OUT_DIR}/${TAG}"/*.flare.anc.vcf.gz )
  if [[ ${#ANC[@]} -eq 0 ]]; then
    echo "----- ${TAG} : no FLARE outputs, skipping arm aggregation -----"
    continue
  fi
  echo "----- ${TAG} : arm aggregation (${#ANC[@]} chromosome(s)) -----"
  run_step "04 arm aggregate" python3 "${HERE}/src/04_arm_aggregate.py" \
    --anc-glob "${OUT_DIR}/${TAG}/*.flare.anc.vcf.gz" \
    --centromeres "${CENTROMERES}" \
    --out-prefix "${OUT_DIR}/${TAG}/${TAG}" || overall_fail=1
done

echo
if [[ $overall_fail -eq 0 ]]; then
  echo "DONE (no failures). Per-sample outputs in ${OUT_DIR}/<tag>/:"
else
  echo "DONE, BUT SOME STEPS FAILED (see !! lines above). Re-running run_all.sh"
  echo "will skip completed chromosomes and retry the rest."
fi
echo "  <tag>.arm_ancestry.long.csv   <- LOCAL ancestry by chromosome arm"
echo "  <tag>.global.csv              <- GLOBAL ancestry (all chromosomes done so far)"
