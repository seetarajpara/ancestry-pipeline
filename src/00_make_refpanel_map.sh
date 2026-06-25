#!/usr/bin/env bash
# =============================================================================
# 00_make_refpanel_map.sh
# Builds the FLARE "ref-panel" file: a 2-column whitespace table
#     <sampleID>  <panel>
# where <panel> is one of $PANELS (superpopulations). Reference samples NOT
# listed here are excluded from the panel, so this also subsets which 1000G
# populations FLARE uses as ancestry sources.
# =============================================================================
set -euo pipefail
source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../config.sh"

# Map 1000G 26 populations -> 5 superpopulations.
declare -A SUPER=(
  [CHB]=EAS [JPT]=EAS [CHS]=EAS [CDX]=EAS [KHV]=EAS
  [CEU]=EUR [TSI]=EUR [FIN]=EUR [GBR]=EUR [IBS]=EUR
  [YRI]=AFR [LWK]=AFR [GWD]=AFR [MSL]=AFR [ESN]=AFR [ASW]=AFR [ACB]=AFR
  [MXL]=AMR [PUR]=AMR [CLM]=AMR [PEL]=AMR
  [GIH]=SAS [PJL]=SAS [BEB]=SAS [STU]=SAS [ITU]=SAS
)

# Keep only requested panels.
declare -A WANT=(); for p in ${PANELS}; do WANT[$p]=1; done

# The ped/population file is whitespace-delimited with a header; columns include
# SampleID and Population. We detect those columns by name to be robust.
awk -v OFS='\t' '
  NR==1 {
    for (i=1;i<=NF;i++){ h[$i]=i }
    sid = (h["SampleID"]?h["SampleID"]:(h["sampleID"]?h["sampleID"]:h["Sample"]))
    pop = (h["Population"]?h["Population"]:h["population"])
    if (!sid || !pop){ print "ERROR: could not find SampleID/Population columns" > "/dev/stderr"; exit 1 }
    next
  }
  { print $sid, $pop }
' "${SAMPLE_META}" > "${REF_DIR}/1000G/_sample_pop.tsv"

# Emit sample -> superpop, filtered to requested panels.
: > "${REFPANEL_FILE}"
while read -r SAMP POP; do
  SP="${SUPER[$POP]:-}"
  [[ -z "${SP}" ]] && continue
  [[ -z "${WANT[$SP]:-}" ]] && continue
  printf '%s\t%s\n' "${SAMP}" "${SP}" >> "${REFPANEL_FILE}"
done < "${REF_DIR}/1000G/_sample_pop.tsv"

echo "Wrote ${REFPANEL_FILE}"
echo "Panel sample counts:"
cut -f2 "${REFPANEL_FILE}" | sort | uniq -c
echo
echo "NOTE: these sample IDs must exist in the reference VCFs. If FLARE later"
echo "complains about unknown samples, subset the ref VCF to these IDs or vice"
echo "versa (bcftools view -S ${REFPANEL_FILE%.*}.ids ...)."
cut -f1 "${REFPANEL_FILE}" > "${REFPANEL_FILE%.txt}.ids"
