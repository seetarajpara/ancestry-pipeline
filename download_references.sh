#!/usr/bin/env bash
# =============================================================================
# download_references.sh
# Pulls the tools + reference data the pipeline needs. Run once.
# Re-running is safe: it skips files that already exist.
#
# Downloads:
#   1. Beagle 5.5 jar          (phasing)              ~ small
#   2. FLARE jar               (local+global ancestry) ~ small  (built from src)
#   3. PLINK genetic maps      (GRCh38)               ~ small
#   4. 1000G phased panel      (per chromosome)       LARGE (~0.5-2 GB / chrom)
#   5. 1000G sample metadata   (population labels)    ~ small
#
# TIP: for your first test, comment out the autosome loop near the bottom and
# pull ONLY chr22 (already the default in $CHROMS). chr22 is the smallest and
# makes the full pipeline runnable end-to-end in minutes.
# =============================================================================
set -euo pipefail
source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/config.sh"

mkdir -p "${ANC_ROOT}/tools" "${REF_DIR}/maps" "${REF_DIR}/1000G/${BUILD}"

echo "============================================================"
echo " 1. Beagle 5.5 (phasing)"
echo "============================================================"
if [[ ! -s "${BEAGLE_JAR}" ]]; then
  # Beagle ships as a single dated jar; this is the 5.5 release.
  wget -c "https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar" \
       -O "${BEAGLE_JAR}"
else
  echo "  already present: ${BEAGLE_JAR}"
fi

echo "============================================================"
echo " 2. FLARE (local + global ancestry)  — built from source"
echo "============================================================"
if [[ ! -s "${FLARE_JAR}" ]]; then
  ( cd "${ANC_ROOT}/tools"
    rm -rf flare
    git clone --depth 1 https://github.com/browning-lab/flare.git
    javac -cp flare/src/ flare/src/admix/AdmixMain.java
    jar cfe flare.jar admix/AdmixMain -C flare/src/ .
  )
else
  echo "  already present: ${FLARE_JAR}"
fi

echo "============================================================"
echo " 3. PLINK genetic maps (${BUILD})"
echo "============================================================"
if [[ ! -s "${REF_DIR}/maps/plink.chr22.${BUILD}.map" ]]; then
  ( cd "${REF_DIR}/maps"
    wget -c "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.${BUILD}.map.zip"
    unzip -o "plink.${BUILD}.map.zip"
  )
  # The zip contains two variants:
  #   chr_in_chrom_field/    CHROM col = chr1..chrX, files named plink.chrchrN.*
  #   no_chr_in_chrom_field/ CHROM col = 1..X,       files named plink.chrN.*
  # Pick the one whose CHROM matches CHR_PREFIX, and normalize filenames to
  # plink.<chrom>.${BUILD}.map so MAP_TMPL resolves. No column rewriting needed.
  if [[ "${CHR_PREFIX}" == "chr" ]]; then
    SRC="${REF_DIR}/maps/chr_in_chrom_field"
  else
    SRC="${REF_DIR}/maps/no_chr_in_chrom_field"
  fi
  for f in "${SRC}"/plink.chr*."${BUILD}".map; do
    base=$(basename "$f")                 # plink.chrchr22.GRCh38.map  OR  plink.chr22.GRCh38.map
    rest=${base#plink.chr}                # chr22.GRCh38.map           OR  22.GRCh38.map
    newchr=${rest%%.*}                    # chr22                      OR  22
    cp "$f" "${REF_DIR}/maps/plink.${newchr}.${BUILD}.map"
  done
  echo "  normalized maps into ${REF_DIR}/maps/ (e.g. plink.chr22.${BUILD}.map)"
else
  echo "  already present: ${REF_DIR}/maps/"
fi

echo "============================================================"
echo " 4. 1000G sample metadata (population labels)"
echo "============================================================"
if [[ ! -s "${SAMPLE_META}" ]]; then
  # Pedigree/panel file with Sample -> Population -> Superpopulation.
  # Lives one level ABOVE the working/ phased dir.
  wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt" \
       -O "${SAMPLE_META}"
else
  echo "  already present: ${SAMPLE_META}"
fi

echo "============================================================"
echo " 5. 1000G phased reference panel (${BUILD}) — per chromosome"
echo "    Source: NYGC high-coverage 3202-sample phased SNV/INDEL release"
echo "============================================================"
BASE="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
for CHR in ${CHROMS}; do
  OUT="${REF_DIR}/1000G/${BUILD}/${CHR}.1000G.phased.vcf.gz"
  if [[ -s "${OUT}" ]]; then echo "  already present: ${CHR}"; continue; fi
  # File naming on the FTP uses the form below; verify once in a browser if a
  # chromosome 404s (release file names occasionally get revised).
  FN="1kGP_high_coverage_Illumina.${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  echo "  downloading ${CHR} ..."
  wget -c "${BASE}/${FN}"     -O "${OUT}"
  wget -c "${BASE}/${FN}.tbi" -O "${OUT}.tbi"
done

echo
echo "All requested downloads complete."
echo "Next: bash src/00_make_refpanel_map.sh"
