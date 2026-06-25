#!/usr/bin/env bash
# =============================================================================
# config.sh  —  central configuration for the VCF ancestry pipeline
# Edit THIS file only. Every other script reads its paths/params from here.
# =============================================================================

# ---- Project root (auto-detected; no need to edit) --------------------------
export ANC_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# ---- Genome build -----------------------------------------------------------
# MUST match the build of your friend's query VCFs AND the reference panel.
# This is the #1 cause of silent failures. GRCh38 is the default here.
export BUILD="GRCh38"

# ---- Chromosome naming in the QUERY VCFs ------------------------------------
# "chr"  -> CHROM looks like chr1, chr2 ... (1000G NYGC GRCh38 panel uses this)
# ""     -> CHROM looks like 1, 2 ...
# Ref panel, genetic map, and query MUST all agree. 01_prep_query.sh harmonizes
# the query to match CHR_PREFIX; set this to whatever the REFERENCE uses.
export CHR_PREFIX="chr"

# ---- Which chromosomes to run -----------------------------------------------
# For a fast first test, set this to just one small chromosome, e.g. "chr22".
# For a real run: "chr1 chr2 ... chr22" (autosomes; arm-level ancestry is
# autosome-only by design).
# export CHROMS="chr22"
# Full autosome set (uncomment when ready):
export CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"

# ---- Reference panel populations (FLARE ancestries) -------------------------
# 1000G superpopulations. Drop any you don't want to model.
# AFR=African, AMR=admixed American, EAS=East Asian, EUR=European, SAS=South Asian.
# NOTE: 1000G "AMR" samples are themselves admixed — they are NOT a clean
# Indigenous-American reference. For true 3-way (AFR/EUR/Indigenous American)
# local ancestry you'll want to add HGDP Native-American samples as a 5th panel.
# See README "Reference panel choices".
export PANELS="AFR EUR AMR EAS SAS"

# ---- Tool locations (filled in by download_references.sh) -------------------
export BEAGLE_JAR="${ANC_ROOT}/tools/beagle.jar"
export FLARE_JAR="${ANC_ROOT}/tools/flare.jar"
export BCFTOOLS="bcftools"     # assumes bcftools on PATH; else absolute path
export JAVA="java"             # needs Java 11+ for FLARE; 8+ for Beagle

# ---- Reference data (filled in by download_references.sh) -------------------
export REF_DIR="${ANC_ROOT}/refs"
# Per-chromosome phased reference VCF, with {CHR} substituted at runtime:
export REF_VCF_TMPL="${REF_DIR}/1000G/${BUILD}/{CHR}.1000G.phased.vcf.gz"
# Per-chromosome PLINK genetic map, {CHR} substituted at runtime:
export MAP_TMPL="${REF_DIR}/maps/plink.{CHR}.${BUILD}.map"
# FLARE ref-panel file (sample -> superpop), built by 00_make_refpanel_map.sh:
export REFPANEL_FILE="${REF_DIR}/1000G/refpanel.${BUILD}.txt"
# 1000G sample metadata (sample -> population -> superpopulation):
export SAMPLE_META="${REF_DIR}/1000G/1000G_samples.tsv"
# Centromere BED for arm assignment (shipped in repo):
export CENTROMERES="${ANC_ROOT}/resources/centromeres.${BUILD}.bed"

# ---- Query (your friend's data) ---------------------------------------------
# Directory holding the query VCF(s). For the 5-sample test, drop them here.
# Can be one multi-sample VCF or several single-sample VCFs (the driver globs).
export QUERY_DIR="${ANC_ROOT}/test"
export QUERY_GLOB="*.vcf.gz"

# ---- Output -----------------------------------------------------------------
export OUT_DIR="${ANC_ROOT}/results"

# ---- Compute ----------------------------------------------------------------
export NTHREADS=4
export JAVA_MEM="8g"          # bump for full autosomes / large panels

# ---- QC thresholds for query prep ------------------------------------------
export MIN_QUERY_MAF=0        # 0 keeps all sites; FLARE handles rare variants
export KEEP_ONLY_BIALLELIC_SNPS=1
