# VCF Ancestry Pipeline (local + global)

Takes one or more **VCF** files of study samples and produces:

- **Global ancestry** — admixture proportions per sample (genome-wide).
- **Local ancestry** — admixture proportions **per chromosome arm** (p/q).

It is built to be portable (mostly Java jars + `bcftools`), so it can be handed
to a collaborator and run on a workstation or an HPC node without compiling
heavy C/C++ tools.

## Pipeline at a glance

```
query.vcf.gz
   │  01_prep_query.sh    normalize, biallelic SNPs, harmonize chr naming
   ▼
   │  02_phase.sh         Beagle 5.5, phased AGAINST the 1000G panel
   ▼
   │  03_flare.sh         FLARE → local ancestry (.anc.vcf.gz)
   │                              + global ancestry (.global.anc.gz)
   ▼
   04_arm_aggregate.py    roll local calls up to p/q ARM fractions
```

Reference: **1000 Genomes** high-coverage phased panel (GRCh38).
Local + global ancestry: **FLARE** (Browning lab). Phasing: **Beagle 5.5**.

## Tools

| tool      | role                    | why it's here                          |
|-----------|-------------------------|----------------------------------------|
| bcftools  | VCF QC / normalize      | standard                               |
| Beagle 5.5| phasing                 | single jar, same lab as FLARE          |
| FLARE     | local + global ancestry | one tool gives BOTH outputs; fast      |
| python3   | arm aggregation         | stdlib only, no extra packages         |

`ADMIXTURE` is intentionally NOT required: FLARE's global file already gives
per-sample admixture proportions. If you want an independent global cross-check
later, ADMIXTURE in supervised mode is the usual companion.

## Setup

```bash
# 0. edit config.sh  (BUILD, CHR_PREFIX, CHROMS, PANELS)  — read it top to bottom
# 1. tools + references (large; for a first test it pulls only chr22)
bash download_references.sh
# 2. build the FLARE reference-panel map (sample -> superpopulation)
bash src/00_make_refpanel_map.sh
# 3. drop your query VCF(s) in ./test/  then:
bash run_all.sh
```

## Fast first test

`config.sh` ships with `CHROMS="chr22"`. chr22 is the smallest chromosome and
lets the whole chain run end-to-end in minutes. Grab ~5 samples as a single
multi-sample VCF (or 5 single-sample VCFs) restricted to chr22 and put them in
`./test/`. Once chr22 looks right, switch `CHROMS` to the full autosome list in
`config.sh` and re-run.

## Outputs (per query file, in `results/<tag>/`)

- `<tag>.arm_ancestry.long.csv` — tidy long format:
  `sample, chrom, arm, ancestry, fraction, bp, n_markers`
  This is the LOCAL ancestry result (admixture by chromosome arm), ready for
  ggplot.
- `<tag>.global.csv` — genome-wide proportions per sample (our cross-check).
- `<chrom>.flare.global.anc.gz` — FLARE's own global estimate (authoritative
  global; should closely match `<tag>.global.csv`).
- `<chrom>.flare.anc.vcf.gz` — raw per-marker local ancestry (AN1/AN2), if you
  want segment-level resolution rather than arm-level.

## Reference panel choices (read before a real run)

The default panel is the five 1000G superpopulations (AFR/EUR/AMR/EAS/SAS).
Important caveat for admixed-American work: **1000G "AMR" samples are themselves
admixed**, so they are not a clean Indigenous-American reference. For a proper
3-way African / European / Indigenous-American local ancestry decomposition you
want unadmixed Indigenous-American reference haplotypes (e.g. HGDP Maya, Pima,
Karitiana, Surui, Colombian, or a curated panel). Add them as a separate panel
label in `src/00_make_refpanel_map.sh` and include their VCF in the `ref=`.
FLARE2 (Oct 2025, same repo) is designed for exactly the case where a source
population is poorly represented by the reference — worth a look if the
Indigenous-American reference is thin.

## The things that actually break these pipelines

1. **Genome build.** The query, the reference panel, and the genetic map must
   all be the same build. Mixing GRCh37 and GRCh38 fails silently (0 overlapping
   sites or garbage calls). Set `BUILD` and confirm your friend's VCFs match.
2. **Chromosome naming.** `chr1` vs `1`. `CHR_PREFIX` in config controls
   harmonization; the map re-prefix step in `download_references.sh` keeps the
   genetic maps consistent.
3. **Phasing with small N.** Always phase against the reference panel (this
   pipeline does). Phasing 5 samples on their own would give poor haplotypes.
4. **Acrocentric p-arms.** chr13/14/15/21/22 have essentially no callable
   variants on the p arm; those arm rows will be empty/NA. That's expected.
5. **Sample ID collisions.** If a query sample ID also appears in 1000G, rename
   it first; otherwise phasing/FLARE can behave oddly.

## Provenance / not-yet-tested note

These scripts were written against the current FLARE (0.6.0), Beagle 5.5, and
the 1000G high-coverage GRCh38 release. They have been syntax-checked and the
arm-aggregation logic unit-tested on synthetic calls, but the full chain has not
been run against real downloads in this environment. Expect to shake out a path
or a file-naming detail on first run — the most likely spots are the 1000G FTP
filenames (step 5 of the downloader) and chromosome-naming harmonization.
