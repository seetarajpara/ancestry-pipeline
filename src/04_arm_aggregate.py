#!/usr/bin/env python3
# =============================================================================
# 04_arm_aggregate.py
# Roll FLARE per-marker local ancestry up to chromosome-ARM admixture fractions.
#
# Method (piecewise-constant integration of the HMM calls):
#   - For each haplotype, between consecutive markers the ancestry is taken as
#     the call at the left marker. The genomic interval length is attributed to
#     that ancestry.
#   - Each marker is assigned to the p or q arm by comparing POS to the
#     centromere span (markers inside the centromere are skipped).
#   - arm fraction(ancestry A) = bp assigned to A over both haplotypes
#                                / (2 * total bp on that arm)
#
# Inputs come from config-driven paths but can be overridden on the CLI.
# Output:
#   <out_prefix>.arm_ancestry.long.csv   sample,chrom,arm,ancestry,fraction,bp,n_markers
#   <out_prefix>.global.csv              sample,ancestry,fraction  (genome-wide cross-check)
#
# Dependency-light on purpose: standard library + bcftools on PATH.
# =============================================================================
import argparse, glob, gzip, os, subprocess, sys
from collections import defaultdict

def read_centromeres(path):
    cen = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c, s, e = line.split()[:3]
            cen[c] = (int(s), int(e))
    return cen

def read_ancestry_codes(anc_vcf):
    """Parse the ##ANCESTRY=<NAME=int,...> meta line -> {int_code: name}."""
    op = gzip.open if anc_vcf.endswith(".gz") else open
    with op(anc_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##ANCESTRY="):
                body = line.strip().split("=", 1)[1].strip("<>")
                codes = {}
                for kv in body.split(","):
                    name, val = kv.split("=")
                    codes[int(val)] = name
                return codes
            if not line.startswith("##"):
                break
    raise SystemExit(f"ERROR: no ##ANCESTRY line in {anc_vcf}")

def samples(anc_vcf, bcftools):
    out = subprocess.run([bcftools, "query", "-l", anc_vcf],
                         capture_output=True, text=True, check=True)
    return out.stdout.split()

def stream_calls(anc_vcf, bcftools):
    """Yield (chrom, pos, [an1,an2 per sample...]) using bcftools query."""
    fmt = r"%CHROM\t%POS[\t%AN1\t%AN2]\n"
    p = subprocess.Popen([bcftools, "query", "-f", fmt, anc_vcf],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        chrom, pos = f[0], int(f[1])
        yield chrom, pos, f[2:]
    p.wait()
    if p.returncode != 0:
        raise SystemExit(f"ERROR: bcftools query failed on {anc_vcf} "
                         f"(are AN1/AN2 present? check FLARE output).")

def arm_of(pos, cen):
    s, e = cen
    if pos < s:  return "p"
    if pos > e:  return "q"
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--anc-glob", required=True,
                    help="glob for FLARE *.anc.vcf.gz files, e.g. results/TAG/*.flare.anc.vcf.gz")
    ap.add_argument("--centromeres", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--bcftools", default=os.environ.get("BCFTOOLS", "bcftools"))
    args = ap.parse_args()

    anc_files = sorted(glob.glob(args.anc_glob))
    if not anc_files:
        raise SystemExit(f"No files matched {args.anc_glob}")
    cen = read_centromeres(args.centromeres)
    codes = read_ancestry_codes(anc_files[0])
    samp = samples(anc_files[0], args.bcftools)
    ncode = len(codes)

    # bp[(sample, chrom, arm)][ancestry_name] = float bp ; mk = marker counts
    bp = defaultdict(lambda: defaultdict(float))
    mk = defaultdict(int)

    for af in anc_files:
        # gather markers per (chrom) so we can integrate intervals in order
        # store per sample/hap the (pos, ancestry) then integrate
        rows = list(stream_calls(af, args.bcftools))
        if not rows:
            print(f"  WARNING: no records in {af}", file=sys.stderr); continue
        # rows are already position-sorted within a chromosome VCF
        for i in range(len(rows) - 1):
            chrom, pos, calls = rows[i]
            chrom2, pos2, _ = rows[i + 1]
            if chrom2 != chrom:
                continue                      # don't span across chromosomes
            if chrom not in cen:
                continue
            arm = arm_of(pos, cen[chrom])
            arm2 = arm_of(pos2, cen[chrom])
            # only attribute an interval that lies wholly within one arm; this
            # drops the centromere-spanning gap between the last p and first q
            # marker, which otherwise inflates the p arm.
            if arm is None or arm != arm2:
                continue
            seglen = pos2 - pos
            if seglen <= 0:
                continue
            for si, s in enumerate(samp):
                an1 = calls[2 * si]
                an2 = calls[2 * si + 1]
                key = (s, chrom, arm)
                mk[key] += 1
                for an in (an1, an2):
                    if an in (".", ""):
                        continue
                    name = codes.get(int(an), f"anc{an}")
                    bp[key][name] += seglen

    # ---- write long per-arm table ------------------------------------------
    long_path = f"{args.out_prefix}.arm_ancestry.long.csv"
    glob_bp = defaultdict(lambda: defaultdict(float))   # (sample)->anc->bp
    with open(long_path, "w") as out:
        out.write("sample,chrom,arm,ancestry,fraction,bp,n_markers\n")
        for (s, chrom, arm), anc_bp in sorted(bp.items()):
            total = sum(anc_bp.values())
            if total == 0:
                continue
            for name in (codes[c] for c in sorted(codes)):
                frac = anc_bp.get(name, 0.0) / total
                out.write(f"{s},{chrom},{arm},{name},{frac:.6f},"
                          f"{anc_bp.get(name,0.0):.0f},{mk[(s,chrom,arm)]}\n")
                glob_bp[s][name] += anc_bp.get(name, 0.0)

    # ---- write genome-wide global cross-check ------------------------------
    glob_path = f"{args.out_prefix}.global.csv"
    with open(glob_path, "w") as out:
        out.write("sample,ancestry,fraction\n")
        for s in samp:
            tot = sum(glob_bp[s].values())
            for name in (codes[c] for c in sorted(codes)):
                frac = (glob_bp[s].get(name, 0.0) / tot) if tot else 0.0
                out.write(f"{s},{name},{frac:.6f}\n")

    print(f"Wrote {long_path}")
    print(f"Wrote {glob_path}  (compare against FLARE's own *.global.anc.gz)")
    print(f"Ancestries: {', '.join(codes[c] for c in sorted(codes))}")
    print(f"Samples: {len(samp)}  |  arm-records: {len(bp)}")

if __name__ == "__main__":
    main()
