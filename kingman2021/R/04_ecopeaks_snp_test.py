#!/usr/bin/env python3
"""
Kingman et al. 2021 SNP-based EcoPeak test, re-implemented from the hub's
ecoPeakDescription.html:

  "we analyzed the distribution of allele counts between marine and freshwater
   populations at every base in the genome with two alleles present at >10%
   frequency in the combined analysis metapopulation. At each such base, after
   conditioning on the observed number of homozygous reference, heterozygous, and
   homozygous non-reference calls, we used a multivariate hypergeometric
   generalization of Fisher's Exact Test to compute a two-sided p-value for the
   probability of an imbalance in allele counts at least as extreme as observed
   occurring by chance"

So: hold the genotype-class totals (n0,n1,n2) and the group sizes (m marine,
f freshwater) fixed, ask how surprising the marine ALT-allele count is when the
N called genotypes are dealt at random into the two groups.

  P(a0,a1,a2) = C(n0,a0) C(n1,a1) C(n2,a2) / C(N,m),   a0+a1+a2 = m
  k = a1 + 2*a2                    (marine ALT allele count)
  p = P(|k - E[k]| >= |k_obs - E[k]|)

The null distribution of k depends only on (n0,n1,n2,m), so it is computed once
per distinct key and cached -- that is what makes this tractable genome-wide.

Usage:
  04_ecopeaks_snp_test.py <cohort> [min_maf] [min_callrate] [out_prefix]
Reads  gts/<cohort>.gt.tsv.gz + <cohort>.samples.txt, and the Table S2 ecotypes.
Writes <out_prefix>.snp_p.tsv.gz  (Chr Pos n0 n1 n2 m f k_obs p)
"""
import sys, gzip, csv, math
from math import comb
from collections import defaultdict

COH   = sys.argv[1] if len(sys.argv) > 1 else "c151_nEur"
MINMAF= float(sys.argv[2]) if len(sys.argv) > 2 else 0.10
MINCR = float(sys.argv[3]) if len(sys.argv) > 3 else 2/3
OUT   = sys.argv[4] if len(sys.argv) > 4 else None

DATA  = "/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021"
gtf   = f"{DATA}/gts/{COH}.gt.tsv.gz"
smf   = f"{DATA}/gts/{COH}.samples.txt"
OUT   = OUT or f"{DATA}/ecopeaks/{COH}"

samples = [l.strip() for l in open(smf) if l.strip()]
meta = {r["seq_id"]: r for r in csv.DictReader(open(f"{DATA}/meta/tableS2_samples.tsv"), delimiter="\t")}
is_marine = [meta[s]["ecotype"] == "M" for s in samples]
M_IDX = [i for i, v in enumerate(is_marine) if v]
F_IDX = [i for i, v in enumerate(is_marine) if not v]
NS = len(samples)
sys.stderr.write(f"{COH}: {NS} samples, {len(M_IDX)} marine, {len(F_IDX)} freshwater\n")

DOSE = {"0/0":0,"0|0":0,"0/1":1,"0|1":1,"1/0":1,"1|0":1,"1/1":2,"1|1":2}

_cache = {}
def pvals_for_key(n0, n1, n2, m):
    """two-sided p for every attainable k, given class totals and marine size m"""
    key = (n0, n1, n2, m)
    hit = _cache.get(key)
    if hit is not None:
        return hit
    N = n0 + n1 + n2
    denom = comb(N, m)
    Ek = m * (n1 + 2*n2) / N
    dist = defaultdict(float)
    for a1 in range(max(0, m-n0-n2), min(n1, m) + 1):
        c1 = comb(n1, a1)
        for a2 in range(max(0, m-n0-a1), min(n2, m-a1) + 1):
            a0 = m - a1 - a2
            if a0 < 0 or a0 > n0:
                continue
            dist[a1 + 2*a2] += c1 * comb(n2, a2) * comb(n0, a0)
    ks = sorted(dist)
    dev = {k: abs(k - Ek) for k in ks}
    out = {}
    for k in ks:
        d = dev[k]
        out[k] = min(1.0, sum(dist[j] for j in ks if dev[j] >= d - 1e-9) / denom)
    _cache[key] = out
    return out

import os
os.makedirs(os.path.dirname(OUT), exist_ok=True)
nin = nout = 0
with gzip.open(gtf, "rt") as fh, gzip.open(f"{OUT}.snp_p.tsv.gz", "wt") as out:
    out.write("Chr\tPos\tn0\tn1\tn2\tm\tf\tk_obs\tp\n")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        nin += 1
        gts = f[4:]
        n0 = n1 = n2 = 0
        km = 0; m = 0
        for i, g in enumerate(gts):
            d = DOSE.get(g)
            if d is None:
                continue
            if d == 0: n0 += 1
            elif d == 1: n1 += 1
            else: n2 += 1
            if is_marine[i]:
                m += 1; km += d
        N = n0 + n1 + n2
        if N < MINCR * NS or m == 0 or m == N:
            continue
        af = (n1 + 2*n2) / (2*N)
        if af < MINMAF or af > 1 - MINMAF:
            continue
        p = pvals_for_key(n0, n1, n2, m)[km]
        out.write(f"{f[0]}\t{f[1]}\t{n0}\t{n1}\t{n2}\t{m}\t{N-m}\t{km}\t{p:.6g}\n")
        nout += 1
        if nout % 250000 == 0:
            sys.stderr.write(f"  {nout} sites tested ({len(_cache)} cached keys)\n")
sys.stderr.write(f"{COH}: {nin} sites read, {nout} tested, {len(_cache)} distinct null keys\n")
