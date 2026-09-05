# module_sim_bgs2 -- SUPERSEDED, frozen 2026-09-05

Everything in this directory was built reading `Nemo_out_bgs`/`Nemo_out_nobgs`
as the raw simulation source (`module_sim/R/00_config.R`'s `raw_nemo_bgs5`/
`raw_nemo_nobgs` at the time). That was the wrong data.

**What went wrong:** `Nemo_out_bgs`'s archives carry essentially no
deleterious mutational load -- checked directly, `adapt_bgs_chr1_V0.5_c1_env1`
has **3 deleterious loci in the whole genome**. That is far below even the
already-documented-as-unmeasurable `bgs2` parameterization (200 loci,
predicted Q1/Q5 contrast 1.3%, under the ~5% detection bar -- see memory
`nemo-bgs-unmeasurable-settings`), which is why this archived copy is named
after it. The correctly-parameterized dataset ("scenario B", 1000 deleterious
loci per chromosome) is `/Volumes/Nemo/Nemo_sim/bgs5/` -- 1720 deleterious
loci confirmed in the same spot-check -- covering 4 cells (V0.5_c1, V0.5_c2,
V1_c1.5, V2_c1), not the 7-9 `Nemo_out_bgs`/`Nemo_out_nobgs` covered.

**Consequence:** every BGS-effect result built on this data (genome-wide
contrast, Fst-vs-recombination, n_snp-vs-recombination) is a numerically
correct measurement of an arm that was never parameterized to show an
effect -- not a code bug, but not informative about real BGS either. The
detection-performance results (PR/Recall/F_beta, FP-vs-cluster-size) do not
depend on the bgs arm having a real deleterious-selection signal in the same
way and were not necessarily wrong, but were built against this same raw
source for consistency and should be re-verified against bgs5 if they matter
downstream.

**Live pipeline:** `module_sim/` now reads `/Volumes/Nemo/Nemo_sim/bgs5/`
directly (both tags live in that one directory, unlike the old separate
`Nemo_out_bgs`/`Nemo_out_nobgs`), restricted to the 4 cells bgs5 currently
covers, designed to pick up the remaining 3 cells automatically once PK adds
their archives to the same directory (grid drivers skip missing cells rather
than fail; nothing needs editing to extend the grid later).

Kept here as a permanent, non-regenerating record -- do not update the code
or figures in this directory again.

**Exception, added 2026-09-06 (PK):** the 7-cell nobgs grid here is exactly
the "not necessarily wrong" detection-performance data described above, and
is the ONLY full 7-cell grid available while bgs5 still covers just 4 -- PK
asked to keep improving the analysis against bgs2's nobgs arm in the meantime
so it "easily transfers to the other simulations later," using an external
audit of the manuscript's simulation section as the worklist. New,
audit-response-only scripts and outputs live in `R_audit/` and are additive:
they read the existing `results/pooled_pr.rds` (nobgs rows only) and do not
touch anything in `R/`, `R_figures/`, `results/`, or `figures/` that predates
this note. Any fix proven here first still needs to be (and, so far, has
been) ported into the live `module_sim/R*/` scripts by hand -- this directory
itself still does not regenerate.

**Known gap, 2026-09-06:** audit item 1's core fix (single-SNP arms score
each significant marker as its own size-1 region, not the enclosing
multi-marker Stage-1 unit -- PK: "Single SNPs go all the way to 1 not 2")
landed in the live `R/04_score.R` + `R/05_pool.R` AFTER this archive's
`cluster_detail` was generated. It cannot be recomputed from the existing
`pooled_pr.rds` (needs a real 04_score.R rerun against raw score objects that
no longer exist for bgs2) and is out of scope for `R_audit/`'s
recompute-only fixes. The archived and `R_audit/`-corrected single-SNP
FP-by-size numbers therefore both still lack any n_loci==1 rows and
understate single-SNP arms' true FP rate at the smallest size. Will be
correct by construction once/if bgs2 nobgs is ever rerun under the current
pipeline.
