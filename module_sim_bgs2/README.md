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
