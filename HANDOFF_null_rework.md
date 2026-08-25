# Handoff — null rework

Task brief for a fresh Claude Code session working in `~/gitlab/LDscnR-paper`.
Written 2026-08-25. Supersedes `notes.txt`, which is stale.

## Read these first, in this order

1. **`LDscnR_framework.html`** — the pipeline spec. §4 is the rework you are
   implementing; §8 is the task list this brief expands. Open it in a browser,
   it is a self-contained page.
2. **`manuscript/ld_aware_outlier_regions.tex`** — `\label{sec:null}` and
   `\label{sec:threenull}` are the sections that change.
3. `kingman2021/README.md` — the external truth set and its headline numbers.

## State

- Branch **`null-rework-two-basis`**, off `explore-c2-grid-stability` at `67f6984`.
- **8 commits are unpushed.** A previous session could not push (no GitHub
  credentials in its environment). Push before doing anything else.
- Working tree is clean apart from long-standing untracked scratch at the repo
  root (`Nemo_sim/`, `OR_performance/`, `Rplot*.png`, …). Leave it alone.

## The change, in one paragraph

The null previously did two jobs — calibrate `tau_C` via a pooled count-FDR, and
assign per-region significance via a location-matched test — across four
surrogate bases. The count-FDR job is removed: it only works when the null is
already silent, and on real data it returned `tau_C = NA, l_min = 1550` (spatial)
and `l_min = 157` (regional permutation) while the location-matched test on the
same regional null put all 17 EMMAX regions at `q_R <= 0.0149`. What remains is
**one instrument** (location-matched region p → BH → `q_R`), **one gate** (median
per-surrogate background counts, checked before any p-value is read), and **two
bases** (home-field specificity per engine; within-locality permutation for
attribution). `tau_C` is no longer calibrated — it is integrated away by the
second-tier C-score over the `(tau_C, l_min)` grid.

---

## Task 1 — Recompute the gate as median-per-surrogate

**Blocks everything else.** The two engines are currently reported on different
statistics, so the numbers in the manuscript are not comparable.

- EMMAX nulls: `module_sticklebacks_LDscnR/results/null_{uncapped,popperm,regionperm,spatial,latent}_3sp.rds`
  (`null_uncapped_3sp` is the genetic MVN basis). Each holds `$C_obs`,
  `$C_surr` (list of named C vectors) and `$universe`.
- LFMM nulls: `module_sticklebacks_LDscnR/results/lfmm_b1_fournulls_Cs.rds`.
- For each (engine × basis) compute, over surrogates, the **median** of
  `n(C > 0)` and of `n(regions >= l_min)` at `tau_C = 0.05, l_min = 3,
  rho_ld = 0.60, d_cap = 5e5`. Reuse the `ld_edges()` / `ld_regions()` pattern
  in `null_sig_landscape.R` (lines 32–42) — build the edge graph once.
- Write `module_sticklebacks_LDscnR/results/gate_background.csv` with observed
  counts alongside.

**Do not trust the `792` peaks-per-surrogate figure in the `.tex`** — it is a
median-per-surrogate number for EMMAX/spatial that has never been computed the
same way as LFMM's pooled `n(C>0)`. Reproduce it or replace it.

**Acceptance:** one table, one statistic, both engines, five bases.

## Task 2 — LFMM regional permutation to B = 200

`null_lfmm_region_perm_3sp.rds` currently has B = 96, so the p-floor is
`1/97 = 0.0103` and all 90 regions are tied at it. EMMAX runs at B = 200
(floor 0.005). They are not comparable.

- Runner: `lfmm_b1_fournulls_3sp.R` / `lfmm_permutation_null_3sp.R`.
- Cost: `lfmm_b1_fournulls.log` records 520.9 min for 5 whole-genome scans
  (~104 min each). 200 surrogates is a cluster job, not a laptop job — plan it
  as one.
- Budget comes from **not** running the LFMM spatial and global-permutation
  nulls, which the rework drops.

**Acceptance:** both engines share the 0.005 floor.

## Task 3 — Repoint the global permutation null to the regional one

Wider than it looks — **six files** reference `null_popperm_3sp.rds`:

| File | Line | What |
|---|---|---|
| `module_sticklebacks_LDscnR/null_sig_landscape.R` | 26 | hardcoded `NULLF`; output names at 20, 110, 126 also say `popperm` |
| `module_sticklebacks_LDscnR/region_empirical_pvals.R` | 35 | arg-4 default |
| `module_sticklebacks_LDscnR/manhattan_emp_pvals.R` | 26 | arg-4 default |
| `module_C2/R/00_helpers.R` | 24 | config default |
| `kingman2021/R/12_overlap_emmax17.R` | 19 | absolute `~/gitlab/...` path — fix the portability while you are there |
| `module_sticklebacks_LDscnR/null_fdr_landscape.R` | 20 | the four-null map; becomes two |

Leave `permutation_null_3sp.R:33` alone — it derives `TAG` from `STRATA` and
generates *both* nulls; that is the generator, not a consumer.

**This invalidates every `C⁽²⁾` currently in the repo.**
`results/region_c2_anchored.csv` and `results/region_stability2.csv` are against
the global null and must be regenerated. Do this before Task 4.

**Acceptance:** no consumer reads `null_popperm_3sp.rds`; every reported `C⁽²⁾`
is against the regional null.

## Task 4 — Does Kingman corroboration track `C⁽²⁾`?

The decisive test for how the LFMM claim gets written. Run it after Task 3.

- Rank the 90 LFMM regions (`region_emp_pvals_lfmm_region_perm_tau0.05_lmin3_rho0.60.csv`)
  by `C⁽²⁾` and ask whether EcoPeak corroboration is monotone in that rank.
- Truth sets: `kingman2021/data/peaks/` — `c155_global` (Global-specific) and
  `c150_pacNW` (Pacific-specific).
- **Coordinates differ.** 3sp is gasAcu1 (`Chr1..Chr21`, arabic); everything in
  `kingman2021/` is gasAcu1-4 (`chrI..chrXXI`, roman). Lift with
  `kingman2021/R/08_liftover.sh`; the chains are in `data/liftover/`.
- Follow the existing patterns: `kingman2021/R/18_cscore_enrichment.R` (per-SNP
  binned enrichment) and `R/17_overlap_sweep.R` (region-level, with the
  within-chromosome rotation null at B = 2000).
- Report Spearman of `C⁽²⁾` rank against the overlap indicator, plus the
  rotation-null fold enrichment per `C⁽²⁾` bin.

**What the answer decides.** ~21 loci are hit by LFMM at `l_min = 3` and by a
Kingman set while EMMAX misses them entirely (15 at n >= 3 sets) — that is real,
externally corroborated sensitivity. 78 loci are unique to the LFMM set and
unadjudicated. If corroboration is monotone in `C⁽²⁾`, that simultaneously
validates `C⁽²⁾` against an external truth set and resolves the tie the
permutation null cannot (every region sits at the p-floor). If it is flat, the
declining per-SNP enrichment (782 LFMM SNPs at C > 0.5 with 3.8× enrichment,
against EMMAX's 13 SNPs at 23×) has propagated to region level and the text must
say so.

## Tasks 5–7 — Text, once the numbers are in

5. **Reconcile the abstract.** It asserts LFMM "inflates on the panel's
   continental-scale structure". `lfmm_b1_fournulls.log` shows LFMM's nulls
   quieter than EMMAX's on *every* basis, spatial included (4 regions vs a
   median 792 peaks/surrogate). One of these is wrong; Task 1 decides which.
6. **Cite Bourgon, Gentleman & Huber (2010, PNAS)** in `sec:cscore` for
   independent filtering — filtering before multiple-testing correction
   preserves type-I error control when the filter statistic is null-independent
   of the test statistic. The neutral-SNP panel is the required diagnostic. It
   is currently presented as an empirical observation.
7. **Re-anchor the "pooled count is blunter than location-matched" argument** on
   the regional permutation (`l_min = 157`) rather than the global one.

---

## Task 0 — ATTRIBUTION: the C-score is Fang et al. (2021), not this paper

**Do this before writing any Methods prose.** `sec:cscore` currently reads as
though the consistency score is introduced here. It is not. Fang et al. (2021),
section "Assessing Sensitivity of Association Analyses to Parameter Settings",
on this same dataset:

> "With three parameters for defining LD clusters (|E|min, SNPmin, and Corth)
> and four methods to correct for multiplicity and P value inflation, the three-
> and nine-spined stickleback data sets were subjected to a total of 144 tests
> each. [...] we calculated a consistency score C for each putative outlier
> region, denoting the proportion of tests where a given genomic region was
> found significant, with C = 1 indicating that a given region was significant
> in all 144 tests. [...] We deemed outlier regions with C < 0.05 to be too
> sensitive to parameter settings to be considered further."

Already present in that paper, by name: the **C-score**, **tau_C = 0.05**,
**l_min = 10** ("at least ten unique loci"), **d_cap = 500 kb single-linkage**,
and the region-as-cluster-of-significant-loci definition. One of the four
corrections is the **permutation of Li et al. (2018)** -- the MVN(0, s2g*A +
s2e*I) surrogate -- feasible only because complexity reduction had shrunk the
test count.

Write `sec:cscore` as **following** Fang et al. (2021), and state what is
actually new here:

1. **What is integrated over.** Theirs: 3 clustering parameters x 4 correction
   methods, a menu of analysis choices. Ours: continuous nuisance parameters of
   the LD *filter* (rho, q*), natural [0,1] quantities, so the integration is
   assumption-free rather than over an arbitrary menu.
2. **What is tested.** Theirs: cluster-level units -- PCs in Li et al., SMLAs in
   Fang et al. ("a modified version of EMMAX that allowed us to test for
   associations between SMLAs rather than a single bi-allelic SNP at a time").
   Ours: SNP-level statistics preserved; LD filters *which* SNPs are tested.
3. **Two nested scales.** Theirs: one tier, at region level. Ours: marker C over
   (rho, q*), then region C2 over (tau_C, l_min) -- and tau_C and l_min are
   swept rather than fixed at 0.05 and 10.
4. **The null's role.** Theirs: permutation is one correction among four inside
   the score. Ours: the structure-aware surrogate is the inferential instrument
   (location-matched region p -> BH q_R), with bases chosen to separate
   specificity from attribution.
5. **Benchmarking.** The consistency principle has **never been tested against
   simulations with known causal variants**. Doing so is a genuine first here,
   and is probably the single cleanest novelty claim the paper has.

Getting this wrong is a reviewer-fatal risk: PK is an author on Fang et al.
(2021), so presenting the C-score as new reads as self-plagiarism rather than
oversight.

## Task 4b — Fang et al. 2021 is a method comparison, NOT a positive control

**The example data IS the Fang et al. (2021) panel.** Same individuals, same
four geographic strata (Baltic Sea, North Sea, Norwegian Sea, White and Barents
Seas -- exactly the strata in `regionperm_null.log`), and that paper ran **EMMAX
with relatedness as a random effect and FDR** on it, reporting 2,996 SNPs in
**26 outlier regions** for *G. aculeatus*.

So Fang 2021 cannot corroborate anything -- it is the same data and the same
association engine. What it *is*, is the tightest possible **method comparison**:
identical data, identical engine, different downstream treatment (LDna
clustering with parameter selection vs. C-score, integration, and a
structure-aware null). Run it that way:

- How many of Fang's 26 regions does LDscnR recover, at the reported operating
  point and across the grid?
- What does LDscnR find that the LDna-plus-selection approach missed, and is
  *that* set corroborated by Kingman?
- Fang et al. chose the parameter cell that "detected the most significant
  regions" (|E|min = 10, Corth = 0.5, SNPmin = [10, 20]). Their own reported
  numbers therefore sit at a selected optimum; C2 is the direct answer to that,
  and the comparison should say so rather than treating 26 as a fixed target.

**Consequence for the truth sets.** Kingman is now the *only* independent
control, and it is geographically mismatched: Fang et al. (2020) show
parallelism is heterogeneous, mostly Eastern-Pacific-specific, with trans-oceanic
sharing "restricted to a limited number of shared genomic regions, including
three chromosomal inversions", while Kingman's panels are global and
North-Eastern Pacific and our data is 47--75 deg N Atlantic. The two available
comparison sets therefore have complementary weaknesses -- Kingman is
independent but mismatched, Fang 2021 is matched but not independent -- and
neither alone is a clean positive control. Say that plainly in the text; it is
the honest position and it pre-empts the obvious reviewer objection to each.

It also weakens the case against LFMM's 78 unique regions: absence from a
Pacific-weighted set is weak evidence when low overlap is what Fang 2020
predicts.

## Task 8 — Name the simulator, and check the argument that leaned on "coalescent"

The `.tex` described the simulations as *coalescent* in four places (lines 48,
771, 802, 830) and never named the simulator. They are the **current NEMO** — forward-in-time, individual-based, Frederic
Guillaume's simulator; he is a co-author, so he can supply the exact version and
citation. Do **not** cite quantiNemo: that is a different program, Fred is no
longer part of its development, and it was the simulator used by Fang et al.
(2020), not here. `neuenschwander2008quantinemo` stays in `references.bib` only
for that Fang citation and is annotated as such.
The four occurrences are now "forward-in-time individual-based", but two things
remain:

1. **Methods must name NEMO**, with version and citation, and add it to
   `references.bib`. There is currently no mention of it anywhere in the
   manuscript.
2. **Re-check the argument at line 771.** It reads "On the [...] simulations,
   where population structure is simple and well captured by a handful of latent
   factors, LFMM is at least as powerful as EMMAX". That premise was implicitly
   doing work as *coalescent = simple structure*. With a forward-in-time
   individual-based model the structure is whatever was simulated, so the claim
   needs to rest on the actual simulated demography (dispersal regime, number of
   demes, migration) rather than on the class of simulator. Verify it against
   the NEMO configuration before it goes to review.

This is also a reason to be careful with the `l_min`/detectability text, which
contrasts a "gene-flow regime" against a "low-dispersal, structure-dominated
regime" — those are NEMO parameters and should be stated as such.

## Environment

- Run everything from the repo root.
- **Do not load `.RData`** (7.6 GB) or `.RDataTmp` (6.4 GB) — gitignored session junk.
- `module_sticklebacks_LDscnR/{data,figures,results}/` and the `module_sim_LDscnR`
  equivalents are **gitignored**. Outputs are not versioned; the scripts are the
  record. Anything you regenerate must be reproducible from a committed script.
- `kingman2021/data/*.rds` (540 MB) are gitignored; rebuild via
  `kingman2021/R/run_all.sh`, which needs the joint VCF outside the repo at
  `~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/`.
- `GTs` in the kingman bundles is **mean-imputed and numeric**, not integer
  0/1/2 — ~12% of calls are missing at `MAF>=0.05 / F_MISSING<=0.20` at 5.5×
  coverage. Expect r²/`ld_w` attenuation. `gcta_grm()` / `emmax_setup()` have no
  NA handling.
- Deps: `LDscnR`, `data.table`, `igraph`, `parallel`, `ggplot2`, `wesanderson`,
  `PRROC`, `patchwork`; parsing also `SNPRelate`, `LEA`, `factoextra`.

## Ground rules

- Check column names before using them. Prefer vectorised, R-idiomatic code —
  the codebase is `data.table` throughout, stay with it.
- Look syntax up or ask rather than guess.
- Double-check your work and show the check. Say "I don't know" when that is the
  honest answer.

## Do not

- **Re-open the design debate.** §4 of the framework is settled. If a number
  contradicts it, say so explicitly and stop — do not quietly revert to the
  four-null design.
- **Commit to `main`.** Work on `null-rework-two-basis`.
- **Trust numbers already in the `.tex`.** Several are known stale: the `792`
  peaks/surrogate (Task 1), the abstract's LFMM claim (Task 5), and an earlier
  `λ ≈ 1.0` for LFMM that was measured *after* genomic control and therefore
  could not diagnose inflation (the uncorrected value is ≈ 1.7).

## First action

`git push -u origin explore-c2-grid-stability && git push -u origin null-rework-two-basis`,
then Task 1.
