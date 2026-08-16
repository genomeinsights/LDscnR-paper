# module_sim — simulated-data TP/FP validation of `ld_outlier_clusters()`

Ground-truth benchmark of the candidate-first LD-aware outlier caller
(`LDscnR::ld_outlier_clusters()`) on Nemo forward simulations, where the causal
QTNs are known. Companion to `module_sticklebacks/` (the empirical 3-species
analysis).

## Data

One `.rds` per replicate chromosome-pair: `list(GTs[160 × ~30k], map, env,
LD_decay, ld_ws)` with **precomputed** `emx_p/emx_F` (EMMAX, GRM on LD-pruned)
and `lfmm_p/lfmm_F` (LFMM K=5, genomic control), plus ground truth (`type`,
`Va`, `MAF`, `max_LD_with_QTN`, …).

- Full grid (external): `/Volumes/Nemo/Nemo_sim/parsed_sim_data2`,
  `adapt_bgs_chr{1..10}_V{0.5,1,2}_c{1,1.5,2}_env{1..10}.rds`; `ld_ws` columns
  are `rho_0.95`-prefixed.
- Local V2 subset: `parsed_sim_data`; `ld_ws` columns are bare `0.95`.

Set `LDSCNR_SIM_DATA` to point elsewhere; `_config.R` falls back to the local
subset automatically. `_config.R::pool_group()` handles both column layouts.

**`V` = inverse selection intensity** (V1 strongest, V2 collapses). **`c` =
dispersal** (c1 some gene flow → tractable; c2 very limited → structure-dominated).

## Pooling

`pool_group()` merges the 10 replicates sharing `(V, c, env)` into one
**20-chromosome genome** (relabel `Chr`/`marker` `R{i}_*`, `cbind` genotypes,
concatenate `ld_w`/decay). Same `env` ⇒ ~independent chromosomes, which
stabilises the RMSC threshold and gives **10 QTN-bearing + 10 neutral**
chromosomes as a TP/FP bed. A cluster is a **TP** if it is LD- and
distance-matched (`r2 > r2min`, `dist < dmax` from the decay fit) to a
Va-qualified `true_pos_QTN` — not by containment.

## Directory layout

```
module_sim/
  R/          live scripts (_config.R + 02..05, 10..14); run from LDscnR-paper/ as module_sim/R/<script>.R
  R/archive/  superseded single-q* RMSC chain (01, 06a, 06b, 07, 08, 09) — kept for reference
  data/       caches / scored results (*.rds) — regenerable, git-ignored
  figures/    plots (*.png)                    — regenerable, git-ignored
  doc/        this README
```
**Current direction:** the single-q\* RMSC caller/comparison (01, 06a–09) is superseded by the
**C-score** pipeline (10–14), which sweeps ρ/q\*/α and gates on consistency (tau_C) instead of
trusting one fragile threshold. The archived chain still holds the legacy-vs-ours comparison,
the ≥2-SNP result, and the RMSC-instability diagnostic that motivated the switch.

**Replicate-average by default.** A single env is unreliable — env1 has repeatedly produced
flukes (rosy candidate-first snapshot, a "low-ρ sweet spot", an "EMMAX≫LFMM α-effect") that all
vanished when averaged over 5 env. Run any C-score / performance analysis across env1–5 and report
mean ± SE; never make a claim from one env. (Confirmed robust: C-score beats single-SNP at every α
for both methods, ΔPR ≈ 0.05–0.16; *not* robust: α-stringency trend and EMMAX-vs-LFMM ordering.)
`R/_config.R` sets `dir_data`/`dir_fig` (created on demand); every script reads/writes there.

## Scripts (in `R/`, run from `LDscnR-paper/`, e.g. `Rscript module_sim/R/06a_run_caller.R 2 1 1`)

| | script | what |
|---|---|---|
| — | `_config.R` | paths, engine source, `pool_group()`, `score_thresholds()`, `cluster_regions()`, `evaluate_ORs_qtn()` |
| 01 | `01_score_pooled.R [V c env]` | candidate-first + background null, TP/FP + FP-location + diagnostics (default V1 c2 env3) |
| 02 | `02_ldw_vs_rawF.R [V c env]` | **mechanism figure**: raw F vs ld_w per chromosome (EMMAX flat/calibrated; LFMM inflated + LD-correlated on all chr incl. neutral) |
| 03 | `03_ldw_qtn_manhattan.R [V c env]` | ld_w Manhattan coloured by `max r²` with a QTN — why c2 fails (77% of high-ld_w is structure LD, not QTN linkage) |
| 04 | `04_structure_null.R` | structure-aware null on the stickleback demo → **over-corrects** (loses Eda); no single null spans both regimes |
| 05 | `05_distance_cap.R [V c env]` | decay-derived vs fixed-500kb cap → cap fragments clusters but hurts precision; structure diagnostics are cap-invariant |

### Legacy (single-SNP) vs ours, on a shared region frame

The real baseline is **plain single-SNP EMMAX / LFMM** (BH-FDR outliers, *no* null —
permuting a genome-wide single-SNP scan is impractical; the cluster null is feasible
for us only because candidate-first clustering collapses the test to a few clusters).
To compare fairly we build **one method-agnostic region frame**: union every method's
outlier SNPs, cluster once (`cluster_regions()`), label each region TP/FP by the same
QTN rule, and score each method by the regions it hits — so a single-SNP hit and an
`ld_w` hit at the same locus are the *same* detection. Two support axes matter:
**method support** (is a region shared across methods?) and **SNP support** (a
*single-SNP outlier* is a region a method tags with only one SNP — mostly false).

This is split into an expensive cache step and cheap scoring/plot steps:

| | script | what |
|---|---|---|
| 06a | `06a_run_caller.R [V c env]` | **expensive, run once**: pool → outlier sets (single-SNP + `ld_w`+null) → shared region frame → QTN-LD table → `cache_*.rds` |
| 06b | `06b_score.R [V c env] [MIN_SNP]` | **cheap, no genotypes**: score all/`>=MIN_SNP`-SNP/`>=2`-method views from the cache (MIN_SNP default 2) |
| 07 | `07_consensus_manhattan.R [V c env] [MIN_SNP]` | **cheap**: 4-panel Manhattan on the shared frame (colour = shared region identity; `+` = single-SNP outlier) + a `>=MIN_SNP` figure |
| 08 | `08_sweep_aggregate.R` | aggregate `consensus_*.rds` across env replicates → mean±SE Precision/Recall/PR/F1 per (condition, filter, method) + `figures/sweep_PR.png` |
| 09 | `09_rmsc.R` | RMSC (discoveries vs `ld_w` quantile) for every dataset × method → `figures/rmsc_all.png` (interior peak ⇒ `ld_w` filter adds power; needs only `ld_w`+p, no genotypes) |
| 10 | `10_cscore.R [V c env] [method]` | q\*-sweep **consistency (C-score)** as a q\*-robust alternative to a single RMSC q\*; C vs QTN-linkage diagnostic + C-gate PR sweep vs single-q\*/single-SNP |
| 11 | `11_rho_qstar_heatmap.R [V c env]` | **poster heatmap**: total outlier regions across the (ρ, q\*) LD-filtering grid (ρ = precomputed `ld_ws` column), EMMAX + LFMM |
| 12 | `12_rho_qstar_pr.R [V c env]` | **truth-aware** version: same grid coloured by PR (precision×recall, dedup-neutral) → the sweet spot. NB the best single (ρ,q\*) cell is UNSTABLE across env (ρ 0.05–0.95) → don't pick one, use the C-score |
| 13 | `13_cscore_2d.R [V c env] [method]` | 2-axis (ρ,q\*) **C-score** (poster's consistency def) → C-Manhattan + C-gated PR path; beats single-SNP, fixed ρ=0.95, and even the oracle best single cell |
| 14 | `14_alpha_cscore.R [V c env]` | fold **α** into the sweep; test that C-benefit grows as α tightens & is larger for EMMAX (conservative) than LFMM (inflated=lenient, saturates). l_min=2 as a post-C filter |
| 15 | `15_pr_auc.R [V c env]` | **headline**: standard trapezoidal **PR-AUC** — C-score (sweep tau_C, ρ/q\*/α folded in) vs single-SNP (sweep α); clustering fixed decay-relative ρ_ld=0.75/ρ_d=0.95≤500kb; l_min ∈ {1,2,4,8} as linetype. Retires the random-search AUC-PR\* (only needed for the raw unordered grid) |
| 15b | `15b_pr_auc_aggregate.R [V c]` | replicate-average 15 → mean±SE PR-AUC vs l_min + mean PR curves faceted by l_min (C-score coloured by tau_C) |
| 15c | `15c_alpha_need.R [V c]` | is the α sweep needed? PR-AUC of α-swept vs α=0.05-fixed C-score (mean 5 env, l_min=2) → EMMAX Δ≈0, LFMM marginal ⇒ **fix α=0.05 empirically, keep the sweep on sim as reviewer evidence** |
| 16 | `16_besttau_manhattan.R [V c env] [l_min]` | joint-OR-frame Manhattan at best tau_C vs best α — which ORs are shared; l_min=1 shows the single-SNP FPs C removes |
| 17 | `17_perm_null.R [V c env] [B]` | naive permutation null for tau_C → too permissive (breaks env↔structure confounding); FDR-tau_C≈0.06 ≪ PR-optimum |
| 18 | `18_structured_null.R [V c env] [B]` | **structured** orthogonal-spatial null (same autocorrelation, orthogonal to env) → prices structure confounding back in; empirically portable; fast EMMAX + per-marker-counter C-sweep (~8 s/surrogate) |
| 18b | `18b_structured_null_aggregate.R [V c]` | replicate-average 18 → **pooled** ratio-of-means FDR–tau_C over env1–5 (+ per-env spread) |
| 18d | `18d_null_bundle.R [V c env] [B]` | **null bundle**: one structured-null run saves per-surrogate sparse C + an edge cache over the union universe → FDR at ANY l_min later (l_min never enters null *generation*) |
| 18e | `18e_null_fdr.R [V c env]` | extract FDR(tau_C, l_min) from a bundle (18d) as a pure lookup; per-SNP + region curves; the calibrated tau_C falls as l_min rises |

**Ad-hoc cache + analysis layer** (touch genotypes once, then everything is a lookup — the 06a/06b split applied to the C-score chain):

| # | Script | Purpose |
|---|--------|---------|
| 20 | `20_build_cache.R [V c env…]` | **build the cache**: per-(V,c,env) `{map, LDW, decay_sum, th, qtab, edges}`; `cluster_from_cache` validated identical to `cluster_regions`. Downstream C-score/tau_C/PR-AUC need no genotypes |
| 21 | `21_estimate.R [V c]` | ad-hoc estimators from the cache (`prauc_cscore`, `load_cache`); reproduces the α-need PR-AUC in seconds. `source()` for the helpers |
| 22 | `22_tauc_by_method.R [V c]` | best tau_C per method (EMMAX ≈0.53 / LFMM ≈0.33 PR-optimum) + the C-distribution shift; a **shared tau_C≈0.35 costs <0.012** (flat plateau) |
| 23 | `23_calibrate_against_emmax.R [V c]` | **genomic control on the C-score**: quantile-map other methods onto EMMAX's null-anchored tau_C (one null calibrates all methods); LFMM-referenced tau_C is higher (controls its inflation) |
| 24 | `24_winning_range.R [V c]` | the tau_C range where the C-score PR beats single-SNP's best — ~half of [0.05,0.95] for both methods; the win is a **broad plateau, not a knife-edge** |
| 25 | `25_tauc_lmin_surface.R [V c]` | tau_C × l_min power surface: PR frontier + TP-vs-FP by l_min + the mechanism (high-C **singletons ~92% FP**, ≥10-SNP clusters ~73% TP). l_min=2 optimal on sim; l_min≥5 hurts (kills small true regions) |
| 26 | `26_lmin_truthfree.R [V c env]` | truth-free l_min diagnostic (single env): `est_TP = obs − null` vs real TP; the null estimates **yield well, precision only as a structure-FDR** |
| 26b | `26b_lmin_truthfree_aggregate.R [V c]` | **replicate-averaged** 26 (env1–5): the null-based objective picks the **same l_min optimum as the true PR-AUC** ⇒ truth-free l_min selection validated at the region level (pooled cor(est_TP,real_TP)=0.64) |
| 27 | `27_clustersize_vs_C.R [V c]` | cluster size vs C coloured TP/FP: TP-fraction rises monotone with size (**singletons ~99.9% FP**, 20+ ~60% TP) — the l_min mechanism visualised |
| 28 | `28_fold_validation.R [V c]` | **folded (poster) vs per-SNP (sim) C**: PR-AUC → per-SNP ≥ folded everywhere (folding neutral EMMAX / worse LFMM) ⇒ **keep per-SNP C + post-l_min, don't fold clustering/l_min into C** |
| 29 | `29_lfmm_gc_sim.R [V c]` | **validate the empirical calibration on the sim**: EMMAX structured-null tau_C → quantile-map (genomic control) to LFMM → score vs truth. LFMM finds **more true ORs** (5.2 vs 4.4, higher recall) at comparable precision (0.63 vs 0.61) ⇒ GC-map realises LFMM's power under EMMAX-anchored FP control (LFMM FP higher + more variable = power-vs-stability trade) |

**tau_C guidance & calibration.** Empirical PR-optimum tau_C ≈ 0.35–0.5 ("called in ≳half the analyses"),
stable across methods/env. **tau_C is not intrinsic** — it trades off with l_min and the α-grid (all reduce
FP), so the optimal tau_C *shifts* with l_min (weaker filter → higher optimal tau_C; PR-AUC itself is
l_min-invariant). This is why l_min stays a post-filter (out of C) and why any calibrated tau_C must be
reported *with* its (l_min, α-grid). Truth-free calibration = the structured null (18): naive permutation
(17) breaks the env↔structure confounding (→ tau_C≈0.06, too permissive); the orthogonal-spatial surrogate
keeps the spatial autocorrelation but removes the true signal, so its empirical FDR prices in structure FPs.
Structured-null result — **replicate-averaged over env1–5** (18b, pooled ratio-of-means, l_min=2, α∈{.001–.1}):
FDR≤0.05 at tau_C≈0.46, FDR≤0.10 at ≈0.24; the PR-optimum band 0.35–0.5 maps to pooled FDR 0.03–0.06, so the
performance-optimal and FDR-controlled tau_C **converge** there. **Aggregate with pooled ratio-of-means, not
mean-of-ratios** (per-env FDR = n_null/n_obs with n_obs 10–70 goes non-monotone; one low-discovery env inflates
the mean). The **per-env spread is large and mechanistic, not noise**: weak-signal replicates make no high-C
discoveries (nothing to calibrate); a high-structure replicate needs tau_C≈0.7 while clean ones need ≈0.02–0.18,
and the pooled value is set by the worst high-structure replicate that still makes discoveries. **Practical
consequence: do not transfer a universal tau_C constant — run the structured null on the target dataset itself**
(cheap now, ~13 min); the sim only certifies that the null yields a sensible, PR-optimum-aligned tau_C.
Uses **fast EMMAX** (`fast_emmax_setup`/`fast_emmax_p` in `_config.R`: eigendecompose K + rotate
genotypes once → whitened per-SNP F per phenotype; identical to `emmax()`, ~25× faster/phenotype). Surrogates
are Gaussian-kernel MVN over coords + Gram-Schmidt; **Moran-spectral (MSR) surrogates were tried and were
worse** here (0.69 vs 0.82 smoothness) because the env is near rank-1 (one dominant x-gradient MEM), so no
orthogonal surrogate can match its autocorrelation exactly — the Gaussian approach is near the achievable ceiling.

**Sim vs empirical validity notes.** Per-replicate GRMs correlate 0.31 (V2_c1, weak structure) vs 0.94
(V1_c2) vs 0.95 (3sp sticklebacks) — GRM correlation reads out structure *strength*, not segregation
independence. Per-simulation vs pooled-GRM EMMAX diverge (r≈0.65–0.75), so the pipeline correctly uses
per-simulation EMMAX and pools only at the C-score level; a pooled-GRM EMMAX would be wrong (10 independent
realizations sharing the same 160 grid-point individuals).

**PR-AUC, not AUC-PR\*.** The C-score collapses the unordered ρ×q\*×α grid into one *ordered* knob (tau_C),
so a standard monotonic PR-AUC (trapezoidal, `pr_auc()` in `_config.R`) applies — the random-search
cummax AUC-PR\* (`auc_cummax_PR` in the paper engine) was only needed because the raw method had many
unordered nuisance parameters. Clustering/TP-match thresholds are fixed decay-relative (ρ_ld=0.75 = the
match r², ρ_d=0.95 capped at 500 kb to avoid merging a whole chromosome in low gene flow); l_min is a
post-C filter, swept only to trace/compare curves.

TP/FP counting is **dedup-neutral** (`evaluate_ORs_qtn`): a region matching an
already-claimed true-positive QTN is dropped (neither TP nor FP), so performance is
robust to clustering-parameter-driven fragmentation. Outputs (`data/*.rds`,
`figures/*.png`) are regenerable and git-ignored.
Only 06a touches genotypes, so iterating on scoring or plots costs seconds.

## Headline findings

1. **Precision comes from the cluster-level null, not the ld_w filter** (marker
   precision 0.12→0.17 EMMAX / 0.23→0.28 LFMM with the filter alone; →0.67 both
   with the null).
2. **`ld_w` pairs with EMMAX, not LFMM.** EMMAX F is calibrated and flat in
   `ld_w`; LFMM F is baseline-inflated and rises with `ld_w` on **every**
   chromosome including the 10 neutral ones, so the filter concentrates LFMM's
   spurious inflation into high-LD regions → neutral false positives (fig. 02).
3. **`rho_d` is dataset-dependent** (sim decay ≈100× slower than the 3sp data;
   `rho_d=0.95` ≈ 700kb–2.2Mb here vs ≈50kb on 3sp).
4. **The dispersal-limited c2 regime is structure-dominated and no null rescues
   it** — correctly flagged by the two diagnostics baked into
   `ld_outlier_clusters()$diagnostics`: `ldw_tracks_structure` (F rises with
   `ld_w` on all chromosomes incl. neutral) and `whole_chr_clusters` (clusters
   span whole chromosomes). Tightening the distance cap does not help (05); the
   structure-aware null over-corrects the good regime (04). The honest output in
   c2 is "unreliable", not a forced call.
5. **The C-score is per-SNP; clustering + l_min stay post-filters.** Folding
   clustering+l_min *into* C (poster style) buys nothing on the sim — per-SNP C ≥
   folded C at every l_min (neutral for EMMAX, worse for LFMM) and costs more
   compute (28). Keep the per-SNP C, apply l_min afterward.
6. **l_min is a dataset-dependent post-filter with a clear mechanism.** High-C
   singletons are ~92–99% false; TP-fraction rises monotonically with cluster
   size (27), so l_min prunes FP preferentially — but too high an l_min kills small
   true regions (l_min=2 optimal on the sparse sim, l_min=10 on the dense 3sp; 25).
   The optimum is **selectable without truth**: the structured null estimates
   yield well enough that its l_min optimum matches the true PR-AUC's (26b), and
   SNP-density-per-LD-block sets the scale.
7. **α is fixed at 0.05 empirically** (the null prices in anticonservativeness), but
   **swept on the sim** as reviewer evidence that the sweep is redundant (15c). The
   sim keeps the α-sweep; empirical (`module_sticklebacks/`) fixes α=0.05.

## Empirical application — `module_sticklebacks/` (3sp)

The C-score machinery applied to the three-spine stickleback data (no ground
truth; deliverable = which loci, and comparison to the poster's own pipeline):
`09` observed per-SNP C; `10` structured null (genetic-structure surrogates
`y ~ MVN(0, GRM) ⊥ ecotype`); `11` per-SNP C Manhattan; `12` l_min=10 Manhattan;
`13` cluster size vs C (clustering **r²=0.1 / 0.5 Mb split** → Chr1 inversion = one
cluster, Chr4 = ~5 loci); `14` **reverse-engineers the poster Manhattan** with the
poster's own functions (folded cluster-membership C, `summarise_stability`); `15`
the same figure from the **sim per-SNP C** — LFMM ≈ identical to the poster,
EMMAX more conservative (the folded C borrows strength across cluster members).
The sim machinery transfers once its defaults are re-set for the dense data
(α=0.05, direct-r² clustering, l_min=10); the earlier "fails on 3sp" was
default-transfer, not the method.

**Full coherent run + engine reconciliation** (`16`–`18`). `16` runs the *complete*
sim machinery on 3sp with a single **reconciled engine** — fast EMMAX with one GRM for
both observed and null. An all-markers GCTA GRM over-corrects (absorbs the concentrated
ecotype signal); the coherent choice is a **GCTA GRM from low-ld_w markers (ld_w_095<0.05)**
— the recipe's neutral-background GRM — which recovers Chr4/Eda. Pipeline: per-SNP C (α=0.05)
→ genetic structured null → **l_min=10 region-FDR** (structure surrogates ≈never make ≥10-SNP
clusters, so the region filter is a strong clean FP control) → tau_C → ORs. EMMAX → **5 ORs
incl Eda/Chr4**. LFMM has no fast null, so it inherits EMMAX's calibration via **genomic
control** (quantile-map EMMAX's tau_C onto LFMM's C scale) → 29 ORs across 14 chr — its
higher power on the polygenic stickleback genome (validated on the sim, `29`), not inflation.
`17` regenerates the two-panel C-Manhattan from the saved run (no null re-run); `18` draws the
**poster-style −log10(q) Manhattan** with each OR a distinct colour (poster palette).
