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
  R/        scripts (_config.R + 01..08); run from LDscnR-paper/ as module_sim/R/<script>.R
  data/     caches / scored results (*.rds) — regenerable, git-ignored
  figures/  plots (*.png)                    — regenerable, git-ignored
  doc/      this README
```
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
| 11 | `11_rho_qstar_heatmap.R [V c env]` | **poster heatmap**: outlier regions across the (ρ, q\*) LD-filtering grid (ρ = precomputed `ld_ws` column), EMMAX + LFMM |

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
