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

## Scripts (run from `LDscnR-paper/`)

| | script | what |
|---|---|---|
| — | `_config.R` | paths, engine source, `pool_group()`, `score_thresholds()` |
| 01 | `01_score_pooled.R [V c env]` | candidate-first + background null, TP/FP + FP-location + diagnostics (default V1 c2 env3) |
| 02 | `02_ldw_vs_rawF.R [V c env]` | **mechanism figure**: raw F vs ld_w per chromosome (EMMAX flat/calibrated; LFMM inflated + LD-correlated on all chr incl. neutral) |
| 03 | `03_ldw_qtn_manhattan.R [V c env]` | ld_w Manhattan coloured by `max r²` with a QTN — why c2 fails (77% of high-ld_w is structure LD, not QTN linkage) |
| 04 | `04_structure_null.R` | structure-aware null on the stickleback demo → **over-corrects** (loses Eda); no single null spans both regimes |
| 05 | `05_distance_cap.R [V c env]` | decay-derived vs fixed-500kb cap → cap fragments clusters but hurts precision; structure diagnostics are cap-invariant |

Outputs (`*.rds`, `*.png`) are regenerable and git-ignored.

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
