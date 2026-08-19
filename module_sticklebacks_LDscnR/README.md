# module_sticklebacks_LDscnR

Clean, self-contained re-run of the 3-spine stickleback (3sp) LD-aware
outlier-region analysis using **only exported `LDscnR` functions** — kept
separate from the exploratory clutter in `module_sticklebacks/`.

## The point

3sp is the **saturated-genome** case: the whole genome carries real outlier
regions, so there is no neutral floor.

- **Single-SNP EMMAX is dead** here — nothing passes BH `q ≤ 0.05` (min `q ≈ 0.057`).
- The **LD-aware C-score** nonetheless resolves the known marine↔freshwater
  architecture (Chr1 inversion, Chr4 *Eda*, Chr7, Chr20, Chr17, …) as clean peaks.
- The **structure-aware null** certifies these regions: region-level FDR ≈ 0
  across the whole `tau_C × l_min` grid. There is *no* null-danger corner — the
  problem is not signal-vs-structure separation but **ranking many genuine
  regions**.
- So regions are ranked **threshold-free** by cross-parameter stability
  (`ld_region_stability()`), avoiding any arbitrary `maxC`-style cutoff.

## Run (two steps)

From the `LDscnR-paper` root (needs the `LDscnR` package installed). **First**
regenerate the input bundle, **then** run the analysis:

```bash
Rscript module_sticklebacks_LDscnR/regen_3sp_data.R    # -> module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds
Rscript module_sticklebacks_LDscnR/run_3sp_LDscnR.R    # -> figures/ + results/
```

## Data regeneration (`regen_3sp_data.R`)

Rebuilds the 3sp inputs from the **raw parsed genotypes**, using the **same
LDscnR chain as the simulations** (`module_sim/R/regen_stats.R`) so the empirical
and simulated analyses are methodologically identical. It fixes two ad-hoc
choices in the earlier exploratory 3sp code:

- **GRM pruning** — was "all SNPs with `ld_w_0.95 < 0.05`"; now the sim chain
  `ld_complexity_reduction(rho = 0.5)` → `ld_prune_and_eMLG` (rep-marker set).
- **LD-decay** — was the outdated `ld_decay()`; now `compute_LD_decay()` with the
  sim `DECAY_ARGS`.

The observed C-score and the structure-aware null both run through `emmax_fast()`
on the **one saved GRM** (the observed phenotype is scored alongside the surrogate
"permuted" phenotypes by the same fast engine — that shared kinship is the thing
that must be identical, and the bundle keeps it). LFMM is **unchanged** and reused.

Raw input: `~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData`
(`GTs_3sp` 117 × 881,786, `map_3sp` Chr/Pos/marker/maf, `pheno_3sp`) +
`.../3sp/lfmm_F.rds` (previous LFMM F, full map). MAF > 0.1 → 790,578 markers.

Output bundle `module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds`:
`GTs`, `map` (+`emx_*`, `lfmm_*`, `ld_w_095`), `eco`, `ld_ws`, `LD_decay`,
`complexity_reduction` (`stage1` + `pruned_markers`), `GRM`, `emx`, `emx_gif`,
`settings`.

## Outputs

`figures/`
- `fig1_manhattan_C_vs_q.png` — C-score vs −log₁₀(q) Manhattan, regions highlighted
  (the core demonstration: single-SNP flat, C-score resolved).
- `fig2_heatmaps_tau_lmin.png` — observed / structured-null / FDR region counts over
  `tau_C × l_min` (FDR ≈ 0 everywhere → saturated genome).
- `fig3_region_ranking.png` — regions recoloured by the two threshold-free ranks
  (`persist_tau` and grid-survival `stability`).

`results/`
- `region_stability_3sp.csv` — the ranked region table.
- `operating_point_3sp.csv` — the zero-null max-cluster region set.
- `fdr_grid_3sp.csv` — the full `tau_C × l_min` observed/null/FDR grid.

## Manhattan figures (`manhattan_regions.R`)

Genome-wide −log₁₀(q) Manhattans per method (EMMAX + LFMM) with the C-score
outlier regions coloured, in the original "sim-machinery" plotting style
(single-row chromosome facets, grey background, no x-axis, bold strips, the
`col_vec` palette). Reuses the migrated pipeline (bundle `ld_w`/GRM, structured
null for EMMAX, `ld_cscore` for LFMM, `ld_edges`/`ld_regions`) and the cached null
(`results/null_uncapped_3sp.rds`), so it renders in seconds. Parameterised on the
three clustering knobs:

```bash
Rscript module_sticklebacks_LDscnR/manhattan_regions.R 0.05 10 0.60   # reference look (r2~=0.4, l_min=10)
Rscript module_sticklebacks_LDscnR/manhattan_regions.R auto auto 0.75 # null-calibrated tau_C/l_min
```

`tau_C`/`l_min` accept `"auto"` (null-calibrated via `calibrate_tauc`/`calibrate_lmin`);
`rho_ld` is the decay-relative r² link (0.60≈r²=0.4, 0.75≈r²=0.27 on 3sp). Writes
`figures/manhattan_{EMMAX,LFMM,both}_tau*_lmin*_rho*.png`.

## Key `LDscnR` functions used

`emmax_setup` / `emmax_fast` → `structured_null` (C-score + surrogate bundle) →
`ld_edges` → `null_fdr` (grid) → `ld_regions` (operating point) →
`ld_region_stability` (threshold-free ranking) → `ld_manhattan` (figures).

`ld_region_stability()` and its data-prep helper `ld_cscore_scan()` were added to
the package for this analysis. The EMMAX path here gets its C-score and null from
`structured_null()`; `ld_cscore_scan()` is the no-null builder for methods without
a cheap null (e.g. LFMM — supply `pvals =` and transfer `tau_C` with
`gc_map_tauc()`).
