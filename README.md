# LDscnR-paper — analysis code for the LDscnR manuscript

Code for every analysis, figure, table and reported value in the LDscnR manuscript (LD-aware
outlier regions for genotype–environment association scans):

| folder | what |
|---|---|
| `module_3sp/` | three-spined stickleback (*Gasterosteus aculeatus*) empirical analysis: Stage-1 LD units, EMMAX (consensus and Simes), LFMM, structure-aware permutation null, EcoPeak relocation tests, sensitivity, floor-stability labels, LD–recombination validation |
| `module_9sp/` | nine-spined stickleback (*Pungitius pungitius*) empirical analysis: same workflow, ecotype adjusted for lineage, within-lineage permutation, sensitivity and variance-component diagnostics |
| `module_sim_3sp53/final_analysis/` | the simulation study (1,400 Nemo runs: 7 selection/migration cells × 2 BGS treatments × 10 maps × 10 environments): parsing, Stage 1, EMMAX/LFMM, truth scoring, pooled precision/recall, floor and null-calibration analyses |
| `module_sim_3sp53/` (rest) | earlier simulation diagnostics on the same runs (structured nulls, BGS windows); not used for the reported results except where a script header says so |
| `3sp_data/rec_maps/` | stickleback recombination map from three F2 crosses (Glazer et al. 2015; Roesti et al. 2013), with its build scripts |
| `kingman2021/` | EcoPeak intervals (Roberts Kingman et al. 2021) lifted to gasAcu1 (`data/liftover/`) and the liftover script (`R/08_liftover.sh`) |
| `manuscript_scripts/` | in the archived copy only: the manuscript repository's value/table generators and the four figure scripts that live there |
| `R/00_fetch_inputs.R` | downloads and verifies the empirical input files |

Each module is a numbered pipeline. Every stage writes `out/<stage>/` plus a `_receipt.rds`
recording its inputs (SHA-256), parameters and the LDscnR version, so a stage reruns only when
something it depends on changed (`FORCE=1` reruns anyway). Every parameter is set in the module's
`R/00_config.R`.

## Software

- R 4.6.1 (aarch64-apple-darwin23); Bioconductor 3.23.
- **LDscnR v0.9.1** (tag `v0.9.1`, commit `2477619`; MIT). The scripts call
  `devtools::load_all()` on its source, because a few use internal functions. Get the source
  from the archived copy or `git clone --branch v0.9.1 https://github.com/genomeinsights/LDscnR`,
  install it as well (`remotes::install_github("genomeinsights/LDscnR@v0.9.1")`), and point
  `LDSCNR_SRC` at the source folder (default `../LDscnR`, beside this repository). Each module
  checks the source against a SHA-256 of `R/` (`check_ldscnr()` in `00_config.R`) and stops on a
  mismatch.
- CRAN/Bioconductor packages, versions used: data.table 1.18.2.1, ggplot2 4.0.3, patchwork 1.3.2,
  SNPRelate 1.46.0, gdsfmt 1.48.0, LEA 3.24.0 (LFMM), digest 0.6.39, ggrastr 1.0.2, pROC 1.19.0.1,
  ggnewscale 0.5.2, devtools 2.5.2, jsonlite (for `R/00_fetch_inputs.R`).
- Simulations: Nemo 2.4.3 with the macOS large-population patch, run from the
  [LDscnR-NEMO](https://github.com/genomeinsights/LDscnR-NEMO) parameter set (see below).
- Optional, only to rebuild the EcoPeak liftover: UCSC `liftOver`, `bigBedToBed`, `chainSwap`;
  Python 3 for `3sp_data/rec_maps/build_glazer_map.py`.

## Paths

Nothing is hard-coded to one machine. Run scripts from the repository root, or set:

| variable | default | meaning |
|---|---|---|
| `LDSCNR_PAPER_ROOT` | `.` | this repository (the unpacked archive) |
| `LDSCNR_SRC` | `<root>/../LDscnR` | LDscnR v0.9.1 source |
| `LDSCNR_NEMO_PARAMS` | `<root>/../LDscnR-NEMO/params_3spC_53cM` | Nemo maps, environments, `rds/rec_map*.rds` |
| `LDSCNR_NEMO_RAW` | `<root>/nemo_raw/prod_out` | raw Nemo production output (not archived) |
| `LDSCNR_SIM_PARSED` | `<root>/nemo_raw/module_sim_3sp53_parsed_final_v1` | parsed simulation bundles |
| `LDSCNR_LDPP` | `<root>/../LD-pruning-paper` | inputs of Table 6 only (see below) |
| `LDSCNR_MANUSCRIPT` | `.` | where `manuscript_scripts/R_figures/*.R` write figures |

The shell drivers (`module_sim_3sp53/**/run_*.sh`) set `LDSCNR_PAPER_ROOT` themselves.

## Inputs

**Empirical genotypes.** `module_3sp/data/3sp_data.RData`, `module_3sp/data/lfmm_F.rds`,
`module_9sp/data/9sp_data.rds` and `module_9sp/data/lfmm_F.rds` derive from Fang et al. (2021,
*Mol Biol Evol* 38:4205, doi:10.1093/molbev/msab144; raw data doi:10.5281/zenodo.4722879) and are
archived byte for byte in doi:10.5281/zenodo.18466057 (`3sp.zip`, `9sp.zip`). To fetch, verify
against `PROVENANCE.csv` and place them:

```bash
Rscript R/00_fetch_inputs.R --install
```

**EcoPeaks** (`kingman2021/data/liftover/pv_c155.specific.bed`, `pv_c150.specific.bed`) and the
**recombination map** (`3sp_data/rec_maps/stickleback_recomb_3crosses_gasAcu1.tsv`) are included.

**Simulations.** The raw Nemo output (1,540 run directories, ~125 GB) is not archived. It is
regenerated from LDscnR-NEMO: `params_3spC_53cM/` holds the maps, QTN/deleterious positions,
allelic values, environment surfaces and dispersal matrices, and `make_prod.sh` records the exact
`R/make_run.R` invocation that writes the 1,540 `.ini` files (Nemo seeds are embedded in them:
`110101 + 1000·chromosome + 100·cell + environment`). Run them with Nemo 2.4.3, set
`LDSCNR_NEMO_RAW` to the output and `LDSCNR_NEMO_PARAMS` to the parameter folder, then run
`module_sim_3sp53/final_analysis/run_full_grid.sh`.

**Large intermediate outputs** (the Stage-1 bundles, scan objects and per-run simulation outputs)
are deposited on Dryad in the same folder layout as this repository; unpack that archive into the
repository root to reproduce figures and values without rerunning the scans. The small summaries
(`module_*/out/**/*.tsv|csv`, `module_sim_3sp53/final_analysis/results/`) are already here.

## Running

Three-spined (`module_3sp/R/`, from the repository root, in this order):

```
01_inputs  02_bundle  03_EMMAX  04_lfmm  09_ld_recombination  10_sensitivity
11_table_sensitivity  12_table_regions  13_goldilocks_units  13b_goldilocks_regions
14_floor_tiers  15_marker_postfilter  16_ldw_peaks  17_span_vs_simulation  05_manuscript
R_figures/figure_manhattan  R_figures/fig_ld_recombination  R_figures/fig_ldw_vs_q
```

e.g. `Rscript module_3sp/R/03_EMMAX.R`. `02_bundle` (LD decay, Stage 1, GRM) takes hours; the
scans take minutes. Nine-spined (`module_9sp/R/`): `01_inputs  02_bundle
02b_env_structure_alignment  03_EMMAX  03b_EMMAX_sensitivity_ELAdmixed  03c_EMMAX_structurednull
03d_diagnose_variance_components  04_floor_tiers`, then `R_figures/*.R`.

Simulations (`module_sim_3sp53/final_analysis/`): `run_full_grid.sh` (per run: `R/run_combo.R` =
`01_parse_nemo`, `02_build_ld_units`, `03_emmax`, `05_score_truth`), `run_lfmm_grid.sh`
(`04_lfmm`), then `R/06_summarise.R` … `R/30_truth_sensitivity_by_group.R` in numeric order
(`run_rerun_sequence.sh` chains the core ones), then `R_figures/*.R`. `final_analysis/README.md`
describes every stage.

Manuscript values and tables: `Rscript manuscript_scripts/generate_empirical_values.R <manuscript dir>`
and `Rscript manuscript_scripts/generate_simulation_manuscript.R <manuscript dir>` write
`values_empirical.tex`, `values_simulation.tex` and the `tables/tableS_*.tex` files;
`module_3sp/R/05_manuscript.R` writes `values_panel_module3sp.tex`.

**Seeds.** LD-decay fits `seed = 1` (`SEEDS[["bundle"]]`); permutation surrogate *b* uses
`set.seed(b)`, b = 1…B (B = 1,000 consensus, 200 Simes); EcoPeak relocations `seed = 1`
(10,000 draws); simulation bootstraps `SEEDS[["bootstrap"]] = 41`; the remaining fixed seeds are
set at the top of the script that uses them (`grep set.seed`).

## Figures and tables

| manuscript file | script (writes) |
|---|---|
| **Main text** | |
| `figures/figure1_ld_units.pdf` | `manuscript_scripts/R_figures/figure1_ld_units.R` |
| `figures/figure_simulation_performance.pdf` | `module_sim_3sp53/final_analysis/R_figures/figure_simulation_performance.R` |
| `figures/figure4_stickleback_validation.pdf` | `module_3sp/R_figures/figure_manhattan.R` → `module_3sp/R/05_manuscript.R` |
| `figures/figure6_permutation_structure.pdf` | `manuscript_scripts/R_figures/figure6_permutation_structure.R` |
| **Supplementary figures** | |
| `figureS_decay_vs_size.pdf` | `manuscript_scripts/R_figures/figureS_decay_vs_size.R` |
| `figureS5_chr1_stage1_stage2.png`, `figureS6_chr4_stage1_stage2.png` | `module_3sp/R_figures/fig_ld_recombination.R` (`Chr1/Chr4_stage1_vs_combined_high.png`) |
| `figureS_simulation_analysis_flow.pdf` | `manuscript_scripts/R_figures/figureS_simulation_analysis_flow.R` |
| `figureS_simulation_environment_sampling.pdf` | `manuscript_scripts/R_figures/figureS_simulation_environment_sampling.R` |
| `figureS_qtn_detection.pdf` | `final_analysis/R_figures/figureS_qtn_detection.R` |
| `figureS_simulation_manhattan_unrestricted.pdf`, `…_stage2.pdf` | `final_analysis/R_figures/figureS_simulation_manhattan_tpfp_unrestricted.R`, `…_tpfp_stage2.R` |
| `figureS_simulation_absolute_performance.pdf` | `final_analysis/R_figures/figureS_simulation_absolute_performance.R` |
| `figureS_simulation_bgs_qc.pdf` | `final_analysis/R_figures/figureS_simulation_bgs_qc.R` |
| `figureS_simulation_truth_sensitivity.pdf` | `final_analysis/R_figures/figureS_simulation_truth_sensitivity.R` |
| `figureS_simulation_floor_decomposition_c1.pdf` | `final_analysis/R_figures/figureS_simulation_floor_decomposition_c1.R` |
| `figureS_simulation_floor_decomposition_by_cell.pdf` | `final_analysis/R_figures/figure_simulation_floor_decomposition_by_cell.R` |
| `figureS_simulation_stage2_size_truth.pdf` | `final_analysis/R_figures/figureS_simulation_stage2_size_truth.R` |
| `figureS_simulation_null_truth_calibration.pdf` | `final_analysis/R_figures/figure_simulation_null_truth_calibration.R` |
| `figureS1_decay_recombination.pdf` | `module_3sp/R_figures/fig_ld_recombination.R` (`figure5_ld_recombination_reproduced.pdf`) |
| `figureS2_roc_low_recombination.pdf` | `module_3sp/R_figures/fig_ld_recombination.R` (`figS_roc_low_recombination.pdf`) |
| `figureS_nwin100_chr1.png` | `module_3sp/R_figures/fig_ld_recombination.R` (`figS_lddecay_nwin100_illustrative.png`) |
| `figureS_empirical_env_structure_alignment.pdf` | `module_9sp/R_figures/figure_env_structure_alignment.R` |
| **Tables** | |
| `tables/table4_reported_regions.tex` | `module_3sp/R/12_table_regions.R` |
| `tables/table5_sensitivity.tex` | `module_3sp/R/10_sensitivity.R` → `11_table_sensitivity.R` |
| `tables/table6_goldilocks_units.tex` | `module_3sp/R/13_goldilocks_units.R` (inputs: see below) |
| `tables/tableS_tier_ecopeak.tex`, `tableS_marker_postfilter.tex` | `manuscript_scripts/generate_empirical_values.R` (from `14_floor_tiers`, `15_marker_postfilter`) |
| `tables/tableS_simulation_performance.tex`, `tableS_chromosome_split.tex` | `manuscript_scripts/generate_simulation_manuscript.R` |
| `values_panel_module3sp.tex`, `values_empirical.tex`, `values_simulation.tex` | `module_3sp/R/05_manuscript.R`, `generate_empirical_values.R`, `generate_simulation_manuscript.R` |

Figure scripts write to their module's `figures/` folder; files are copied into the manuscript under
the names above (the renamed ones are given in brackets).

**Table 6** (`13_goldilocks_units.R`, `13b_goldilocks_regions.R`) reproduces an analysis from the
LD-pruning-paper repository on that repository's own Stage-1/Stage-2 partitions, and reads its
`cache_3sp/` (`map_only.rds`, `K_eco.rds`, `GTs_s1units.rds`, `stage2_Chr*_d1e+05.rds`; 58 MB) and
`results_3sp/3sp_cmin_loci.csv` from `LDSCNR_LDPP`. These inputs are deposited on Dryad with the
other intermediates.

## Licence and citation

Code: MIT (`LICENSE`). Cite the manuscript and this archive; cite Fang et al. (2021) and
doi:10.5281/zenodo.4722879 for the stickleback genotypes, Roberts Kingman et al. (2021,
*Sci Adv* 7:eabg5285) for the EcoPeaks, and Glazer et al. (2015, *G3* 5:1463) and Roesti et al.
(2013, *Mol Ecol* 22:3014) for the recombination map.
