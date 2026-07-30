# LDscnR-paper

Analysis code for the LDscnR paper: an LD-aware outlier-region (OR) pipeline
for genotype–environment association scans, benchmarked on simulated data with
known causal loci.

The pipeline builds on the [LDscnR](https://github.com/genomeinsights/LDscnR)
package for LD-decay estimation and per-SNP local LD support (`ld_w`), and adds
a paper-specific layer that clusters association signals into outlier regions,
scores their consistency across a filtering grid (the *C-score*), and — for
simulated data — evaluates detection performance (precision/recall, PR-AUC)
against the true QTNs.

## The `R/` folder

All analysis code lives in `R/`. It is sourced/run from the repository root.

| File | Role |
|------|------|
| `Parse_sim_data.R` | Parses raw simulation output into per-file analysis objects. Reads each simulated replicate, keeps one individual per population, computes MAF, flags true QTNs and their allelic effects, runs the EMMAX (`emx_p`) and LFMM (`lfmm_p`) association scans, and calls `LDscnR::compute_LD_decay()` / `LDscnR::compute_ld_w()` for the LD structure. Writes one `.rds` per replicate containing `GTs`, `map`, `env`, `LD_decay`, and `ld_ws`. |
| `define_ORs_functions.R` | The generic OR engine — works on any data, no ground truth required. FDR + local-LD (`th_ldw`) filtering of association p-values, LD-clustering of the retained SNPs into outlier regions (via an `r2`/distance edge list), and the C-score: how consistently each SNP/region is called across a grid of window size (`rho`), filtering threshold (`th_ldw`), minimum cluster size (`l_min`), and significance (`alpha`). |
| `Outlier_regions_simulation.R` | The simulation layer on top of the engine. Uses the known true QTNs to assign focal QTNs to outlier regions, compute TP/FP/FN, precision/recall/PR and PR-AUC, and provides the genome-wide summaries and plotting helpers (Manhattan, C-score, AUC trajectories). Source this **after** `define_ORs_functions.R`. |
| `analyse_sim.R` | Driver script. Runs the per-file pipeline across a set of replicates and produces (1) AUC-PR* trajectories comparing full LD-filtering against the unfiltered baseline for EMMAX and LFMM, and (2) genome-wide C-score and −log10(q) panels coloured by LD to the nearest true QTN. |
| `examples.R` | Companion driver. Produces a PR heatmap over the LD-filtering grid (`rho` × `th_ldw`, faceted by `l_min` × method), per-genome Manhattan plots (−log10(q) and `ld_w`), and genome-wide C-score plots. |

## Usage

Run from the repository root. A driver loads the engine and its dependencies:

```r
library(LDscnR)
source("R/define_ORs_functions.R")        # generic OR / C-score engine
source("R/Outlier_regions_simulation.R")  # simulation layer (ground truth, PR/AUC)

source("R/analyse_sim.R")                 # or R/examples.R
```

Set the input location and the file selection at the top of the driver
(`parsed_folder` and the file pattern) before running.

## Dependencies

- [LDscnR](https://github.com/genomeinsights/LDscnR) (LD decay, `ld_w`)
- `data.table`, `igraph`, `parallel`, `ggplot2`, `wesanderson`, `PRROC`, `patchwork`
- Parsing (`Parse_sim_data.R`) additionally uses `SNPRelate`, `LEA`, `factoextra`
