# Handoff — sweep the GRM construction to bring lambda into [1, 1.1]

For the session that regenerates simulation bundles (`module_sim_LDscnR/regen_sim_data.R`).
Written 2026-08-25 against `regen_sim_data_nobgs`.

## Read this first: it is NOT the LD-clustering thresholds

The request was phrased as "sweep the LD-clustering thresholds so lambda is
appropriate". **LD-clustering cannot change lambda**, and sweeping `rho_ld` or
`dcap` will produce a flat, uninformative surface.

Lambda is fixed at `regen_sim_data.R:80`:

```r
gif <- stats::median(map$emx_F) / stats::qf(0.5, 1, n - 2, lower.tail = FALSE)
```

That is a property of the **EMMAX scan**, which depends on the **GRM**. Clustering
(`ld_edges` / `ld_regions`, governed by `rho_ld` and `dcap`) happens strictly
downstream of the p-values — verified: no clustering parameter appears anywhere
in `regen_sim_data.R`. Changing how markers are grouped into regions cannot move
a statistic computed before grouping happens.

**The knob that does move lambda is GRM marker selection**, `regen_sim_data.R:43-47`.
That is genuinely LD-based, which is probably where the phrasing came from.

## Baseline: only 42 of 100 files are in target

Measured on the current `regen_sim_data_nobgs` bundles (`emx_gif`, i.e. lambda
BEFORE genomic control), 10 replicates x 10 env:

| | |
|---|---|
| in [1, 1.1] | **42 / 100** |
| deflated (< 1) | **35** |
| inflated (> 1.1) | 23 |

By replicate (median lambda, and how many of its 10 env cells are in target):

| replicate | median | in target | range |
|---|---|---|---|
| chr8 | 0.999 | 1/10 | 0.848 – 1.167 |
| chr6 | 1.002 | 5/10 | 0.967 – 1.094 |
| chr3 | 1.008 | 5/10 | 0.933 – 1.104 |
| chr2 | 1.020 | 2/10 | 0.896 – 1.226 |
| chr1 | 1.022 | 5/10 | 0.952 – 1.103 |
| chr7 | 1.025 | 5/10 | 0.928 – 1.163 |
| chr10 | 1.027 | 8/10 | 0.988 – 1.055 |
| chr4 | 1.060 | 5/10 | 0.918 – 1.159 |
| chr9 | 1.067 | 5/10 | 0.912 – 1.224 |
| **chr5** | 1.110 | **1/10** | **0.629 – 1.526** |

chr5 is the outlier and is expected to be: its decay `b` is 0.0585 against ~0.029
everywhere else — roughly twice the background LD, consistent with long
low-recombination stretches. It also has the smallest GRM marker sets
(~7,900–8,250 vs ~10,200 typical). Do not tune the global default to rescue chr5;
report it separately.

## Two things that change what you should target

**1. The upper bound is already enforced; the lower bound is not.**
`regen_sim_data.R:81` applies genomic control when `gif > 1.1`:

```r
if (gif > 1.1) { map[, emx_F := emx_F / gif]
                 map[, emx_p := stats::pf(emx_F, 1, n - 2, lower.tail = FALSE)] }
```

So the 23 inflated files are already corrected downstream, and `emx_gif` records
the pre-correction value. **Deflation below 1 is neither corrected nor flagged** —
it silently costs power, and 35 files have it. Prioritise the lower bound.

**2. Do not chase the marker count.** Across the 100 files,
`cor(lambda, fraction of markers kept in the GRM) = -0.018` — no relationship.
The GRM currently keeps a median 33% of markers (range 26–46%). Keeping more or
fewer will not fix this; *which* markers are kept is what matters.

## What to sweep

All in `regen_sim_data.R:43-47`:

```r
GRM_METHOD <- Sys.getenv("SIM_GRM", "complexity_chain")   # or "ld_w_threshold"
CR_RHO     <- 0.5                                          # ld_complexity_reduction rho
PRUNE_ARGS <- list(ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
                   score_threshold = 0.80, ...)
GRM_LDW_THRESHOLD <- 0.02                                  # ld_w_threshold method only
```

Suggested grid, but adapt it — the point is coverage of both arms, not these
exact values:

- `GRM_METHOD`: `complexity_chain` and `ld_w_threshold` (both arms; the README's
  earlier comparison found chain gif ~1.04 vs ld_w<b gif ~0.93, so the second arm
  is the deflating one and should be included precisely to confirm that)
- `CR_RHO`: 0.3, 0.4, 0.5, 0.6, 0.7
- `PRUNE_ARGS$ld_w_threshold`: 0.015, 0.025, 0.04
- `PRUNE_ARGS$score_threshold`: 0.7, 0.8, 0.9
- `GRM_LDW_THRESHOLD` (second arm only): 0.01, 0.02, 0.03

You do not need to re-run the full regeneration per cell. Lambda needs only:
GRM markers -> `snpgdsGRM` -> `emmax` -> the one-line gif. Decay, `ld_w` and LFMM
are unaffected by these parameters and can be reused from the existing bundles.
That makes the sweep cheap; a full regeneration per grid point would not be.

## How to report

- Lambda **before** genomic control, always. The stored `emx_gif` is already this;
  keep it that way.
- **Replicate-average over all 10 env cells, mean +/- SE.** A threshold chosen on
  one env will fluke — this is a standing rule in `module_sim_LDscnR/README.md`.
- Per grid point: fraction of files in [1, 1.1], fraction below 1, fraction above
  1.1, and median lambda. Break chr5 out separately.
- Prefer the setting that **minimises deflation**, then narrows spread. A setting
  with median 1.02 and range 0.98–1.06 beats one with median 1.00 and range
  0.63–1.53.

## The trap to check for

`module_sim_LDscnR/README.md` records it, and it matters more than the lambda
target itself:

> gif ~ 1 is a good check against *global* deflation (the sim failure mode) but
> blind to *regional* signal absorption (the 3sp failure mode).

You can drive lambda to exactly 1 by absorbing the signal into the GRM, which
destroys power while looking perfectly calibrated. So for every candidate
setting, also report:

- how many **true QTNs** are in `grm_markers` (should be none, or as few as possible);
- how many markers within ~500 kb of a QTN are in `grm_markers`;
- the number of single-SNP BH hits at q < 0.05, as a power proxy — if lambda
  improves while hits collapse, the GRM is eating the signal.

A setting that reaches [1, 1.1] by absorbing QTNs is worse than the status quo.

## Do not

- Do not sweep `rho_ld` or `dcap` for this. They cannot affect lambda.
- Do not change `regen_sim_data.R`'s defaults until the sweep is
  replicate-averaged and the QTN-absorption check is clean.
- Do not tune the global default to fix chr5. Two outliers in twenty is expected
  from its LD structure; report it as a caveat instead.
