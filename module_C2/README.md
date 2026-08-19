# module_C2 — operating-grid stability (exploratory)

Focused methodological evaluation of the manuscript proposal currently written up as
*"A second-tier consistency score over the operating grid"* (`manuscript/ld_aware_outlier_regions.tex`,
paragraph at l.348, Eq. `cscore2`).

**Status: exploratory.** Nothing here is wired into `LDscnR` or into the manuscript
pipeline. All outputs are confined to `module_C2/`. Branch `explore-c2-grid-stability`.

Working name used throughout: **operating-grid stability**. The name
"second-tier C-score" is deliberately avoided during exploration — see
`notes/recommendation.md`.

## What is actually being measured — four distinct things

These are routinely conflated. They are not the same quantity and they do not
control the same error.

### 1. `ld_region_stability()` — *survival* stability (in the package)

`LDscnR/R/ld_region_stability.R`. For a fixed base region set, sweeps
`(tau_C, l_min)` and asks, in each cell, whether the region is still **called as a
region** — via `maxsz`, the size of the largest re-clustered region covering any of
its members, with survival iff `maxsz >= l_min`.

* Null enters only through `stability_null`, and only as a **per-cell gate**: a cell
  is "clean" when its genome-wide region-level FDR (mean surrogate region count /
  observed region count) is `<= fdr`. The gate is a property of the *cell*, identical
  for every region in it. It is **not** location-matched.
* Denominator is **all** grid cells.
* Answers: *is this region robust to the region-calling knobs?*

### 2. Proposed C2 — *significance* stability (this module)

For each cell, each observed region gets a **location-matched** empirical p-value
against the surrogates, BH-adjusted **within that cell** to `q_R`; the locus scores
the fraction of cells where it is significant.

* The null enters **per region and per location**, not as a blanket per-cell gate.
  This is the substantive difference from (1) — it is why a region can be *called*
  in a cell but not be *significant* there.
* Prototype denominator is the **usable** cells `|U|` only.
* Answers: *is this region's null-confirmation robust to the region-calling knobs?*

### 3. Genome-wide count FDR (`null_fdr()`, `calibrate_tauc()`)

`E[#null regions] / #observed regions` at a cell. A single genome-wide ratio used to
**pick an operating point**. It is not per-region, carries no location information,
and is the quantity the manuscript notes is too conservative here (ecotype is nested
in population, so a range-restricted true signal is charged to structure).

### 4. Location-matched empirical `p_R` and BH `q_R` (`region_empirical_pvals.R`)

The per-region test at **one** operating point. This is what actually carries formal
error control in the manuscript: BH over the tens of regions present at that point.
C2 is an aggregation *over* these, and inherits nothing automatic from them —
in particular, **C2 is not itself FDR-controlled**.

### The relationship

C2 (2) is best read as decomposable into detection (1) and null-confirmation (4):

```
detected in cell  (region-calling robustness)   ==  what ld_region_stability measures
        x
significant given detected  (null confirmation) ==  what the region p-value adds
        =
significant in cell                             ==  what C2 counts
```

`04_grid_sensitivity.R` / `02_compare_denominators.R` report these separately
(`D_l`, `Q_l`, `Q_l/D_l`) precisely because collapsing them hides which of the two a
low score comes from.

## Headline finding

**The significance layer is vacuous on these data.** Across 17 anchors x 175 cells x
4 matching rules there are **zero** cells where an anchor is detected as a region but
fails the location-matched significance test. `C2` is therefore numerically identical
to detection stability — what `ld_region_stability()` already computes. Same in the
simulations (`Q/D` = 1.00 true, 0.99 false). See `notes/recommendation.md`.

## Known defects in the prototype (`null_sig_landscape.R`)

Established in `R/01_reproduce_current_C2.R` (which reproduces the committed output
exactly, max |delta| = 0); see `results/01_reproduce.log`.

1. **The "loci" are chromosomes.** The script computes a 10 kb coordinate-consensus
   cluster id `cl` and then aggregates `by = chr`, never using `cl`. Every significant
   region on a chromosome is merged into one row, so the reported
   `Chr4:4.69-31.26 (Eda)` spans 26.6 Mb of chromosome 4 and the `(inversion)` /
   `(Eda)` annotations attach to whole chromosomes rather than to those loci.
2. **`B` is capped at 100** (`BCAP <- 100L`) although the `pop_perm` bundle carries
   **200** surrogates — an unnecessarily coarse p-value floor.
3. **Grid-defined hypotheses.** Loci are constructed *from* the significant cells, so
   the grid defines both the hypotheses and their stability. A locus that exists only
   in permissive cells cannot, by construction, be seen to fail in strict ones.
4. **Conditional denominator.** `|U|` (usable cells) hides how much of the grid is
   dead: `S = 1` is compatible with a single usable cell. On this dataset it is also
   uninformative -- every one of the 134 productive cells is usable, so `|U|` just
   counts cells that found anything.
5. **The grid omits its own operating point.** `0.05 %in% seq(0.02, 0.50, 0.02)` is
   `FALSE`, so the script's operating-point diagnostic evaluated to `character(0)` and
   silently printed nothing (0 occurrences in `sig_landscape.log`).

## Layout

```
R/00_helpers.R                 constants, loader (B auto-detected), clustering,
                               location-matched emp p + BH, anchor match rules
R/01_reproduce_current_C2.R    Q1  reproduce prototype; build shared grid cache
R/02_compare_denominators.R    Q2  S^U vs S^G vs availability vs D/Q decomposition
R/03_anchor_loci.R             Q3  fixed 17-region anchor set + marker/LD matching
R/04_grid_sensitivity.R        Q5  original vs coarse vs breakpoint grids
R/05_null_resolution.R         Q6  B subsampling, p floor, BH feasibility
R/06_make_figures.R            all figures
R/run_all.R                    driver
cache/                         git-ignored intermediates (grid_core.rds etc.)
```

`cache/grid_core.rds` holds the one expensive pass (per-tau observed regions **with
member markers** + all surrogate regions); every other script reads it.
