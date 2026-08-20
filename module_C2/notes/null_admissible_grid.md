# A null-informed admissible operating grid

Second iteration of `module_C2`. Question: instead of choosing one operating point,
can the permutation null be used to **exclude** `(tau_C, l_min)` cells that
frequently manufacture spurious LD-supported regions, and can genomic signals then
be assessed across the remaining admissible grid?

Grid: `tau_C = seq(0.02, 0.50, 0.02)` x `l_min = {1,2,3,5,10,15,20}` = **175 cells**.
Null: population permutation, **B = 200** (all available; not capped at 100).
`tau_C = 0.05` is a **reference coordinate only** — it is not on the grid
(`0.05 %in% seq(0.02, 0.50, 0.02)` is `FALSE`) and never enters a denominator.
Every admissibility rule below is a function of null quantities alone.

---

## 1. Null behaviour over the grid (`results/null_grid_diagnostics.csv`)

`p_null_any` — the fraction of surrogates producing at least one region — ranges
**0.025 to 0.755** and is **never zero anywhere on the grid**. There is no
null-silent corner.

**The lenient axis is `l_min`, not `tau_C`.**

| | Spearman with `p_null_any` |
|---|---|
| `l_min` | **-0.970** |
| `tau_C` | -0.216 |

Across the entire `tau_C` range at `l_min = 1`, `p_null_any` falls only from 0.755 to
0.395 — it never approaches any candidate tolerance. Per `l_min` the ranges are:

| `l_min` | 1 | 2 | 3 | 5 | 10 | 15 | 20 |
|---|---|---|---|---|---|---|---|
| `p_null_any` | .40–.76 | .21–.48 | .13–.32 | .08–.18 | .055–.09 | .035–.085 | .025–.055 |

So the expected "lenient corner" is really a **horizontal band**: admissibility is
almost entirely a minimum-region-size requirement, with `tau_C` trimming only the
boundary rows.

**There is no natural break to cut at.** Sorted `p_null_any` has a median gap of
**0.0000** and a largest gap of 0.045 — a smooth continuum. Any tolerance is
therefore an arbitrary choice, not the identification of a real transition.

## 2. Candidate rules (`results/admissible_grid_rules.csv`, `figures/fig12`)

| rule | cells retained | `l_min` retained | max `p_null` retained |
|---|---|---|---|
| `p_null_any <= 0.20` | 116 / 175 | 3,5,10,15,20 | 0.200 |
| `p_null_any <= 0.10` | 86 | 5,10,15,20 | 0.100 |
| `p_null_any <= 0.05` | 42 | 15,20 | 0.045 |
| upper CP bound `<= 0.10` | 61 | 10,15,20 | 0.055 |
| **upper CP bound `<= 0.05`** | **0** | (none) | — |
| mean null coverage `<= 0.001` | 107 | all | 0.470 |
| mean null regions `<= 0.5` | 81 | 5,10,15,20 | 0.085 |

Every retained set is a single connected component and a contiguous high-`tau_C`
suffix within each retained `l_min`, so the shapes are well behaved.

**Finite-B resolution is the binding constraint.** At `eps = 0.05` the point-estimate
rule keeps 42 cells and the Clopper-Pearson upper-bound rule keeps **none** — not one
cell on the grid can be *shown* at B = 200 to have `p_null_any <= 0.05`. At
`eps = 0.10`, 25 of the 86 retained cells (29 %) are statistically indistinguishable
from failing. A cell at `k = 20/200` has CI [0.062, 0.150]: separating `eps = 0.05`
from `eps = 0.10` is beyond this null's resolution.

## 3. Consequence for the reported region set (`results/support_by_grid_threshold.csv`)

Detection support `D_r` = admissible cells recovering `r` / |G_adm| (denominator =
**all** admissible cells, including empty ones), `rec50` marker matching:

| grid | cells | regions with **zero** support | mean `D_r` |
|---|---|---|---|
| full | 175 | 0 / 17 | 0.087 |
| `P20` | 116 | **4 / 17** | 0.036 |
| `P10` | 86 | **9 / 17** | 0.027 |
| `U10` | 61 | **12 / 17** | 0.011 |
| `P05` | 42 | **17 / 17** | 0.000 |
| `C001` | 107 | **17 / 17** | 0.000 |

Tightening the tolerance does not merely re-rank the regions; past `eps = 0.10` it
**erases the entire reported region set**. The mechanism is mechanical: 9 of the 17
reference regions contain ≤ 4 SNPs, so once admissibility requires `l_min >= 10` they
can only ever be recovered by merging into a larger component.

Rank agreement with the full grid degrades correspondingly: Spearman 0.654 (`P20`),
0.477 (`U10`), 0.404 (`P10`); top-5 overlap 0.6 throughout.

## 4. The null rule does what it claims (`results/lenient_cell_validation.csv`)

Corroboration, computed **after** the rule was fixed and never used to build it:

* the 4 regions that lose all `P20` support are detected only in cells of mean
  `p_null_any` = **0.488**, and **0.0 %** of their detecting cells are in `P20`;
* the 13 regions retaining support are detected in cells of mean `p_null_any` =
  **0.361**, with 5–57 % of their detecting cells inside `P20`.

So the excluded cells are exactly where the lenient-only regions lived. Note the
diagnostic association is strong genome-wide (Spearman(`p_null_any`, observed region
count) = 0.925), which is why the empirical burden was kept out of the rule.

## 5. Detection vs BH significance

Across **all 19** grid x threshold configurations, the number of cells where a region
is detected but not BH-significant is **0**. Significance support is exactly
redundant with detection support and must not be presented as independent evidence.
This reproduces the first iteration's finding on the admissible grids too.

## 6. Matching rule and reference point

* `rec25` / `rec50` / `rec75` agree closely (Spearman 0.78 / 1.00 / 0.97 against
  `rec50`); the reciprocal variant gives 0.880. **The matching threshold is not a
  sensitive knob** — contrary to the first iteration's concern, which was driven by
  permissive "any-overlap" matching rather than by the retention fraction.
* Reference-point sensitivity is **material for the inventory, mild for the ranking**:
  neighbouring reference points yield 8–23 regions and 0–9 zero-support regions, but
  `Chr5:2.48–2.50 Mb` is the top-supported region at 4 of 5 reference points, with an
  identical `max D_r` = 0.121.

## 7. Anchor-free marker support

`S_m` correlates with the marker C-score at Spearman 0.69–0.85, so it is largely a
re-expression of `C_m` and is reported as *support*, never as new evidence.

On the admissible grids the anchor-free view finds **nothing the reference point
misses**: every locus family (r²-connected components of the support-gated markers,
500 kb gap-split — grouped by chromosome **and** component, not by chromosome alone)
overlaps a reference-point region, `n_novel = 0` at every support gate. Maximum `S_m`
outside the reference set falls to 0.069 under `P20` (0.269 on the full grid, where
41 "novel" families appear but with a median of **1** marker each — singleton noise).
Region- and marker-level support agree at Spearman **0.992** under `P20`.
