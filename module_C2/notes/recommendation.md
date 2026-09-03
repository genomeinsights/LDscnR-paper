# Operating-grid stability — findings and recommendation

Exploratory evaluation of the manuscript's *"second-tier consistency score over the
operating grid"* (`ld_aware_outlier_regions.tex` l.348, Eq. `cscore2`).
Branch `explore-c2-grid-stability`; all work under `module_C2/`.

Dataset: 3sp sticklebacks, `pop_perm` null, **B = 200** (the prototype capped at 100;
the bundle carried 200). Grid `tau_C = seq(0.02, 0.50, 0.02)` x
`l_min = {1,2,3,5,10,15,20}` = 175 cells. Anchors: the 17 primary EMMAX regions at
`tau_C = 0.05, l_min = 3, rho_ld = 0.60`, verified identical to the committed
`regions_tau0.05_lmin3_rho0.60.csv` (17/17, exact coordinates and sizes).
Simulations: Nemo V2_c1, env1–5, spatial null B = 100.

---

## The headline result

**On both datasets, the location-matched significance layer never fires
independently of detection.** Across all 17 anchors x 175 cells x four matching
rules, the number of cells where an anchor is *detected as a region but not
location-matched significant* is **zero** (`results/03_anchor.log` [3b];
`figures/fig2_anchor_grid_tiles.png` has no "detected only" tiles). In the
simulations the same ratio `Q/D` is 1.00 for true regions and 0.99 for false ones.

So `C^{(2)}` as currently defined is **numerically identical to detection
stability** — i.e. to what `ld_region_stability()` already computes — on the data it
was proposed for. The "second null layer" that motivates it contributes nothing.

The mechanism (Q4/Q6, `results/05_null.log`):

* the `pop_perm` null rarely rebuilds any observed locus, so **27.8 %** of all
  region-cells sit exactly at the p-value floor `1/(B+1) = 0.00498`;
* BH is a step-up procedure, so a set of `k` regions tied at the floor is rejected
  *en bloc* as soon as `p0 <= alpha·k/m`. Ties rescue each other;
* where `l_min` is large, cells test only 1–3 regions, and BH is trivially permissive
  (`m_median = 1` at `l_min >= 10`);
* either way significance is near-automatic: among detected cells, the *raw*
  empirical p is already `< 0.05` in **100 %** of cases and BH demotes **none**
  (`anchor_bh_behaviour.csv`: `frac_raw_p_sig = frac_BH_sig = 1`, `lost_to_BH = 0`).

Note this corrects a bound I initially derived: a naive "distinct-p" argument says a
cell testing `m > alpha(B+1) = 10.1` regions cannot reject anything. That is wrong
under ties — all 32 such cells do produce significant regions.

---

## Defects found in the prototype (`null_sig_landscape.R`)

The prototype reproduces **exactly** (max |Δ n_cells_sig| = 0 against the committed
`region_stability2.csv`), so these are definitional, not numerical.

1. **The reported "loci" are chromosomes.** The script computes a 10 kb consensus
   cluster id `cl` and then aggregates `by = chr`, never using `cl`. Hence
   `Chr4:4.69-31.26 (Eda)` — a 26.6 Mb span of chromosome 4 labelled as the *Eda*
   locus, whose actual anchor is `Chr4:12.809–12.812 Mb` (4 SNPs). The intended
   `by = chr, cl` merge gives **126** loci, not 20.
   Consequence: the published ranking is essentially a per-chromosome C-mass
   statistic — Spearman(stability, max C-mass) = **0.988** — and it cannot separate
   the 7 of 17 anchors that share a chromosome with another anchor.
2. **The grid does not contain the operating point.** `seq(0.02, 0.50, 0.02)` steps
   over `tau_C = 0.05`; `0.05 %in% seq(0.02, 0.50, 0.02)` is `FALSE`. The script's
   `op <- res[tau == OP_TAU & lmin == OP_LMIN]` is therefore empty, and the
   `sprintf` reporting it returns `character(0)` — the whole "MAX-significant cell
   vs a-priori operating point" diagnostic **silently vanished** (0 occurrences in
   `sig_landscape.log`). The heatmap's `x` marker is drawn off the panel.
3. **`B` capped at 100** although 200 surrogates were available.
4. **Grid-derived hypotheses.** Loci are built *from* the significant cells, so the
   grid defines both the hypotheses and their stability; a locus that exists only in
   permissive cells cannot be observed to fail in strict ones.

Fixed-anchor scores are **uncorrelated** with the prototype's published ranking under
the recommended matching rule (Spearman **+0.065**; +0.372 under permissive matching).
The prototype's ranking should not be carried into the manuscript in any form.

---

## Answers to the twelve questions

**1. Is operating-grid stability scientifically useful?**
Partly. As a *descriptive* statement of how sensitive a reported region is to
`(tau_C, l_min)` it is useful and honest, and the simulations show it separates true
from false regions (true mean `Q` 0.437 vs false 0.094; Wilcoxon `p <= 7e-4` in all
four informative envs). As a *second-tier significance* construct it is not useful on
these data, because the significance layer is degenerate.

**2. Retain, revise, or drop?**
**Revise and demote.** Drop the "second-tier C-score" framing and the `C^{(2)}`
notation; keep a plainly-named descriptive robustness statistic. Do not present it as
a nested analogue of the marker C-score — the marker C-score integrates a genuinely
varying candidacy criterion, whereas this one integrates a criterion that is constant
given detection.

**3. Recommended definition.**
> The **operating-grid stability** of a *prespecified* region `l` is the fraction of
> prespecified `(tau_C, l_min)` cells in which `l` is recovered as a region and that
> region passes the location-matched region-level significance criterion
> (`q_R < 0.05`):
> `S_l = |{g in G : l recovered and significant in g}| / |G|`,
> reported together with the grid `G`, the null basis, `B`, and the matching rule.

It is a robustness ranking. It is **not** a p-value and confers **no** FDR control
across the grid.

**4. Denominator: all cells or usable cells?**
**All cells (`|G|`).** Two reasons. (a) On this dataset the conditional denominator is
a pure rescaling — `S^U == S^G · |G|/|U|` exactly, for every locus (verified), so it
changes reported magnitudes (by 1.31x here) without changing any ranking, while
making `S = 1` mean "significant wherever anything was found" rather than
"significant everywhere". (b) `|U|` turns out to carry no information about null
danger: **every one of the 134 cells that produces any region produces a significant
one**, so `U` is just the set of productive cells and `|U|/|G| = 0.77` measures only
how much of the grid is too strict to find anything. If availability is of interest,
report `A = |U|/|G|` separately, never fold it into the denominator.

**5. How should loci be fixed and matched?**
Fix loci **before** the grid search, as the prespecified region set at the primary
operating point, retained as **member-marker vectors** (not coordinate spans).
For matching, **require `|A ∩ R| / |A| >= 0.5`** ("recover ≥ 50 % of the anchor's
markers"). Rationale: `any`, `overlap-coefficient >= 0.5` and physical-span overlap
are **numerically identical** here (all give 700 detected cells) because they accept a
match in both directions — a called region that has swallowed the anchor into a
chromosome-scale component counts just as much as one that reproduces it. That
permissive behaviour is exactly what reproduces the prototype's chromosome-level
statistic (Spearman vs prototype +0.372 for all three, +0.065 for `recover`). The
`recover` rule instead asks whether the *anchor itself* is still there, which is the
biologically meaningful question and the only rule of the four that is not silently
chromosome-level.
Caveat: the rule matters a lot — Spearman between `recover` and the permissive rules
is 0.806 with only 3/5 top-5 overlap, and the top-ranked locus changes
(`Chr4:27.55` vs the Chr1 inversion). Any published ranking must state its rule.
Reassuringly, chaining ambiguity was **not** a practical problem: mean matched
regions per detected cell = 1.00 for 16 of 17 anchors.

**6. Should detection and significance be separate?**
Yes — and on these data reporting both is what exposes the degeneracy. Report `D`
(detection) and `Q` (significance) side by side; `Q/D` is worth computing precisely
because a value pinned at 1.0 is the diagnostic that the null layer is inactive.

**7. Is the score stable to reasonable grid changes?**
**No.** This is the second-most important negative result.
* Restricting to `tau_C > 0.26` sends **every** anchor to `S = 0` (`G5`); restricting
  to `tau_C <= 0.26` reproduces the full-grid ranking **exactly** (Spearman 1.000).
  All information lives in the permissive half; the strict half contributes only
  denominator, so lengthening the `tau_C` range mechanically dilutes every score.
* A coarser but entirely reasonable grid (5 `tau_C` x 5 `l_min`) drops 13 of 17
  anchors to exactly 0 and gives Spearman 0.537, top-5 overlap 0.6.
* The breakpoint grid changed nothing, because at this resolution **all 25** `tau_C`
  values alter the observed region set — there is no redundancy to remove.
* Sparse sampling of the same range: Spearman 0.894.
Named loci move accordingly: *Eda* ranks 10th (original), 4th (coarse), 12th
(breakpoint); the Chr1 inversion 4th, 4th, 5th; small 3–4 SNP anchors occupy ranks
12–16 under the original grid and collapse to ties elsewhere.

**8. Is the permutation count adequate?**
For the **ranking**, yes trivially: subsampling to B = 25 leaves Spearman >= 0.999 and
top-5 overlap 1.00. But that stability is itself a symptom — the ranking is carried by
detection, not by the null. For the **p-values**, B = 200 is thin: the floor is
0.00498 and 27.8 % of region-cells are tied at it, so those regions are
rank-indistinguishable no matter how strong they are. Use all 200 (not 100). If the
significance layer is ever to do real work, B would need to grow by an order of
magnitude.

**9. Does it separate true and false regions in simulations?**
Yes, but it does not beat a trivial baseline. Replicate-averaged over 5 envs
(`sim_validation_prauc.csv`):

| ranking | PR-AUC (mean ± SE) |
|---|---|
| significance stability `Q` | 0.901 ± 0.052 |
| detection stability `D` | 0.901 ± 0.052 |
| **max C-score `maxC`** | **0.943 ± 0.025** |
| region size | 0.807 ± 0.088 |
| `Q/D` | 0.392 ± 0.166 |

`Q` and `D` are identical to three decimals in every env. A single number requiring no
grid at all — the region's maximum marker C-score — ranks true regions **better**.
Caveats: env5 is degenerate (1 anchor, base rate 1.0) and inflates all means;
excluding it, `Q` = 0.876 and `maxC` = 0.929, so the ordering is unchanged.

**10. Main manuscript, supplement, or development notes?**
**Supplement**, as a parameter-sensitivity analysis — not the main text, and not as a
second inferential tier. The main text should keep the single prespecified operating
point with its location-matched `q_R` as the reported inference. If a robustness
ranking is wanted in the main text, `ld_region_stability()`'s detection stability is
the honest version and is already implemented; the significance wrapper adds nothing
here. Do not migrate into `LDscnR` — the package function already covers what this
measure actually computes on these data.

**11. What claims are safe?**
* High stability means the locus is insensitive to the stated `(tau_C, l_min)` range.
* Low stability means the locus sits near a calling boundary — **not** that it is false.
* The statistic exposes parameter sensitivity and reduces reliance on one chosen cell.
* Formal error control comes from the region-level empirical tests at the operating
  point, not from the stability fraction.
* Reporting requires the grid, null basis, `B`, and matching rule; the number is
  meaningless without them.

**12. What would be overstated?**
* That grid cells are independent, or that the fraction is a probability. They are
  strongly dependent — `tau_C <= 0.26` alone reproduces the full ranking exactly.
* That integrating over the grid solves the multiple-testing problem, or confers FDR
  control. It does not; and on these data it does not even add a second criterion.
* That the score is invariant to grid design. A defensible coarser grid zeroes 13 of
  17 loci.
* That `S = 1` is strong evidence. Under the conditional denominator it can mean one
  usable cell.
* That physical overlap identifies two regions as the same LD signal.
* That it validates the method in simulation *beyond* what `maxC` already achieves.

---

# Iteration 2 — the null-informed admissible grid

Full detail in [`null_admissible_grid.md`](null_admissible_grid.md). Scripts
`R/08`–`R/13`. This iteration asked whether the null can *exclude* cells rather than
the analyst *choosing* one, and then assess signals across what remains.

## Answers to the eight iteration-2 questions

**1. Is there a defensible null-derived admissible grid?**
A rule can be written down cleanly and applied honestly, but it does not identify a
real feature of the data. `p_null_any` is a **smooth continuum** from 0.025 to 0.755
with a median gap between sorted values of 0.0000 and no break anywhere; it is
**never zero** on any of the 175 cells. So there is no natural boundary to discover —
only a tolerance to impose. The rule is also not the "lenient corner" that motivated
it: null-region production is governed by `l_min` (Spearman -0.970), not `tau_C`
(-0.216), so admissibility is essentially a minimum-region-size requirement wearing a
different name.

**2. Is the boundary stable at B = 200?**
**No, at the tolerances of interest.** At `eps = 0.05` the Clopper-Pearson
upper-bound rule retains **zero cells**: not one cell can be *shown* at B = 200 to
have `p_null_any <= 0.05`. At `eps = 0.10`, 29 % of retained cells are statistically
indistinguishable from failing. Separating `eps = 0.05` from `eps = 0.10` is beyond
this null's resolution (a cell at 20/200 has CI [0.062, 0.150]). Only `eps = 0.20` is
reasonably resolved (13 % ambiguous).

**3. Does excluding null-prone cells materially improve interpretation?**
It materially *changes* interpretation; whether that is an improvement is not
established. Zero-support regions go 0/17 (full) → 4/17 (`P20`) → 9/17 (`P10`) →
12/17 (`U10`) → **17/17** (`P05`, `C001`). Past `eps = 0.10` the filter does not
sharpen the region set, it deletes it — mechanically, because 9 of the 17 reference
regions have ≤ 4 SNPs and admissibility demands `l_min >= 10`. The genuine gain is at
`eps = 0.20`, where the 4 regions that drop out are demonstrably the ones living in
null-prone cells (mean `p_null_any` 0.488 vs 0.361, with 0 % of their detecting cells
admissible). That validation is real and was computed after the rule was fixed.

**4. Which signals are robust across reasonable tolerances?**
Three, and only three: **`Chr5:2.48–2.50 Mb`** (top-supported at 4 of 5 reference
points), the **Chr1 inversion**, and **`Chr20:0.18–0.26 Mb`**. These are the only
regions retaining support at `U10`, the only stable-core families in the anchor-free
analysis, and the only ones whose rank is unaffected by the tolerance. `Eda` survives
`P20` (rank 10 → 5) but not `P10`. Everything else is tolerance-dependent.

**5. Does BH significance add anything beyond detection?**
**No.** Zero detected-but-not-significant cells across all 19 grid × threshold
configurations, extending iteration 1's result to every admissible grid. Do not
present significance support as independent evidence.

**6. Is the anchored Manhattan an adequate presentation device?**
Yes, with one caveat now fixed: region support and marker support live on different
scales (max 0.12 vs 0.65), because a region must be recovered ≥ 50 % intact while a
marker need only fall in any region. A shared colour scale flattens the anchored
panel, so both a shared-scale version (for side-by-side reading) and an own-scale
version are provided. The figure is now genuinely per-region, unlike the legacy
`region_c2_manhattan.png`, which is per-chromosome.

**7. Does the anchor-free analysis reveal signals the reference point omits?**
**No.** On every admissible grid, every locus family overlaps a reference-point
region (`n_novel = 0` at all support gates), region- and marker-level support agree
at Spearman 0.992, and the maximum `S_m` outside the reference set falls to 0.069.
The 41 "novel" families on the full grid have a median of **1** marker — singleton
noise. The reference point is not hiding anything.

**8. Hard rule, sensitivity across tolerances, or no filtering?**
**Show sensitivity across tolerances; do not adopt a hard rule.** Admissibility
depends entirely on an arbitrary tolerance (no break in the null continuum), the
tolerance most people would reach for (`eps = 0.05`) is unresolvable at B = 200 and
annihilates the region set, and the axis actually being filtered is `l_min` — which
is already an explicit, reportable parameter of the method. Recommended framing: keep
the single prespecified operating point with its location-matched `q_R` as the
reported inference, and present the null landscape (`fig7`) plus support at
`eps ∈ {0.20, 0.10}` as a supplementary sensitivity analysis, stating plainly that
only three regions are tolerance-robust.

## What changed from iteration 1

* Iteration 1 concluded the matching rule dominated the ranking. That was driven by
  permissive **any-overlap** matching; among genuine retention thresholds
  (`rec25/50/75`) agreement is 0.78–1.00, so the threshold is **not** a sensitive
  knob. The distinction that matters is retention-based vs overlap-based matching.
* Iteration 1's headline (significance ≡ detection) is confirmed and extended.
* The `|U|` conditional denominator is again shown to be a pure rescaling
  (Spearman 1.000 with the full-grid score).
