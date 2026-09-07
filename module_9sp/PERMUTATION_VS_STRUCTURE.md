# When permutation tells you the truth: signal vs. structure in 9sp vs. 3sp

**Status:** worked example for Supplementary materials + the Discussion's treatment of
permutation. Numbers below are all reproducible from committed scripts (paths given
throughout); nothing here is drawn from a scratch/ad hoc calculation that isn't also saved.

## The headline

3sp's permutation test worked cleanly: 97 significant consensus units, a permutation null
mean of 3.476, add-one p = 0.0080. 9sp's did not: 55 significant consensus units against a
permutation null mean of **67.36** (i.e. the null typically produces *more* apparent
discoveries than we observed), p = 0.3117.

The point of this write-up is that **9sp's null result is not a failure of the method — it
is the method correctly reporting that ecotype and population structure cannot be separated
in this dataset.** A structure-preserving permutation test is supposed to go silent exactly
when the phenotype of interest is too confounded with the structure it's permuted within to
say anything about it specifically. That is what happened. 3sp's permutation test being able
to speak at all is the thing that needed an explanation, not 9sp's being unable to.

## Why 9sp and 3sp differ

Both are stickleback panels run through the same pipeline (LD decay → Stage-1 LD clustering
→ EMMAX consensus/Simes → within-group phenotype permutation), but the population structure
underneath them is not comparable:

| | 3sp (*Gasterosteus aculeatus*) | 9sp (*Pungitius pungitius*) |
|---|---|---|
| Individuals | 117 | 149 |
| Populations / regions / lineages | a few regional localities, ecotype present in most | 30 populations, 6 localities, **4 lineages** |
| Ecotype ~ structure confound | modelled with no covariate | ecotype ~ lineage R² = **0.722** |
| Background LD | 0.055 | **0.301** |
| GRM off-diagonal (mean, SD) | tighter | mean −0.009 to −0.012, SD **0.26–0.32** |
| Marker-wise scan λ<sub>GC</sub> | **1.094** (healthy) | **0.868** (deflated) |

The lineage confound is severe enough that two of 9sp's four lineages (WL, WA — 43 of 149
individuals) are *entirely monomorphic for ecotype*: every WL and WA individual is
freshwater. There is no within-lineage ecotype contrast in 43% of the sample to permute in
the first place.

## The four diagnostic checks

All in `module_9sp/R/03d_diagnose_variance_components.R` (committed, rerunnable).

**1. The real tested phenotype hits a REML boundary.** `eco_resid` (ecotype, residualised on
lineage — see `EMMAX_COVAR` in `00_config.R`) is fit by REML as `vg=0.0560, ve=0.0000,
h²=1.0000` under the canonical GRM. That is a boundary solution: the model attributes *all*
phenotype variance to the genetic random effect and none to individual-level noise.

**2. That's not a generic numerical artefact.** Eight independent random phenotypes
(`rnorm(n)`, uncorrelated with the GRM by construction) all land at sensible interior
estimates under the identical GRM — median h² = 0.000. The boundary in check 1 is specific
to a phenotype that is actually structured like the kinship matrix, not a property of the
REML machinery itself.

**3. It is not unique to 9sp.** The identical check on 3sp's own real ecotype phenotype,
under 3sp's own canonical GRM, hits the same boundary: `vg=0.1565, ve=0.0000, h²=1.0000`.
So a population-differentiated binary phenotype hitting a REML boundary against a GRM built
to capture that same structure is not, by itself, a 9sp-specific pathology — it is a known
hard case for mixed models generally, and 3sp's own trusted, already-published-quality scan
sits in exactly the same regime.

**4. Where 9sp and 3sp actually diverge: genome-wide calibration.** The marker-wise EMMAX
scan's genomic inflation factor is healthy for 3sp (λ<sub>GC</sub> = 1.094, the mild
inflation expected from real polygenic signal) but **deflated** for 9sp (λ<sub>GC</sub> =
0.868) — i.e. 9sp's mixed-model correction is genuinely over-conservative genome-wide, not
merely sitting on a shared boundary quirk.

A leaner, more conventional GRM basis was tried as a mitigation: greedy LD-pruning
(19,452 markers) instead of the canonical Stage-1-pruned basis (753,625 markers).
It measurably improves the marker-wise λ<sub>GC</sub> (0.957) but does **not** rescue the
consensus-arm bottom line: 44 vs. 55 significant units, both still non-significant under
permutation (p = 0.2837 vs. 0.3117). A partial diagnostic improvement, not a fix — the
confound is in the phenotype-structure relationship itself, not merely in which markers
build the kinship matrix.

## A second, independent null agrees that something is off — for a different reason

`module_9sp/R/03c_EMMAX_structurednull.R` replaces the discrete within-lineage label
permutation with a *continuous* null: a surrogate phenotype drawn as `s ~ MVN(0, K)` (same
covariance as the real GRM), Gram-Schmidt-orthogonalised against the tested phenotype. Under
this null, the surrogate mean collapses to essentially zero (0.07/4,547 units for consensus,
0.03 for Simes) against the same 55/18 observed — flipping p from 0.31/0.75 to **0.001/0.005**.

This is not a "the structured null found the real signal" result to take at face value. Given
checks 1–4, a structured surrogate is drawn under the *same* over-corrected regime the real
data is fit in, and — being pure noise by construction — collapses to virtually nothing under
it. The real data, carrying genuine (if suppressed) signal, still produces some discoveries.
The two nulls disagreeing this violently (p=0.31 vs p=0.001) on the *same* observed data is
itself the diagnostic: it says the apparent significance is highly sensitive to which null
construction is used, which is a reason for caution, not a reason to report whichever number
looks better.

## What a user should actually look at

Distilled into a checklist — the numbers that tell you whether a result survives its null,
roughly in the order to check them:

1. **Observed vs. surrogate mean, and the permutation p-value.** The primary test:
   `ld_outlier_perm()`'s `observed`, `mean(surrogates)`, and `p`. 3sp: 97 vs 3.476, p=0.008.
   9sp: 55 vs 67.36, p=0.31.

2. **The null-to-observed ratio** (`mean(surrogates) / observed`, `ld_outlier_perm()`'s
   `realised_fdr` field — misleadingly named, it is *not* a false-discovery proportion, see
   `module_3sp`'s `00_config.R`/`11_table_sensitivity.R`). Read it as a magnitude, not just a
   pass/fail: low (single digits to low tens of %) means the observed count is a real
   multiple of what the null typically produces. Approaching or exceeding 100% means the null
   typically matches or beats the observed count outright — no count-level evidence of signal,
   regardless of the exact p-value. 3sp: 3.6–4.0%. 9sp: 122.5–328.4% (full sample),
   1257.4–674.8% (EL+Admixed subset — *worse*, not better, ruling out "small unbalanced
   permutation strata" as the explanation).

3. **Genomic inflation, λ<sub>GC</sub>, on the marker-wise scan.** Compute *before* trusting
   any p-value from the pipeline: `median(qchisq(1-p, df=1)) / qchisq(0.5, df=1)` over the
   genome-wide marker-wise EMMAX scan. Near 1, with mild inflation typical for real polygenic
   signal (3sp: 1.094) — trust the machinery. Substantially below 1 (9sp: 0.868) — the mixed
   model is over-correcting, and neither a "significant" nor a "non-significant" result from
   it should be taken at face value without understanding why. Substantially above 1 — under-
   correction / residual confounding, the opposite failure mode.

4. **The REML variance-component split for the tested phenotype**, `vg` and `ve` from the
   null-model fit (`emma.REMLE()` internally; exposed here by calling it directly on the
   tested phenotype and the canonical GRM). An interior solution (0 < h² < 1) is unremarkable.
   A boundary solution (h² → 1, ve → 0) is a flag to check further — as checks 2–3 above show,
   it is not disqualifying by itself (3sp hits it too), but it means the "genetic" random
   effect and the phenotype are close to indistinguishable, and downstream calibration checks
   (λ<sub>GC</sub> especially) become more important, not less.

5. **A cheap, pre-modelling check that would have flagged this before any of the above**:
   how well does population/lineage/locality membership predict the phenotype of interest?
   `summary(lm(phenotype ~ structure_covariate))$r.squared`. 9sp: R² = 0.722. A value this high
   is itself a warning that any downstream null-calibration exercise is going to be fighting
   an uphill battle, independent of which specific null construction is chosen.

6. **Agreement between independently-constructed nulls**, where feasible. A discrete label
   permutation and a continuous structure-matched null answering the same question should
   roughly agree in verdict. When they diverge sharply (9sp: p=0.31 vs p=0.001), that
   divergence is informative on its own — it says the significance verdict is unstable under
   reasonable changes to the null's construction, which argues for reporting the disagreement
   honestly rather than picking the more favourable number.

## What this is (and isn't) good evidence for

This is **not** evidence that 9sp has no real ecotype-associated genetic signal. The
consensus arm restricted to EL + Admixed (the two lineages with real within-lineage ecotype
contrast, `module_9sp/R/03b_EMMAX_sensitivity_ELAdmixed.R`) still finds 6 significant units,
and **all 6 are also significant in the full-sample result** — a small, repeatable core that
survives sample restriction even though the aggregate count doesn't clear either null. It is
evidence that, at this sample size and with this degree of structure, the standard toolkit
(EMMAX + GRM correction + structure-preserving permutation) cannot currently tell that signal
apart from a data-generating process built purely from relatedness. That is a real, useful,
and honestly-reported limit — precisely the kind of thing a permutation test is supposed to
be able to say, and precisely what 3sp's much more favourable numbers meant it never had to
say.
