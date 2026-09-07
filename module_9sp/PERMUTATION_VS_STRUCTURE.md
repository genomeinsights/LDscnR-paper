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

## 3sp vs. 9sp, side by side

Both are stickleback panels run through the identical pipeline (LD decay → Stage-1 LD
clustering → EMMAX consensus/Simes → within-group phenotype permutation), so every row below
is a like-for-like comparison, not a difference in method.

| | 3sp (*G. aculeatus*) | 9sp (*P. pungitius*) |
|---|---:|---:|
| **Sample & structure** | | |
| Individuals | 117 | 149 |
| Populations | a few, regional | 30 |
| Regions / localities | 4 | 6 |
| Lineages | — | 4 (Admixed, EL, WA, WL) |
| Phenotype ~ structure covariate R² | — (no covariate used) | 0.722 (ecotype ~ lineage) |
| Background LD | 0.055 | **0.301** |
| GRM off-diagonal: mean (SD) | −0.0098 (**0.0847**) | −0.009 to −0.012 (**0.26–0.32**) |
| **Model calibration** | | |
| Marker-wise scan λ<sub>GC</sub>, canonical GRM | **1.094** (healthy) | **0.868** (deflated) |
| Marker-wise scan λ<sub>GC</sub>, leaner (greedy) GRM | not run | 0.957 (improves, doesn't fix result) |
| REML h² for the real tested phenotype | 1.0000 (boundary) | 1.0000 (boundary) |
| REML h², median over 5–8 random phenotypes | 0.166 | 0.000 |
| **Consensus arm** | | |
| Tested units | 1,356 | 4,547 |
| Observed significant | **97** | **55** |
| Permutation surrogate mean | 3.476 | **67.36** |
| Permutation p | **0.0080** | 0.3117 |
| Null-to-observed ratio | **3.6%** | **122.5%** |
| **Simes arm** | | |
| Observed significant | **74** | **18** |
| Permutation surrogate mean | 2.95 | **59.12** |
| Permutation p | **0.0149** | 0.7512 |
| Null-to-observed ratio | **4.0%** | **328.4%** |

Cells marked "not run" reflect what has actually been computed, not an assumption that 3sp
would look the same — 3sp has never needed these diagnostics run on it, which is itself part
of the contrast.

### 9sp-only follow-ups (no 3sp analogue exists to compare against)

These two checks were built specifically to chase down *why* 9sp's numbers above look the way
they do; 3sp was never confounded enough to need either one.

| | Consensus | Simes |
|---|---:|---:|
| **EL+Admixed subset** (106 of 149 individuals, own GRM) | | |
| Observed significant | 6 | 10 |
| Surrogate mean | 75.44 | 67.47 |
| Permutation p | 0.7762 | 0.7015 |
| Null-to-observed ratio | **1257.4%** | **674.8%** |
| **GRM-structured null** (full 149, MVN(0,K) surrogate) | | |
| Observed significant | 55 | 18 |
| Surrogate mean | 0.07 | 0.03 |
| p | **0.0010** | **0.0050** |

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
