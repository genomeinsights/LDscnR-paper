# Handoff — genomic control: DECIDED

**Status: settled 25 Aug 2026 by PK. Do not re-open, do not re-litigate.**

This note exists because the GC question was argued at length and the argument
went in several directions before landing. If you find reasoning elsewhere — in
`notes.txt`, in commit messages, in script comments — that points toward dropping
GC, calibrating by a surrogate-derived λ, or replacing the within-cell BH, that
reasoning is **superseded by this file**.

Applies to both repos: `~/gitlab/LDscnR-paper` and `~/gitlab/LDscnR`.

---

## The decision

**One rule, stated once, applied to every engine:**

> Apply genomic control when λ > 1.1. Test at α = 0.05.

It fires for LFMM (λ ≈ 4 on the simulations) and not for EMMAX (gif 1.02–1.08).
That is **one rule producing different outcomes because the engines differ in
calibration** — a property of the engines, not an inconsistency in the analysis.
Write it that way everywhere.

## Why (for your own understanding — do not re-derive it)

1. **It demonstrably works.** The headline PR-AUC result — the C-score
   outperforming plain α — was obtained under exactly this scheme.
2. **The bias runs the safe way.** λ_obs ≈ 4 and λ_perm ≈ 2–3, so per-dataset GC
   divides the observed by more than it divides the surrogates. Observed C-mass is
   depressed relative to null C-mass, s_R shrinks relative to the surrogate
   distribution, and **p_R is biased upward — conservative.** Power is lost, not
   validity.
3. **It is a rule a reader can follow.** No per-dataset tuning, no new machinery,
   consistent with standard practice in both engines' ecosystems.

---

## What to actually do

### Task A — harmonise the code so the rule is expressed once

`R/Parse_sim_data.R` currently states the rule twice, differently:

- **line 465–467** — LFMM, unconditional:
  ```r
  pv <- lfmm2.test(..., genomic.control = TRUE, full = TRUE)
  map[, lfmm_F := pv$fscores / pv$gif]
  ```
- **line 532–535** — EMMAX, gated at 1.1:
  ```r
  emx_gif = map[, median(emx_F) / qf(0.5, 1, nrow(GTs) - 2, lower.tail = FALSE)]
  map[, emx_F := if (emx_gif > 1.1) emx_F / emx_gif else emx_F, ]
  if (emx_gif > 1.1) map[, emx_p := pf(emx_F, df1 = 1, df2 = nrow(GTs) - 2, lower.tail = FALSE)]
  ```

Put both behind **one** gate at λ > 1.1, expressed once.

**Hard requirement: the numbers must not move.** LFMM's λ ≈ 4 is above the gate
and EMMAX's 1.02–1.08 is below it, so a correctly harmonised version reproduces
the current results exactly. Verify that before committing — regenerate at least
one replicate and diff the resulting `emx_p` / `lfmm_p` against the stored
version. If anything changes, the refactor is wrong, not the rule.

The same rule should apply in `module_sim_LDscnR/grm_comparison.R:71–73`, which
already implements the 1.1 gate — check it agrees.

### Task B — record λ so it can be reported

The Methods will carry a table of λ_obs and λ_perm per engine. λ_obs is available;
**λ_perm is not** — the null bundles store `C_surr` only, so per-surrogate
statistics are discarded after reduction.

- **Package (`~/gitlab/LDscnR`)**: store the per-surrogate inflation factor in the
  `ld_null` object. One number per surrogate, negligible cost, and it is expensive
  to retrofit once the B = 200 LFMM run has happened. Do this **before** that run.
- **Paper repo**: for the simulations, λ_perm ≈ 2–3 is already known and can be
  reported now.

### Task C — do NOT do any of the following

These were all discussed and rejected. They are recorded so you don't rediscover
them and start implementing:

- Do **not** switch to dividing by a surrogate-derived λ_perm instead of λ_obs.
- Do **not** remove GC from the LFMM arm.
- Do **not** build an empirical null from pooled surrogate statistics.
- Do **not** replace the within-cell BH with a count-calibrated threshold.
  (`module_sticklebacks_LDscnR/cscore_surrogate_threshold.R` is an existing
  prototype of this — leave it as a prototype, do not wire it into the pipeline.)
- Do **not** change α from 0.05.

Each of these is defensible in isolation and each carries its own cost; the
trade was made deliberately. The alternatives table belongs in the manuscript,
not in the code.

---

## What is still open — and it is a Discussion question, not a pipeline one

λ_perm ≈ 2–3 on **signal-free permuted** simulated data means a large share of
LFMM's inflation is engine miscalibration rather than biology. The Discussion
currently argues the opposite for the real data — "higher sensitivity, not
inflation", resting on λ ≈ 1.7 being insensitive to K over 1–12.

Two acceptable resolutions, and one of them must happen:

1. Compute λ_perm from the **within-locality** surrogates on the stickleback data
   and let the number decide; or
2. Soften the Discussion claim to what the current evidence supports.

**Caveat if you measure it:** LFMM estimates the latent factors *jointly* with the
environmental effect, so permuting the environment refits the model — part of any
λ_obs − λ_perm gap is a different estimation problem rather than signal. Use the
within-locality basis, which preserves between-locality environmental structure,
and check the fitted factors are comparable between observed and permuted runs.
This does not arise for EMMAX, where **K** is fixed across all runs.

---

## Two consequences for text already drafted

**The relative-power claim needs the uniform-rule framing.** The Discussion says
LFMM's region-detection PR-AUC is "marginally higher" than EMMAX's. Under the
uniform rule that comparison is legitimate, but the reader must be told that GC
fired for one engine and not the other, and why. Without that sentence it reads
as a confound.

**The EMMAX/LFMM calibration gap is a result, not a nuisance.** ~1.05 versus 2–3
under the null on identical simulated data: the explicit `w_j`-pruned kinship is
doing its job in this demography, the latent-factor correction is not. That
belongs in the Discussion alongside the structure-model-dependence argument.

---

## Where the fuller reasoning lives

The claude.ai project doc `claude/manuscript_notes.md`, §3, carries the decision,
the alternatives table with each option's cost, the caveats to state in the text,
and the provenance of every λ value quoted here. If you need to justify the choice
in prose, take it from there rather than reconstructing it.
