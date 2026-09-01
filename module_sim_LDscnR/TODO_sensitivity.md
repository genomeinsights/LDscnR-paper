# To do: sensitivity analyses, and the C-score reinterpreted

Status: planned, not started. LFMM replication runs first.

## The reframing

The C-score is dead **as a test statistic** — it beats BH at no deployable
threshold (development_history.pdf §2). But an earlier analysis in this project
found its significance layer to be *vacuous because it is identical to detection
stability*. That was recorded as a negative result. As a **sensitivity measure it
is exactly the right property**: C is the fraction of parameter settings under
which a region is called, which is what a stability analysis wants to report.

So C is retained, with its interpretation changed and its inferential claims
dropped:

- **not** a ranking statistic, **not** a significance layer, **no** threshold
  recommended;
- reported as: *of the settings examined, this region was called in a fraction C*;
- the only reading offered is that regions scoring **close to 1 are probably
  safe**. Nothing is claimed about regions scoring low — they may be genuine
  calls that one setting happens to miss.

This is deliberately weaker than anything the C-score was previously asked to do,
and it is the reading the evidence supports.

## Design: vary one at a time

Hold everything else at the canonical values and sweep each parameter alone, so
each region gets a stability score per parameter as well as an overall one.

| Parameter | Canonical | Sweep | What it controls |
|:--|:--|:--|:--|
| `ld_w_threshold` | 0.025 | to be chosen | which stage-1 clusters are flagged for eMLG treatment |
| `min_r2_rho` | 0.5 | to be chosen | the derived per-chromosome `min_r2` for stage-2 merging |
| minimum cluster size | 1 (none) | to be chosen | which units are eligible to be reported |

Report per region: C overall, and C per parameter so a region unstable in only
one dimension is distinguishable from one unstable in all three.

## Cautions carried from this module

- **Vary-one-at-a-time does not explore interactions.** If two parameters
  interact, the one-at-a-time profile understates instability. State this.
- **Do not report a threshold.** Every threshold this module tried on a
  C-like quantity failed, and the one that looked principled (τ = 0.05) was a
  coin flip.
- **Denominator discipline.** Report how many settings were examined alongside
  every C, since C = 1 over three settings is not C = 1 over thirty.
- **Stability is not validity.** A region called under every setting may still be
  a false positive; the sweep varies the analysis, not the truth.

## Also to do

- Same sensitivity scheme in the 3sp panel session, for the empirical analysis.
- LFMM replication of the main comparison (running).
- Close the LD-elevation-by-arm runs, which produced one line per cell and need
  re-checking.
