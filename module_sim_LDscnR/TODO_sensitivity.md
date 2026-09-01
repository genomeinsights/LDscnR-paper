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

## Design: full factorial, few levels

**Superseded:** an earlier version of this plan proposed vary-one-at-a-time. PK
withdrew it — the runs are cheap, so there is no reason to give up orthogonality.
Full factorial with a small number of levels per parameter, rather than many
levels along one axis at a time.

| Parameter | Canonical | Levels (to be finalised) | Stage |
|:--|:--|:--|:--|
| `ld_w_threshold` | 0.025 | 0.0125, 0.025, 0.05 | stage 2 — needs re-clustering |
| `min_r2_rho` | 0.5 | 0.35, 0.5, 0.65 | stage 2 — needs re-clustering |
| minimum cluster size | 1 (none) | 1, 2, 5 | post-filter — free |

3 x 3 x 3 = 27 cells. C = the fraction of the 27 in which a region is called.

### Cost: 27 cells for the price of 9

Only the two stage-2 parameters change the partition; the size floor is applied
afterwards to an existing partition. So the grid needs **9 re-clusterings**, each
yielding 3 cells. At roughly 20 s per bundle that is about 4 hours on mini1 at
full width for all 80 panels — reduce the panel count first if a pilot is wanted.

### What factorial buys that one-at-a-time did not

- **Interactions.** `ld_w_threshold` and `min_r2_rho` both act on stage-2
  merging and are the pair most likely to interact; one-at-a-time cannot see it.
- **Marginal stability per parameter** is still available — average C over the
  levels of the other two — so nothing is lost relative to the earlier design.
- **A region unstable in one dimension** remains distinguishable from one
  unstable in all three, which was the point of reporting C per parameter.

## Cautions carried from this module

- **Report the factorial structure, not just C.** With 27 cells, a region at
  C = 0.5 called in every cell of one `min_r2_rho` level and none of another is a
  different object from one failing at random. The pattern is the result.
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
