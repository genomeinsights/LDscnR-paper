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

## Two axes, not one: what C is a stability measure OF

Raised by the 3sp panel session and it changes the design. All three swept
parameters — `ld_w_threshold`, `min_r2_rho`, minimum cluster size — move the
**selection or the unit definition**. **None moves the test.** So C over that grid
measures how stable the *selection* is, and "called under every setting" reads as
though it were a statement about the finding.

The second association engine is exactly the missing axis:

| Axis | Varied | Held fixed | C measures |
|:--|:--|:--|:--|
| **Selection** | the 27-cell parameter grid | engine | stability of the units and filter |
| **Test** | EMMAX vs LFMM | units, filter, truth, budget | stability of the association test |

**Report them as two numbers, never as one.** Combined into a single C over 54
cells they are not interpretable, because a region can be perfectly stable on one
axis and unstable on the other and the scalar hides which. Reported separately
they answer different questions, and the pair is more informative than either.

On the panel side the selection axis is the one already known to be confounded —
a region called under every `ld_w` floor is stable because the filter keeps
finding it, not because the test keeps rejecting it. That confound does not exist
here, since the arms compared in this module share an identical unit set and
nothing is pre-filtered, which makes this the cleaner place to demonstrate
whether the two axes give different answers at all.

## Cautions carried from this module

- **Report the factorial structure, not just C.** With 27 cells, a region at
  C = 0.5 called in every cell of one `min_r2_rho` level and none of another is a
  different object from one failing at random. The pattern is the result.
- **Do not report a threshold.** Every threshold this module tried on a
  C-like quantity failed, and the one that looked principled (τ = 0.05) was a
  coin flip.
- **Denominator discipline.** Report how many settings were examined alongside
  every C, since C = 1 over three settings is not C = 1 over thirty — and also
  **whether they lie on a grid or a line**. C = 1 over ten settings varying one
  parameter is a much weaker statement than C = 1 over ten varying three, and the
  count alone does not distinguish them.

- **Name the axis.** State whether a reported C is selection stability or test
  stability. Unlabelled, "called under every setting" reads as a claim about the
  finding rather than about the analysis.
- **Firing rate is not impact.** A bug that fires rarely but *selectively on the
  heavy end of a distribution* is worse than one firing often and uniformly. The
  3sp panel session's sampler bug fired on 5.1% of draws and removed ~73% of the
  null's hits, because the classes where it fired were the large clusters
  carrying most of the mass. The question for any defect is not "can it fire" but
  "if it fired, would it fire where the mass is".

- **Print a bounded diagnostic next to its bound.** If any cell does matching or
  sampling, report a statistic that can come back *impossible* rather than merely
  disappointing. The 3sp panel session caught a live sampler bug because
  `max |lspan_case - lspan_control| = 4.619` against a caliper of 0.20 is
  arithmetically impossible for a correct matcher, which localised the fault with
  no judgement call. A standardised mean difference has no such bound and would
  have absorbed the same corruption as a plausible-looking 0.4.

- **Stability is not validity.** A region called under every setting may still be
  a false positive; the sweep varies the analysis, not the truth.

## Also to do

- Same sensitivity scheme in the 3sp panel session, for the empirical analysis.
- LFMM replication of the main comparison (running).
- Close the LD-elevation-by-arm runs, which produced one line per cell and need
  re-checking.
