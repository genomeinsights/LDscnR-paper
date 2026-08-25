# Nemo sim pipeline — raw tarballs to analysis bundles

One script turns a raw Nemo `.tgz` into the bundle `run_sim_LDscnR.R` reads, and
(optionally) a set of population-genetic / local-adaptation / background-selection
summaries. It replaces the old two-step route (`Parse_sim_data.R` → `regen_sim_data.R`).

## Files

| file | what it is |
|---|---|
| `parse_and_regen_sim_data.R` | the pipeline; one file per invocation, or a whole cell with `all` |
| `run_batch.sh` | one tag × one subsample fraction, parallel over (V, c, env) cells |
| `run_all.sh` | all three fractions, sequentially |
| `run_sim_nulls.R` | surrogate/permutation nulls for one cell (5 null types × 2 engines) |
| `run_nulls.sh` | drives the nulls over a set of cells, one stage at a time |

## Running on another machine

Only two things are machine-specific:

```bash
export SIM_ROOT=/path/to/Nemo_sim          # everything else defaults relative to this
export LDSCNR_PATH=/path/to/LDscnR         # omit if LDscnR is installed as a package
```

Needs R with `data.table`, `SNPRelate`, `LEA`, `igraph`, `LDscnR` (plus `devtools`
only if loading LDscnR from source). **LDscnR must be at commit `5da812a` or later** —
earlier versions have no `ld_complexity_reduction(gds =)` and cannot rebuild edge
lists, which the GRM step relies on. Use `757b835` or later if you can: it cuts peak
memory per process from 9.4 GB to 7.3 GB, which is what decides how many workers fit. The script checks this at startup and stops with
instructions rather than failing part-way through the first file. Note that an
*installed* LDscnR may well be older than your source checkout; `LDSCNR_PATH` is the
safer choice. Then:

```bash
SIM_CELLS=V0.5_c2,V1_c1.5,V2_c1 ./run_all.sh 12     # the reduced scope: 900 files, ~1 h
./run_batch.sh nobgs 0.75 12                        # or one stage at a time
```

`SIM_CELLS` restricts to a comma-separated list of (V, c) cells; omit it for all nine.
`SIM_STAGE` (default `$HOME/sim_stage`) is a local scratch folder on the **internal**
SSD: bundles are written there and moved to the volume only once a file is complete,
so the slow external leg is one sequential move rather than the whole write. The bundle
is moved last, so a killed run cannot leave a half-written file that the resume logic
mistakes for finished.

Both are **resumable** — an existing bundle is skipped, so re-run after any
interruption. Expect roughly 3–6 min per file, 900 files per tag per fraction.

## Parallelism and memory

Peak resident memory is ~7.3 GB per process, nearly all of it inside
`compute_LD_decay`. Note that this figure overstates what is actually resident: R's
high-water mark for live objects is about half of it, because the allocator holds
freed pages rather than returning them. Six workers ran fine on a 39 GB machine
despite 6 x 7.3 GB exceeding it.

On 48 GB, 12 workers is roughly 2.3x oversubscription — more than has been tested.
Start there but check within the first ten minutes:

```bash
vm_stat 60          # watch Swapins / Swapouts
```

If swapping climbs steadily, `pkill -f parse_and_regen_sim_data.R` and restart with
8 workers; nothing is lost, the run resumes. On **Linux** start at 8 rather than 12:
memory accounting is stricter there and the OOM killer terminates workers outright
instead of compressing.

## Subsampling

`SIM_SUBSAMPLE` keeps a fraction of the analysis individuals *within each
population*. The base analysis set is 2 individuals per patch × 80 patches = 160
(Nemo samples 80 patches × 4; the pipeline keeps every second individual).

| fraction | individuals | per patch |
|---|---|---|
| `1` | 160 | 2 |
| `0.75` | 120 | half the patches keep 2, half keep 1 |
| `0.5` | 80 | 1 |

0.75 does not divide evenly into 2-per-patch, so the fractional part is spread
across patches rather than rounded away: every patch keeps `floor()`, then a random
subset gets one more, giving exactly `round(0.75 × 160)` individuals overall. The
draw is seeded from the file name, so a given file always yields the same subsample.

Note that at 0.5 every patch has exactly one individual — no within-patch
replication left, which changes what the GRM can see. The MAF filter runs after
subsampling, so smaller samples also keep fewer markers; that is intentional (it is
what a smaller study would actually see), but it means marker sets differ between
fractions.

Each bundle records the fraction it was built with as `subsample`.

## Decay windows: `n_win_decay = 10`

Raised from 5, which was chosen for far smaller SNP sets. At this density 5 left only
~9 windows per chromosome, and `summarize_decay()` needs **≥ 5 windows with a valid
fit** (`!is.na(a) & a > 0 & !is.na(c)`). Chromosome 5 produced 8 windows of which only
4–5 fitted, so it straddled that gate and failed stochastically — 91 of 92 failures in
the first full batch were chr5. At 10 windows it gets ~19, well clear.

This **changes the fitted a/b/c** (see `?compute_LD_decay`, "Comparing decay
estimates"), so 5-window and 10-window bundles are not comparable and must never be
mixed. The earlier set was moved to `archive_w5/` rather than deleted.
`popgen_sim_data/` was kept — those summaries are computed from raw genotypes and never
touch `LD_decay`.

## Nulls (`run_nulls.sh`, `run_sim_nulls.R`)

The sim counterpart of `~/3sp_lfmm_perm/run_paired_nulls_3sp.R`: same `ld_null` bundle
format, same per-`(type, b)` seeding so EMMAX and LFMM see identical surrogate
phenotypes, same checkpoint-and-skip. One cell = 10 pooled chromosomes; one surrogate
phenotype goes through all ten, giving one pooled surrogate C.

| null | construction | engine | role |
|---|---|---|---|
| `genetic` | MVN(0, cell's mean GRM) | **EMMAX only** | does EMMAX invent signal from the kinship it corrects? |
| `latent` | MVN in the top-5 pooled genotype-PC subspace | **LFMM only** | same question for LFMM's latent factors |
| `global_perm` | env shuffled among the 80 patches | both | model-agnostic |
| `env_orth` | 48×48 env field shifted toroidally + rotated/reflected | both | **the arbiter** |
| `spatial` | MVN(0, Gaussian kernel over patch coords) | both | continuity with earlier runs |

The two home-field nulls are **not crossed** between engines. A null drawn from one
engine's own model of structure only answers that engine's specificity question;
scoring EMMAX against LFMM's latent basis (or the reverse) tests neither. The runner
skips those pairs, so each engine gets its own home-field null plus the
method-agnostic ones — which are what the cross-engine comparison rests on.

`env_orth` is the sim-specific one, standing in for 3sp's `region_perm` (the lattice has
no localities). A torus shift is a *relabelling of positions*, so the empirical
autocorrelation function is preserved exactly while the field moves off the true one;
draws are screened to |cor| < 0.3 against the observed env and then residualised, so the
surrogate is exactly orthogonal to it — same spatial structure, none of the signal.
Every null is residualised against the observed env, so none can carry the real signal.

**Measured costs** (`emmax_fast` 0.10 s/scan, `ld_cscore` 0.83 s — the C-score dominates;
LFMM 12.5 s per file):

| work | per cell × type | 30 cells | on 12 cores |
|---|---|---|---|
| EMMAX B=100, 4 types, one stage | 15.6 min | 31 core-h | 2.6 h |
| EMMAX, all three stages | | 94 core-h | 8 h |
| LFMM B=100, 2 types, full sample only | 3.5 h | 210 core-h | 17.5 h |

Dropping the cross pairs takes LFMM from 44 h to 17.5 h, so **B=100 is affordable for
both engines** — no need to run LFMM at a reduced B. Both engines then use the same
per-`(type, b)` seeds, so `global_perm`, `env_orth` and `spatial` surrogates are
identical between them and the comparison is exactly paired.

```bash
./run_nulls.sh 1 12 ; ./run_nulls.sh 0.75 12 ; ./run_nulls.sh 0.5 12   # EMMAX side
SIM_NULL_ENGINES=lfmm ./run_nulls.sh 1 12                             # LFMM, full sample
```

Cells run one at a time with draws forked inside — a loaded cell is 2–4 GB and forked
workers share it copy-on-write, so memory stays at one cell regardless of `nproc`.

## Outputs

| folder | contents |
|---|---|
| `regen_sim_data_<tag>/` | full-sample bundles |
| `regen_sim_data_<tag>_sub75/`, `_sub50/` | subsampled bundles |
| `popgen_sim_data/` | per-file `.rds` (summary, per-class diversity, per-patch table, BGS windows, ini settings) + `_summary.csv` |
| `nulls_nobgs<suffix>/` | `null_<engine>_<type>_V<V>_c<c>_env<e>.rds` (`ld_null` bundles) |
| `archive_w5/` | the superseded 5-window bundles |
| `logs_<tag><suffix>/`, `logs_nulls<suffix>/` | one log per cell |

The popgen summaries are computed from all 320 sampled individuals, *before* any
analysis subsetting, so they describe the simulation and are identical across
fractions. `run_batch.sh` therefore computes them only on the full-sample stage and
sets `SIM_POPGEN=""` for the subsampled ones.

## A caveat on the background-selection columns

`bgs_cor_pi_rec`, `bgs_cor_pi_delet` and friends are **within-run diagnostics, not
BGS measurements**. On a nobgs run — no deleterious selection at all —
`bgs_cor_pi_delet` still came out at +0.33 (+0.46 on the neutral chromosome),
because deleterious-site placement, recombination and SNP density are confounded in
the map. The window grid is fixed and shared across runs precisely so the paired
contrast `B_obs = pi_bgs / pi_nobgs` per window can be collated afterwards; that
form differences the confounding away and is the estimator to trust.

Related: the Chr1-vs-Chr2 contrast is *not* a BGS contrast. Both chromosomes carry
100 deleterious loci; what differs is the QTNs (100 vs 1), so an elevated Chr1 F_ST
is the local-adaptation footprint — and it shows up in nobgs runs too (0.055 vs
0.039).
