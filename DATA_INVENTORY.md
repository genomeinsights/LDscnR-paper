# Untracked data inventory — LDscnR-paper

Generated 2026-08-25. **Nothing here has been deleted or moved.** This is the
report requested before any decision about the ~35 GB of untracked data.

The repo is 35 GB; `.git` is 11 MB; 310 files are tracked. Everything below is
outside git, so **deleting any of it is permanent** — a branch does not protect it.

---

## 1. Gitignored SOURCE — at risk, should be ADDED to git, not deleted

These are scripts, not data. They are currently protected by nothing.

| Path | Size | What | Why it matters |
|---|---|---|---|
| `NEMO/` (`R/`, `cluster/`, `ini/`, `run/`, `README.md`, `PATCHES.md`) | 36 KB + 2.7 MB source | The BGS re-implementation module written for a collaborator to run on the cluster | 741 source files, gitignored by the `NEMO/` rule. Not backed up, not shareable through the repo, and includes the macOS Nemo loader patch notes |
| `3sp_data/run_3sp_perm_global.R` | 4 KB | Generator of the 3sp global-permutation EMMAX output | The only record of how `out_emx_global.rds` was made |
| `Nemo3/*.R` | 2 files | Nemo map parsing | Referenced from `legacy/root/parse_nemo_map.R` |
| `R_080726/*.R` | 3 files | Poster-figure scripts | Referenced from `module_sticklebacks/14_poster_manhattan.R` |

**Recommendation:** narrow the `.gitignore` rules from whole directories to data
extensions, so `NEMO/**/*.R` and friends are tracked while `NEMO/params_V4/`
stays out. This is the single highest-value change in this report.

## 2. Live inputs — needed by the CURRENT modules

| Path | Size | Used by | Regenerable? |
|---|---|---|---|
| `parsed_sim_data/` | 245 MB | `module_sim_LDscnR/regen_sim_data.R` | Only by re-parsing `Nemo_sim/` |
| `NEMO/params_V4/`, `NEMO/grid/` | 248 MB | the NEMO module | Re-runnable on the cluster, at cost |
| `gds_3sp.gds` | 27 MB | `module_sticklebacks_LDscnR` | From the joint VCF |

## 3. Expensive upstream — keep, but belongs outside the repo

| Path | Size | Notes |
|---|---|---|
| `Nemo_sim/` | 14 GB | Raw Nemo simulation output; parsed into `parsed_sim_data/`. Costs cluster time to regenerate. The single largest item |
| `3sp_data/` | 3.6 GB | Five ~569 MB EMMAX permutation outputs (`out_emx_*.rds`, `emmax_perm_reginal.rds`) plus small derived bundles. Derived, but expensive |
| `emmax_perm_reginal.rds` (root) | 596 MB | Duplicate of the copy inside `3sp_data/` — verify before touching |

## 4. Superseded — referenced only from `legacy/` or the old modules

Safe to retire *if* the old modules are retired with them. Check each reference first.

| Path | Size | Only referenced from |
|---|---|---|
| `Nemo3/` | 969 MB | `legacy/root/parse_nemo_map.R` |
| `LFMM_3sp/` | 159 MB | `module_sticklebacks/` (superseded by `module_sticklebacks_LDscnR`) |
| `parsed_sim_data_genomes/` | 57 MB | `legacy/root/R_sim.R`, `legacy/R_sim/score_ORs_sim.R` |
| `OR_performance/`, `OR_performance2/` | 42 MB | `legacy/root/eMLG_complexity_reduction.R` |
| `R_080726/` | 56 KB | `module_sticklebacks/14_poster_manhattan.R` (rescue the 3 scripts first) |

## 5. Session junk — no script reads these

Verified by grep: nothing in the repo loads the root `.RData` or `.RDataTmp`.

| Path | Size |
|---|---|
| `.RData` | 7.1 GB |
| `.RDataTmp` | 5.9 GB |
| `EL_tmp/` | empty |
| `claude_snapshot_summaries/` | 24 KB |

**13 GB, zero references.** The lowest-risk reclaim available, and the only
category I would call unambiguous.

---

## Reproducibility gaps found while surveying

Independent of size, three result files have **no committed generator**:

1. `results/null_emmax_{genetic,global_perm,region_perm,spatial,latent}_3sp.rds`
2. `results/null_lfmm_region_perm_3sp.rds`
3. the `region_emp_pvals_emmax_*` CSVs derived from them

`lfmm_permutation_null_3sp.R` builds a *global* permutation
(`permute_among_pops()`, no `by = loc`) and writes `null_lfmm_perm_3sp.rds` —
not the `lfmm_region_perm` file on disk. Since `results/` is gitignored and the
scripts are the record, these cannot currently be rebuilt from a fresh clone.
Writing those generators is a prerequisite for any deletion pass, because a
result nobody can rebuild is a result nobody dares delete.
