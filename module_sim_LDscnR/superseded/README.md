# Superseded

`analyse_one_dataset.R` — the C-score (`ld_scan()`/Fang et al. consistency score) driver
for the simulation benchmark. Archived 2026-09-03: PK confirmed the sims and the 3sp panel
are on the same method — stage-1 clustering as the test unit, consensus dosage or Simes per
unit, BH across units, rotation/permutation null — and this script is not it.

It calls `ld_scan()` bare, which is now unexported from LDscnR (still reachable as
`LDscnR:::ld_scan()` if ever needed — see LDscnR commit `ddb7b65`).

Not verified against fe's current live driver before archiving; `grm_comparison.R` was
checked and has no C-score references, which is the evidence this one is stale rather than
confirmation from the session that owns this module.
