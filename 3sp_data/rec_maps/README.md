# Threespine stickleback recombination map — three F2 crosses, gasAcu1 coordinates

Every rate here comes from **crossovers counted in pedigrees**. No population LD enters
the estimation, so the map is safe to use alongside LD-based analyses.

## Sources

| cross | parents | F2 | markers | map length | published rate |
|---|---|---|---|---|---|
| FTC | Fishtrap Ck (WA, freshwater) × Little Campbell R. (BC, marine) | 357 | 1,001 | 1,570 cM | — |
| BEPA | Bear Paw Lake (AK, freshwater) × Little Campbell R. (BC, marine) | 360 | 978 | 1,963 cM | — |
| ROESTI | Lake Constance lake × stream | 282 | 1,872 | 1,251 cM | 3.11 cM/Mb |

FTC and BEPA: Glazer, Killingbeck, Mitros, Rokhsar & Miller 2015, G3 5:1463–1472,
doi:[10.1534/g3.115.017905](https://doi.org/10.1534/g3.115.017905) (CC-BY), Files S1–S3.
ROESTI: Roesti, Moser & Berner 2013, Mol Ecol 22:3014–3027,
doi:[10.1111/mec.12322](https://doi.org/10.1111/mec.12322), Appendix S4.

`build_glazer_map.py` and `merge_maps.R` only reshape the published tables — nothing is
re-estimated. As a check on that: recomputing the Roesti map from Appendix S4 returns
**3.11 cM/Mb** genome-wide, matching the published figure exactly.

## Files

- **`stickleback_recomb_3crosses_gasAcu1.tsv`** — the deliverable. 880 bins tiling all 21
  gasAcu1 chromosomes at a median 489 kb.
- `glazer2015_rate_bins.tsv`, `glazer2015_markers.tsv` — the two-cross Glazer layer alone.
- `source/roesti_AppendixS4.txt` — Roesti Appendix S4 verbatim, as retrieved.

### Columns

`old_chr`/`old_start`/`old_end` are **gasAcu1** (BROAD S1, 2006); `new_*` are the Glazer
revised assembly (gasAcu1-4). Bins come from Glazer File S3, which carries both.

| column | meaning |
|---|---|
| `rate_FTC`, `rate_BEPA`, `rate_ROESTI` | local rate, cM/Mb, per cross |
| `rel_*` | same, divided by that cross's own genome-wide mean |
| `n_crosses` | how many crosses scored this bin (1–3) |
| `rel_consensus` | mean of the available `rel_*` |
| `rel_sd` | raw spread of `rel_*` across crosses. Descriptive only — **do not filter on it**, see below |
| `rank_sd` | spread of each cross's genome-wide *percentile* for this bin |
| `disagree` | `rank_sd` ranked **within rate decile**; uniform on [0,1] and orthogonal to rate. 0 = the crosses agree about this bin relative to others of similar rate |
| `concordant` | ≥2 crosses **and** `disagree` ≤ 0.75 (655 of 880 bins). ⚠️ **Redefined in `d74fdc7`** — previously thresholded `rel_sd`; **206 of 880 bins flipped**, so any conclusion recorded against the earlier flag is stale and needs recomputing. [Why](#why-concordant-is-not-built-on-rel_sd) |
| `rate_consensus_cMperMb` | `rel_consensus × 3.11`, i.e. the consensus landscape on Roesti's published absolute scale |

Rates are the slope of cM against Mb over each marker's two flanking markers, taken in
the **revised** assembly order (the order these studies corrected). Roesti's marker-level
rates are carried back to their gasAcu1 positions and averaged into the Glazer bins by
overlap length.

## Re-fetching the source files

`source/` is gitignored, so **a fresh clone will not have it** and neither build script
will run. The four files are all published supplementary material and are recoverable by
hand, but read this first.

> **The obvious way to fetch them fails silently.** Wiley — and to a lesser extent the PMC
> file endpoint — sits behind Cloudflare, which answers a scripted `curl`/`wget` with
> **HTTP 200 and a ~6 KB HTML interstitial**, not an error. That lands on disk under the
> name you asked for and parses as *something*, so the failure presents as a corrupt or
> truncated download rather than as an access block, and a rebuild can half-succeed.
> `roesti_AppendixS4.txt` in particular must be retrieved through a browser.

| file | source |
|---|---|
| `supp_g3.115.017905_FileS{1,2,3}.xlsx` | Glazer et al. 2015, G3 5:1463–1472, doi:[10.1534/g3.115.017905](https://doi.org/10.1534/g3.115.017905), Files S1–S3 (CC-BY). Also mirrored in the Europe PMC supplementary package for PMC4502380, which does allow scripted download. |
| `roesti_AppendixS4.txt` | Roesti, Moser & Berner 2013, Mol Ecol 22:3014–3027, doi:[10.1111/mec.12322](https://doi.org/10.1111/mec.12322), Appendix S4. Wiley only — browser required. |

Both scripts verify size and hash before doing any work and **refuse to proceed** on a
file that is absent, truncated, or not what it claims, naming the Cloudflare case
explicitly when the size looks like an interstitial:

| file | bytes | |
|---|---|---|
| `roesti_AppendixS4.txt` | 120,011 | md5 `255a8f02adc2bed9993c5284edfe3e0c` |
| `supp_g3.115.017905_FileS1.xlsx` | 1,342,078 | sha256 `f8bd177572bf8c2b…` |
| `supp_g3.115.017905_FileS2.xlsx` | 1,857,861 | sha256 `c695fca8dbe02dda…` |
| `supp_g3.115.017905_FileS3.xlsx` | 244,412 | sha256 `d3cc2fb083fd1efc…` |

If a hash mismatches at full size, the file was revised upstream or edited locally — do
not build on it without deciding which you have.

## Agreement between crosses

Local rate, across bins scored by both:

| | BEPA | ROESTI |
|---|---|---|
| **FTC** | r = 0.861 (n = 820) | r = 0.687 (n = 836) |
| **BEPA** | — | r = 0.705 (n = 824) |

The two Glazer crosses agree more closely with each other than either does with Roesti,
which is what you would expect: they share a marker set, a mapping pipeline, and the LITC
marine parent, whereas Roesti is an independent lake–stream cross scored on its own
markers. Treat `r ≈ 0.7` as the honest reproducibility of this landscape across
populations and labs, and `disagree` as the per-bin version of that.

### Why `concordant` is not built on `rel_sd`

`rel_sd` is a standard deviation of untransformed rates, so it scales with the mean:
`cor(rel_sd, rel_consensus) = +0.635`. Thresholding it does not select well-measured
bins, it selects **low-recombination** bins — under that rule 100% of bottom-decile-rate
bins passed against 25% of top-decile, and mean rate among passing bins was 2.03 cM/Mb
against 6.45 among failing ones. Filtering a reference set on it truncates the top of the
rate range, so any correlation computed against recombination attenuates from range
restriction alone, whatever the measurement noise is doing.

`rank_sd` removes the magnitude scaling but is still not orthogonal
(`cor = −0.222`): percentiles are bounded, so bins at the very top or bottom of every
cross's ordering have no room to disagree, and thresholding it keeps the extremes and
drops the middle.

`disagree` ranks `rank_sd` *within* rate decile, which makes it uniform by construction
and carries no rate information: `cor(disagree, rel_consensus) = −0.018`, retention is
73–76% in every one of the ten deciles, and mean rate is 3.15 cM/Mb among concordant bins
against 3.21 among the rest. `merge_maps.R` prints all three correlations and the
by-decile retention on every run, so the check is visible rather than asserted.

> **Changed 2026-09-03.** `concordant` previously thresholded `rel_sd`. 206 of 880 bins
> (23%) change membership. Anything built on the earlier flag should be recomputed.

## Coverage against the 3sp SNP set

790,578 SNPs (gasAcu1, Chr1–Chr21): **99.6%** fall in a bin carrying a consensus rate;
**93.1%** in a bin scored by all three crosses.

```r
b <- fread("stickleback_recomb_3crosses_gasAcu1.tsv")
setkey(b, old_chr, old_start, old_end)
snps[, chr_i := as.integer(sub("Chr", "", Chr))]
hit <- foverlaps(snps[, .(old_chr = chr_i, old_start = Pos, old_end = Pos)], b, type = "within")
```

## Caveats

- **Sex-averaged.** F2 maps cannot separate maternal from paternal crossovers, and
  stickleback heterochiasmy is large (females ~1.6× males; Sardell et al. 2018).
- **~490 kb resolution.** The landscape, not hotspots.
- **Absolute cM/Mb differs threefold-ish between crosses** (bin means 3.92 / 4.91 / 3.22).
  F2 map length inflates with marker number and genotyping error, which is why the
  consensus is built from within-cross *relative* rates and only then put on Roesti's
  published scale. If absolute rate enters a calculation, say which scale you used.
- **Near-zero bins at some chromosome starts are real features of the input maps**, not
  artefacts of the merge: runs of co-segregating markers get epsilon cM offsets, so the
  local slope is ~0 there.
- 98 Roesti markers and 162 Glazer bins sit on scaffolds gasAcu1 left unanchored; they are
  dropped from this file because they have no gasAcu1 coordinate to look up.
- Chr19 (the sex chromosome) is present in the source maps but absent from the 3sp SNP set.
- The Glazer corrigendum adds only an omitted methods citation; it does not touch the data.
