#!/usr/bin/env python3
"""Write the 'Files and variables' section of README_dryad.md for the staged
Dryad folder: a file table (from MANIFEST.tsv) in which files that share a
purpose are collapsed into one row (path pattern, count, size range), and one
column dictionary for the main tables. MANIFEST.tsv itself still lists every
file individually.

Usage: python3 tools/write_dryad_readme.py DRYAD_DIR README_PATH
If README_PATH exists, only its '### Files and variables' section (up to the next
'## ' heading) is replaced.
"""
import csv, re, sys, os

D, OUT = sys.argv[1], sys.argv[2]

## (regex, description, display pattern or None). The first matching rule wins.
## Files matching a rule with a display pattern are listed as ONE row under that
## pattern; files whose rule has display None get one row each.
RULES = [
    # --- simulation ------------------------------------------------------------
    (r"^simulation/nemo_inputs/", "Nemo 2.4.3 inputs from LDscnR-NEMO: `params_3spC_53cM/` (per-map neutral, deleterious and QTN loci maps, QTN allelic values, 10 environmental surfaces `env_*.txt`, 3 dispersal matrices `disp_mat_*.txt`, recombination maps `rds/rec_map*.rds`), `make_prod.sh` (exact regeneration of maps and the 1,540 `.ini` files, seeds embedded), the `R/` scripts and `ini/` templates it calls, and the Nemo large-population patch.", "simulation/nemo_inputs/**"),
    (r"^simulation/results/simulation_performance\.tsv$", "Pooled region-level precision and recall by BGS treatment, cell and method (Figure 2, Table S-performance).", None),
    (r"^simulation/results/simulation_performance_(c1|lfmm|region|diagnostic)\.tsv$", "The same scores for the 600 primary (c = 1) runs, for the LFMM arms, at region level for all methods, and the marker/unit-level diagnostics.", "simulation/results/simulation_performance_{c1,lfmm,region,diagnostic}.tsv"),
    (r"^simulation/results/simulation_(method|diagnostic)_contrasts(_c1)?\.tsv$", "Paired bootstrap differences in precision and recall between methods (primary, c = 1 and diagnostic scores).", "simulation/results/simulation_{method,diagnostic}_contrasts*.tsv"),
    (r"^simulation/results/simulation_stage2_region_details\.(tsv|rds)$", "One row per reported Stage-2 region in every run: method, size floor, span, markers, units, truth status, member markers and units, matched QTN (column `linked_qtn`); TSV and R object.", "simulation/results/simulation_stage2_region_details.{tsv,rds}"),
    (r"^simulation/results/simulation_qc\.tsv$", "QC record of the simulation pipeline (parser gates, grid completeness, marker-order checks).", None),
    (r"^simulation/results/simulation_summary_full\.rds$", "All pooled simulation summaries in one R object.", None),
    (r"^simulation/results/floor_stability_", "Floor-stability analysis of simulated regions: labels at half/double floors under three marker-density schemes, contrasts, region audits and summaries.", "simulation/results/floor_stability_*"),
    (r"^simulation/results/simulation_floor_", "Size-floor sweep and its decomposition into unit exclusion vs aggregation (Figures S-floor-decomposition).", "simulation/results/simulation_floor_*"),
    (r"^simulation/results/figure_simulation_.*_data\.tsv$", "Data plotted in the named figures.", "simulation/results/figure_simulation_*_data.tsv"),
    (r"^simulation/results/simulation_truth_sensitivity", "Precision and recall with the detectable-QTN threshold varied (Figure S-truth-sensitivity).", "simulation/results/simulation_truth_sensitivity*.tsv"),
    (r"^simulation/results/simulation_null_truth", "Null calibration of reported regions: expected null discoveries vs observed and true false-discovery proportion.", "simulation/results/simulation_null_truth*.tsv"),
    (r"^simulation/results/simulation_stage2_(size|opportunity)", "Region size vs false-positive status, and the opportunity null.", "simulation/results/simulation_stage2_{size,opportunity}_*.tsv"),
    (r"^simulation/results/simulation_bgs_", "Background-selection checks: marker retention, heterozygosity, BGS x aggregation interaction (Figure S-BGS-QC).", "simulation/results/simulation_bgs_*.tsv"),
    (r"^simulation/results/simulation_bootstrap_", "Sensitivity of the pooled estimates to the resampling scheme.", "simulation/results/simulation_bootstrap_*"),
    (r"^simulation/results/simulation_chrsplit_", "Analyses with chromosomes split into separate units (Table S-chromosome-split).", "simulation/results/simulation_chrsplit_*_c1.tsv"),
    (r"^simulation/results/simulation_marker_", "Marker-level redundancy: significant markers tagging the same QTN; singleton vs multi-marker units.", "simulation/results/simulation_marker_*.tsv"),
    (r"^simulation/results/simulation_qtn_detection", "Detection of individual QTN by effect size and method (Figure S-QTN-detection).", "simulation/results/simulation_qtn_detection*_c1.tsv"),
    (r"^simulation/results/simulation_env_structure_alignment", "Environment-structure alignment of the simulated environments and its relation to false discoveries.", "simulation/results/simulation_env_structure_alignment*.tsv"),
    (r"^simulation/results/simulation_(lfmm_portability|region_granularity)", "Floor effect on LFMM vs EMMAX; region- vs unit-level scoring.", "simulation/results/simulation_{lfmm_portability,region_granularity}_contrast.tsv"),
    (r"^simulation/results/exploratory_", "Exploratory LD/Fst vs null-burden analysis; not cited as a result.", "simulation/results/exploratory_*.tsv"),
    # --- empirical, both panels ------------------------------------------------
    (r"^empirical_(3sp|9sp)/bundle_slim\.rds$", None, None),
    (r"/_receipt\.rds$", "Stage receipts: inputs (paths, SHA-256), parameters, LDscnR version and run time of the stage that wrote each folder.", "empirical_*/*/_receipt.rds"),
    (r"/01_inputs/inputs\.csv$", "SHA-256 of every input file each module read (paths are those of the original machine).", "empirical_*/01_inputs/inputs.csv"),
    (r"^empirical_3sp/03_EMMAX/scan\.rds$", "Three-spined EMMAX consensus and Simes results: per arm `test` (ld_outlier_test(): `units`, `regions`), `null` (ld_outlier_perm(): observed, surrogates, p; 1,000 / 200 locality-preserving permutations) and `rotation` (EcoPeak relocation test: observed, null_mean, fold, p); plus the EcoPeak intervals and chromosome lengths.", None),
    (r"^empirical_3sp/04_lfmm/lfmm\.rds$", "LFMM Simes arm: `test`, `rotation` and the marker-wise LFMM p-values (`lfmm_p`).", None),
    (r"^empirical_3sp/05_manuscript/", "Generated LaTeX macros (`values_panel_module3sp.tex`) and main-text Figure 4.", "empirical_3sp/05_manuscript/{values_panel_module3sp.tex,figure4_stickleback_validation.pdf}"),
    (r"^empirical_3sp/09_ld_recombination/ld_recombination\.rds$", "LD decay vs pedigree recombination: window-level decay rate, ld_w and map rate, correlations and ROC (Supplementary Figs. S1-S2).", None),
    (r"^empirical_3sp/10_sensitivity/sensitivity\.rds$", "Full results for every point of the rho / size-floor grid (test, permutation null and relocation test, three analyses).", None),
    (r"^empirical_3sp/10_sensitivity/sensitivity_summary\.csv$", "One row per grid point and analysis (Table 5).", None),
    (r"^empirical_3sp/13_goldilocks_units/goldilocks_units\.rds$", "Unit representations compared at matched test counts (Table 6).", None),
    (r"^empirical_3sp/13b_goldilocks_regions/", "The same comparison counted in outlier regions (Table S7): R object, region list, summary and the LaTeX table as written.", "empirical_3sp/13b_goldilocks_regions/*"),
    (r"^empirical_3sp/14_floor_tiers/(ecopeak_by_tier|ecopeak_single_floor)\.tsv$", "EcoPeak relocation test by floor-stability label (Table S-tier) and for regions found at a single comparison floor.", "empirical_3sp/14_floor_tiers/ecopeak_{by_tier,single_floor}.tsv"),
    (r"^empirical_(3sp/14|9sp/04)_floor_tiers/regions_by_tier\.tsv$", "Regions reported at the primary floor with their floor-stability label.", "empirical_{3sp/14,9sp/04}_floor_tiers/regions_by_tier.tsv"),
    (r"^empirical_(3sp/14|9sp/04)_floor_tiers/(regions_all_floors|tier_counts)\.tsv$", "Regions found at each of the three floors; label counts.", "empirical_{3sp/14,9sp/04}_floor_tiers/{regions_all_floors,tier_counts}.tsv"),
    (r"^empirical_3sp/15_marker_postfilter/marker_postfilter_ecopeak\.tsv$", "Size floor applied as a post-filter to marker-wise scans (Table S-postfilter).", None),
    (r"^empirical_3sp/16_ldw_peaks/", "Local-LD (ld_w) peaks, their overlap with reported regions and EcoPeaks, per-unit ld_w, and regions an ld_w filter would affect.", "empirical_3sp/16_ldw_peaks/*.tsv"),
    (r"^empirical_3sp/17_span_vs_simulation/span_comparison\.tsv$", "Three-spined region spans vs simulated true-positive spans.", None),
    (r"^empirical_3sp/goldilocks_cache/stage2_Chr\d+_d1e\+05\.rds$", "Table 6 input: Stage-1 and Stage-2 partitions of each chromosome (100 kb distance threshold), per-individual consensus matrices removed.", "empirical_3sp/goldilocks_cache/stage2_Chr{1..18,20,21}_d1e+05.rds"),
    (r"^empirical_3sp/goldilocks_cache/", "Table 6 inputs from the LD-pruning-paper analysis: marker positions and EMMAX p-values, kinship and ecotype, and the proxy truth loci.", "empirical_3sp/goldilocks_cache/{map_only.rds,K_eco.rds,3sp_cmin_loci.csv}"),
    (r"^empirical_9sp/02b_env_structure_alignment/alignment_observed\.tsv$", "Observed environment-structure alignment per panel and phenotype.", None),
    (r"^empirical_9sp/02b_env_structure_alignment/", "Alignment null draws (permutation and GRM surrogates), their summary, and all results as an R object.", "empirical_9sp/02b_env_structure_alignment/{alignment.rds,alignment_null_draws.tsv,alignment_null_summary.tsv}"),
    (r"^empirical_9sp/03_EMMAX/scan\.rds$", "Nine-spined EMMAX consensus and Simes results (ecotype adjusted for lineage; within-lineage permutation), structured as the three-spined `scan.rds` without relocation tests.", None),
    (r"^empirical_9sp/03[bc]_", "Nine-spined sensitivity analyses: EL and Admixed lineages only (03b); GRM-based surrogate phenotypes (03c).", "empirical_9sp/03{b,c}_*/scan.rds"),
    (r"^empirical_9sp/03d_diagnose_variance_components/diagnostics\.rds$", "Variance-component diagnostics for both panels: lambda_GC, REML h2, random-phenotype h2, alternative GRM basis.", None),
]
BUNDLE = ("Stage-1 bundle without genotypes, one per panel: `n_markers`, `n_individuals`, `map` (marker, Chr, Pos, maf, ld_w_095), "
          "`eco` (0 freshwater, 1 marine), `pheno` (index, ecotype, population, lineage, anonymous locality_group), `LD_decay`, "
          "`stage1` (map_snp, clusters, pruned), `GRM`, `grm_markers`, `settings`; three-spined also `fig1_window` "
          "(r^2 of the 26 chromosome-9 markers in Figure 1B-C).")

## One dictionary: columns shared by several tables are listed once.
DICT = [
    ("`tag`, `cell`, `method`", "simulation", "BGS treatment (nobgs, bgs); selection/dispersal cell V<QTN variance>_c<dispersal>; analysis (e.g. emmax_snp_region = marker-wise, emmax_simes_region, emmax_consensus_region)"),
    ("`TP`, `FP`, `FN`", "simulation", "reported regions matched / not matched to a detectable QTN; detectable QTN not recovered"),
    ("`n_detectable_qtn`", "simulation", "QTN with MAF > 0.1 explaining >= 5% of its chromosome's additive variance"),
    ("`coverage`", "simulation_performance", "fraction of detectable QTN covered by an eligible Stage-1 unit"),
    ("`*_ci_lo`, `*_ci_hi`", "simulation", "95% cluster-bootstrap interval (2,000 resamples of the ten map/burn-in replicates)"),
    ("`diff_precision`, `diff_recall`", "*_contrasts", "method minus `reference_method`"),
    ("`grid`", "sensitivity_summary", "Stage-1 rho and size floor, e.g. rho=0.50_sf=8"),
    ("`tested`, `significant`", "sensitivity_summary", "units meeting the floor; units significant after BH correction"),
    ("`on_ecopeak`, `ecopeak`, `frac`", "empirical", "reported regions overlapping an EcoPeak, and their fraction"),
    ("`null_mean`, `fold` / `rotation_fold`", "empirical", "mean EcoPeak overlaps over 10,000 within-chromosome relocations; observed / null mean"),
    ("`p` / `rotation_p`", "empirical", "relocation p, counting the observed placement as one relocation"),
    ("`perm_p`, `realised_fdr`", "sensitivity_summary", "permutation p (EMMAX only); null/obs ratio = mean significant units under permutation / observed (not an FDR)"),
    ("`occupancy`", "regions_by_tier", "fraction of markers in the region's span that belong to its significant units"),
    ("`floor_lo`, `floor`, `floor_hi`, `found_lo`, `found_hi`, `tier`", "regions_by_tier", "comparison and primary size floors; whether the region is recovered at each comparison floor; floor-stability label (number of floors)"),
    ("`set`, `statistic`", "ecopeak_by_tier", "floor-stability label (all, 1, 2, 3, >=2); unit (consensus) or simes"),
    ("`engine`, `floor`, `snps_kept`", "marker_postfilter_ecopeak", "marker-wise scan (LFMM, EMMAX); minimum Stage-1 cluster size (1 = no filter); significant markers retained"),
    ("`structure_group`", "alignment_observed", "null-defining group (pop_locality or lineage)"),
    ("`grm_axes_r2_population`, `structure_r2_*`", "alignment_observed", "R^2 of the (population-mean) phenotype on the leading GRM axes / on the structure groups"),
]

def rule_for(path):
    for pat, d, disp in RULES:
        if re.search(pat, path): return pat, (BUNDLE if d is None else d), disp
    return None, "", None

def size(b):
    b = float(b)
    for u in ("B", "KB", "MB", "GB"):
        if b < 1000 or u == "GB": return f"{b:.0f} {u}" if u == "B" else f"{b:.1f} {u}"
        b /= 1024

rows = list(csv.DictReader(open(os.path.join(D, "MANIFEST.tsv")), delimiter="\t"))
total = sum(int(r["bytes"]) for r in rows)
groups, order = {}, []
missing = []
for r in sorted(rows, key=lambda r: r["file"]):
    pat, desc, disp = rule_for(r["file"])
    if pat is None: missing.append(r["file"])
    key = disp if disp else r["file"]
    if pat is not None and RULES[[p for p, _, _ in RULES].index(pat)][1] is None and disp is None:
        key = "empirical_{3sp,9sp}/bundle_slim.rds"   # the two slim bundles share one row
    if key not in groups: groups[key] = {"desc": desc, "bytes": []}; order.append(key)
    groups[key]["bytes"].append(int(r["bytes"]))

L = ["### Files and variables", "",
     f"The dataset has {len(rows) + 1} files ({total/1e9:.2f} GB) in three folders that mirror the analysis modules of the "
     "archived code. Set `LDSCNR_DATA` to the unpacked folder and the manuscript generators and figure scripts read "
     "everything from it.", "",
     "No genotype matrices are included; the only individual-level matrices are the kinship matrices (`GRM` in each "
     "`bundle_slim.rds`, and `K_eco.rds`). The Stage-1 bundles are slim copies without genotypes, and the Table 6 cache "
     "lacks the per-individual consensus matrices, so Table 6 is included as output but recomputing it needs the "
     "genotypes. Sample metadata is limited to ecotype, population and lineage labels and an anonymous locality code; "
     "the full sample table is in the cited Zenodo records.", "",
     "#### Column dictionary", "",
     "Columns of the main tables whose meaning is not evident from the name (`n_*`, `Chr`, `from`, `to`, `precision`, "
     "`recall`, `regions` are counts, coordinates or rates as named).", "",
     "| column | in | meaning |", "|---|---|---|"]
L += [f"| {c} | {w} | {m} |" for c, w, m in DICT]
L += ["", "#### Files", "",
      "Rows with a pattern (`*`, `{a,b}`, `{1..18}`) stand for several files; `MANIFEST.tsv` lists every file "
      "individually with its size and SHA-256.", "",
      "| file(s) | n | size | description |", "|---|---:|---:|---|"]
for k in order:
    b = groups[k]["bytes"]
    sz = size(b[0]) if len(b) == 1 else f"{size(min(b))}–{size(max(b))} ({size(sum(b))})"
    L.append(f"| `{k}` | {len(b)} | {sz} | {groups[k]['desc']} |")
L.append("| `MANIFEST.tsv` | 1 | | Every file with its size (bytes) and SHA-256. |")
section = "\n".join(L) + "\n"

if os.path.exists(OUT):
    txt = open(OUT).read()
    m = re.search(r"^### Files and variables.*?(?=^## |\Z)", txt, re.S | re.M)
    txt = txt[:m.start()] + section + "\n" + txt[m.end():] if m else txt.rstrip() + "\n\n" + section
else:
    txt = "# LDscnR manuscript: analysis outputs\n\n" + section
open(OUT, "w").write(txt)
print(f"{len(rows)} files in {len(order)} rows, {total/1e9:.2f} GB; {len(missing)} without description: {missing[:10]}")
