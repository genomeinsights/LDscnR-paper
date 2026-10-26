#!/usr/bin/env python3
"""Write the 'Files and variables' section of README_dryad.md for the staged
Dryad folder: every file (from MANIFEST.tsv) with its size and a one-line
description, plus a column dictionary for the main tables.

Usage: python3 tools/write_dryad_readme.py DRYAD_DIR README_PATH
If README_PATH exists, only its '### Files and variables' section (up to the next
'## ' heading) is replaced.
"""
import csv, re, sys, os

D, OUT = sys.argv[1], sys.argv[2]

RULES = [
    # empirical, both panels
    (r"bundle_slim\.rds$", "Stage-1 bundle without genotypes: list(n_markers, n_individuals, map [marker, Chr, Pos, maf, ld_w_095], eco [0 freshwater / 1 marine], pheno [index, ecotype, population, lineage, locality_group (anonymous code)], LD_decay [per-chromosome decay fits, decay_sum, decay_model, ld_ws], stage1 [map_snp, clusters, pruned], GRM [GCTA kinship on Stage-1 representatives], grm_markers, settings, fig1_window [chromosome-9 window r^2 for Figure 1B-C; three-spined only])."),
    (r"_receipt\.rds$", "Stage receipt: inputs (paths, SHA-256), parameters, LDscnR version and run time of the stage that wrote this folder."),
    (r"01_inputs/inputs\.csv$", "SHA-256 of every input file read by the module (provenance; paths are those of the original machine)."),
    (r"empirical_3sp/03_EMMAX/scan\.rds$", "EMMAX consensus and Simes results: per arm, ld_outlier_test() (units: unit_id, Chr, from, to, n_markers, p, q, significant; regions: Chr, from, to, n_units, n_markers, occupancy), ld_outlier_perm() (observed, surrogates [1,000 consensus / 200 Simes], p) and ld_region_rotation() against the EcoPeaks (observed, null_mean, fold, p); plus the EcoPeak intervals and chromosome lengths used."),
    (r"empirical_3sp/04_lfmm/lfmm\.rds$", "LFMM Simes arm: ld_outlier_test() result, EcoPeak relocation test, and the marker-wise LFMM p-values (lfmm_p, named by marker)."),
    (r"05_manuscript/values_panel_module3sp\.tex$", "LaTeX macros with the three-spined values cited in the manuscript (generated)."),
    (r"05_manuscript/figure4_stickleback_validation\.pdf$", "Main-text Figure 4 (three-spined EMMAX and LFMM Manhattan plots with EcoPeaks)."),
    (r"09_ld_recombination/ld_recombination\.rds$", "LD decay vs pedigree recombination: window-level decay rate a, ld_w and recombination rate, per-chromosome correlations and ROC results (Supplementary Figs. S1-S2)."),
    (r"10_sensitivity/sensitivity\.rds$", "Full results of the rho / size-floor sensitivity grid (every grid point: test, permutation null and relocation test for consensus, Simes and LFMM)."),
    (r"10_sensitivity/sensitivity_summary\.csv$", "One row per grid point and analysis (Table 5); columns below."),
    (r"13_goldilocks_units/goldilocks_units\.rds$", "Unit-representation comparison at matched test counts (Table 6): per partition and representation, tested units, min p, Bonferroni/BH counts, median p and rank."),
    (r"13b_goldilocks_regions/goldilocks_regions\.rds$", "Goldilocks comparison counted in outlier regions (Table S7)."),
    (r"13b_goldilocks_regions/goldilocks_regions_list\.csv$", "Regions behind Table S7, one row per region and representation."),
    (r"13b_goldilocks_regions/goldilocks_regions_summary\.csv$", "Summary counts behind Table S7."),
    (r"13b_goldilocks_regions/table6_goldilocks_units\.tex$", "LaTeX table as written by the stage (copy)."),
    (r"14_floor_tiers/ecopeak_by_tier\.tsv$", "EcoPeak relocation test by floor-stability label (Table S-tier); columns below."),
    (r"14_floor_tiers/ecopeak_single_floor\.tsv$", "EcoPeak relocation test for regions found at only one comparison floor."),
    (r"04_floor_tiers/tier_counts\.tsv$", "Number of reported regions per floor-stability label."),
    (r"floor_tiers/regions_by_tier\.tsv$", "Reported regions at the primary floor with their floor-stability label; columns below."),
    (r"floor_tiers/regions_all_floors\.tsv$", "Regions reported at each of the three floors (primary, half, double)."),
    (r"15_marker_postfilter/marker_postfilter_ecopeak\.tsv$", "Size floor applied as a post-filter to marker-wise scans (Table S-postfilter); columns below."),
    (r"16_ldw_peaks/ldw_peaks\.tsv$", "Peaks of local LD support (ld_w) along the genome."),
    (r"16_ldw_peaks/ldw_peaks_in_regions\.tsv$", "ld_w peaks overlapping reported regions and EcoPeaks."),
    (r"16_ldw_peaks/ldw_units\.tsv$", "Stage-1 units with their ld_w summary."),
    (r"16_ldw_peaks/ldw_low_regions\.tsv$", "Reported regions with low ld_w."),
    (r"16_ldw_peaks/ldw_filter_missed_regions\.tsv$", "Regions an ld_w filter would have removed."),
    (r"17_span_vs_simulation/span_comparison\.tsv$", "Three-spined region spans vs simulated true-positive region spans (bp and relative to d_0.95)."),
    (r"goldilocks_cache/3sp_cmin_loci\.csv$", "Proxy truth set (EcoPeak-derived loci) for the Table 6 Bonferroni cross-check (from the LD-pruning-paper analysis)."),
    (r"goldilocks_cache/K_eco\.rds$", "Kinship matrix K (117 x 117) and ecotype vector used by Table 6 (from the LD-pruning-paper analysis)."),
    (r"goldilocks_cache/map_only\.rds$", "Per-marker positions and EMMAX p-values used by Table 6 (from the LD-pruning-paper analysis)."),
    (r"goldilocks_cache/stage2_Chr\d+_d1e\+05\.rds$", "Stage-1 and Stage-2 partitions of one chromosome (distance threshold 100 kb) used by Table 6; per-individual eMLG consensus matrices removed."),
    (r"02b_env_structure_alignment/alignment\.rds$", "Environment-structure alignment diagnostic, both panels (observed and null draws)."),
    (r"02b_env_structure_alignment/alignment_observed\.tsv$", "Observed environment-structure alignment per panel and phenotype; columns below."),
    (r"02b_env_structure_alignment/alignment_null_draws\.tsv$", "Alignment statistics for each null phenotype draw (permutation and GRM surrogates)."),
    (r"02b_env_structure_alignment/alignment_null_summary\.tsv$", "Summary of the alignment null draws."),
    (r"empirical_9sp/03_EMMAX/scan\.rds$", "Nine-spined EMMAX consensus and Simes results (ecotype residualised on lineage; within-lineage permutation null), structure as the three-spined scan.rds without relocation tests."),
    (r"03b_EMMAX_sensitivity_ELAdmixed/scan\.rds$", "Nine-spined sensitivity analysis restricted to the EL and Admixed lineages."),
    (r"03c_EMMAX_structurednull/scan\.rds$", "Nine-spined analysis against GRM-based (multivariate normal) surrogate phenotypes."),
    (r"03d_diagnose_variance_components/diagnostics\.rds$", "Variance-component diagnostics: genomic inflation (lambda_GC), REML h2, random-phenotype h2, for both panels and GRM bases."),
    # simulation inputs
    (r"nemo_inputs/params_3spC_53cM/map_ntrl\d+\.txt$", "Nemo map: neutral loci positions (cM) for one simulated chromosome/map."),
    (r"nemo_inputs/params_3spC_53cM/map_delet\d+\.txt$", "Nemo map: deleterious loci positions for one map (BGS treatment)."),
    (r"nemo_inputs/params_3spC_53cM/map_QTN\d+\.txt$", "Nemo map: QTN positions for one map."),
    (r"nemo_inputs/params_3spC_53cM/allelic_values_\d+\.txt$", "Nemo QTN allelic effects for one map."),
    (r"nemo_inputs/params_3spC_53cM/env_\d+\.txt$", "Environmental surface (one value per deme) for one environment replicate."),
    (r"nemo_inputs/params_3spC_53cM/disp_mat_.*\.txt$", "Dispersal matrix between demes for one migration setting (c)."),
    (r"nemo_inputs/params_3spC_53cM/rds/rec_map\d+\.rds$", "Recombination map for one simulated chromosome (R object read by the pipeline)."),
    (r"nemo_inputs/params_3spC_53cM/", "Nemo parameter-set file (LDscnR-NEMO params_3spC_53cM)."),
    (r"nemo_inputs/make_prod\.sh$", "Exact invocation that regenerates the maps and the 1,540 production .ini files (seeds embedded)."),
    (r"nemo_inputs/R/.*\.R$", "LDscnR-NEMO script used by make_prod.sh (map building, .ini generation, BGS checks)."),
    (r"nemo_inputs/ini/.*\.tmpl$", "Nemo .ini template (burn-in / adaptation)."),
    (r"nemo_inputs/PATCHES\.md$", "Why the Nemo 2.4.3 large-population patch is needed."),
    (r"nemo_inputs/nemo-macos-bigpop\.patch$", "Patch applied to Nemo 2.4.3 (nemo-release d3e321b)."),
    (r"nemo_inputs/README\.md$", "LDscnR-NEMO README (parameter choices for background selection)."),
    # simulation results
    (r"results/simulation_performance\.tsv$", "Pooled region-level precision and recall by BGS treatment, cell and method (Figure 2, Table S-performance); columns below."),
    (r"results/simulation_performance_c1\.tsv$", "As simulation_performance.tsv for the 600 primary (c = 1) simulations."),
    (r"results/simulation_performance_lfmm\.tsv$", "Precision and recall for the LFMM arms."),
    (r"results/simulation_performance_region\.tsv$", "Region-level performance, all methods."),
    (r"results/simulation_performance_diagnostic\.tsv$", "Marker- and unit-level diagnostic scores (superseded as primary estimand by regions)."),
    (r"results/simulation_method_contrasts(_c1)?\.tsv$", "Paired cluster-bootstrap differences in precision and recall between methods; columns below."),
    (r"results/simulation_diagnostic_contrasts\.tsv$", "Paired contrasts for the diagnostic (marker/unit) scores."),
    (r"results/simulation_qc\.tsv$", "Quality-control record of the simulation pipeline (parser gates, grid completeness, marker-order checks)."),
    (r"results/simulation_truth_sensitivity", "Precision and recall when the detectable-QTN threshold (share of Va) is varied (Figure S-truth-sensitivity)."),
    (r"results/simulation_null_truth", "Null-calibration of reported regions: expected null discoveries vs observed and true false-discovery proportion (Figure S-null-calibration)."),
    (r"results/simulation_stage2_region_details\.(tsv|rds)$", "One row per reported Stage-2 region in every simulation: method, size floor, span, markers, units, truth status, member markers/units, linked QTN."),
    (r"results/simulation_stage2_(size|opportunity)", "Region size vs false-positive status and the opportunity null (Figure S-stage2-size)."),
    (r"results/simulation_floor_", "Size-floor sweep and decomposition: precision/recall by minimum unit size and analysis (Figures S-floor-decomposition)."),
    (r"results/figure_simulation_.*_data\.tsv$", "Data plotted in the named figure."),
    (r"results/floor_stability_", "Floor-stability analysis of simulated regions (labels at half/double floors, density schemes) and its contrasts."),
    (r"results/simulation_bgs_", "Background-selection checks: marker retention, heterozygosity, BGS x aggregation interaction (Figure S-BGS-QC)."),
    (r"results/simulation_bootstrap_", "Bootstrap sensitivity of the pooled estimates (map/burn-in vs crossed resampling)."),
    (r"results/simulation_chrsplit_", "Analyses with chromosomes split into separate units (Table S-chromosome-split)."),
    (r"results/simulation_marker_", "Marker-level redundancy: how many significant markers tag the same QTN; singleton vs multi-marker units."),
    (r"results/simulation_qtn_detection", "Detection of individual detectable QTN by effect size and method (Figure S-QTN-detection)."),
    (r"results/simulation_env_structure_alignment", "Environment-structure alignment of the simulated environments and its relation to false discoveries."),
    (r"results/simulation_lfmm_portability", "Floor effect on LFMM compared with EMMAX."),
    (r"results/simulation_region_granularity", "Region-level vs unit-level scoring contrast."),
    (r"results/simulation_summary_full\.rds$", "All pooled simulation summaries in one R object."),
    (r"results/exploratory_", "Exploratory analysis (LD/Fst vs null burden); not cited as a result."),
]
COLS = {
    "simulation/results/simulation_performance.tsv": [
        ("tag", "BGS treatment: nobgs or bgs"), ("cell", "selection/migration cell, V<QTN variance>_c<dispersal>"),
        ("method", "analysis (emmax_snp_region = marker-wise, emmax_simes_region, emmax_consensus_region, …)"),
        ("n_tested", "hypotheses tested, summed over runs"), ("n_significant", "reported regions"),
        ("TP / FP", "regions linked / not linked to a detectable QTN"), ("FN", "detectable QTN not recovered"),
        ("n_detectable_qtn", "detectable QTN (MAF > 0.1, >= 5% of chromosome Va)"), ("n_recovered", "detectable QTN recovered"),
        ("coverage", "fraction of detectable QTN covered by an eligible Stage-1 unit"),
        ("precision, recall", "pooled ratios"), ("*_ci_lo / *_ci_hi", "95% cluster-bootstrap interval (2,000 resamples of map/burn-in)")],
    "simulation/results/simulation_method_contrasts.tsv": [
        ("tag, cell", "as above"), ("method / reference_method", "compared analyses"),
        ("diff_precision, diff_recall", "method minus reference"), ("*_ci_lo / *_ci_hi", "paired 95% bootstrap interval")],
    "empirical_3sp/10_sensitivity/sensitivity_summary.csv": [
        ("grid", "rho and size floor, e.g. rho=0.50_sf=8"), ("analysis", "consensus, simes or lfmm"),
        ("tested", "units meeting the floor"), ("significant", "units significant after BH"), ("regions", "reported regions"),
        ("on_ecopeak", "regions overlapping an EcoPeak"), ("rotation_fold", "overlap / mean of 10,000 within-chromosome relocations"),
        ("rotation_p", "relocation p (observed placement counted)"), ("perm_p", "locality-preserving permutation p (EMMAX only)"),
        ("realised_fdr", "null/obs ratio: mean significant units under permutation / observed")],
    "empirical_3sp/14_floor_tiers/regions_by_tier.tsv": [
        ("Chr, from, to", "region coordinates (gasAcu1)"), ("n_units, n_markers", "significant units and their markers"),
        ("occupancy", "fraction of markers in the span that belong to these units"), ("statistic", "unit (consensus) or simes"),
        ("floor_lo, floor, floor_hi", "comparison and primary floors"), ("found_lo, found_hi", "region recovered at the low / high floor"),
        ("tier", "floor-stability label (number of floors at which it was found)")],
    "empirical_3sp/14_floor_tiers/ecopeak_by_tier.tsv": [
        ("set", "label (all, 1, 2, 3, >=2)"), ("n", "regions"), ("ecopeak", "overlapping an EcoPeak"),
        ("null_mean", "mean overlaps over relocations"), ("fold", "ecopeak / null_mean"), ("p", "relocation p"), ("statistic", "unit or simes")],
    "empirical_3sp/15_marker_postfilter/marker_postfilter_ecopeak.tsv": [
        ("engine", "LFMM or EMMAX marker-wise scan"), ("floor", "minimum Stage-1 cluster size (1 = no filter)"),
        ("snps_kept, clusters, regions", "retained significant markers, their clusters, assembled regions"),
        ("ecopeak, frac", "regions overlapping an EcoPeak and their fraction"), ("null_mean, fold, p", "relocation test as above")],
    "empirical_9sp/02b_env_structure_alignment/alignment_observed.tsv": [
        ("species, phenotype", "panel and tested phenotype (raw or lineage-adjusted ecotype)"),
        ("structure_group", "null-defining group (pop_locality or lineage)"), ("n_individuals, n_populations, n_structure_groups", "counts"),
        ("grm_axes_r2_population (adjusted)", "R^2 of population-mean phenotype on the leading GRM axes"),
        ("structure_r2_population / _individual_weighted", "R^2 of the phenotype on the structure groups")],
}

def describe(path):
    for pat, d in RULES:
        if re.search(pat, path): return d
    return ""

rows = list(csv.DictReader(open(os.path.join(D, "MANIFEST.tsv")), delimiter="\t"))
missing = [r["file"] for r in rows if not describe(r["file"])]
def size(b):
    b = int(b)
    for u in ("B", "KB", "MB", "GB"):
        if b < 1024 or u == "GB": return f"{b:.0f} {u}" if u == "B" else f"{b:.1f} {u}"
        b /= 1024
total = sum(int(r["bytes"]) for r in rows)

L = ["### Files and variables", "",
     f"{len(rows) + 1} files ({total/1e9:.2f} GB), in three folders that mirror the analysis modules of the code archive "
     "(LDscnR-paper). Point the code at this folder with `LDSCNR_DATA=<this folder>`; the manuscript generators and figure "
     "scripts then read everything they need from it. `MANIFEST.tsv` lists every file with its size and SHA-256.",
     "",
     "**No individual-level genotypes are included.** The stickleback genotypes (Fang et al. 2021) are licensed CC-BY and are "
     "cited, not redistributed: get them from Zenodo doi:10.5281/zenodo.18466057 (exact analysis inputs) or "
     "doi:10.5281/zenodo.4722879 (raw data). The Stage-1 bundles are therefore slim copies (`bundle_slim.rds`) and the "
     "per-individual consensus matrices of the Table 6 cache are removed; the only individual-level matrices are kinship "
     "matrices (GRM, 117 x 117 and 149 x 149, and the 117 x 117 kinship of the Table 6 cache). Table 6 is included as output "
     "(`empirical_3sp/13_goldilocks_units/`, `13b_goldilocks_regions/`), but recomputing it needs the genotypes from Zenodo "
     "doi:10.5281/zenodo.18466057. Sample metadata is limited to a row index, ecotype, population and lineage labels and an "
     "anonymous locality code (`pheno` in `bundle_slim.rds`); the full sample table, with individual IDs, localities and "
     "coordinates, is in the cited Zenodo records.", "",
     "#### Main tables: columns", ""]
for f, cols in COLS.items():
    L += [f"`{f}`", "", "| column | meaning |", "|---|---|"] + [f"| `{c}` | {m} |" for c, m in cols] + [""]
L += ["#### All files", "", "| file | size | description |", "|---|---:|---|"]
for r in sorted(rows, key=lambda r: r["file"]):
    L.append(f"| `{r['file']}` | {size(r['bytes'])} | {describe(r['file'])} |")
L.append("| `MANIFEST.tsv` | | Every file with size (bytes) and SHA-256. |")
section = "\n".join(L) + "\n"

if os.path.exists(OUT):
    txt = open(OUT).read()
    m = re.search(r"^### Files and variables.*?(?=^## |\Z)", txt, re.S | re.M)
    txt = txt[:m.start()] + section + "\n" + txt[m.end():] if m else txt.rstrip() + "\n\n" + section
else:
    txt = "# LDscnR manuscript: analysis outputs\n\n" + section
open(OUT, "w").write(txt)
print(f"{len(rows)} files, {total/1e9:.2f} GB; {len(missing)} without description: {missing[:10]}")
