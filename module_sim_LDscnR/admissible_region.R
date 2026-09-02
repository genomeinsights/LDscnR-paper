## =============================================================================
## module_sim_LDscnR / admissible_region.R
##
## Map the JOINT admissible region of ld_w_threshold x distance_threshold.
##
## The two are coupled by construction (ldw_dcap_coupling.R): dcap splits a run
## at a gap between consecutive FLAGGED stage-1 clusters, and ld_w decides which
## clusters are flagged. So the pair has a joint admissible region rather than
## two independent ranges, and a factorial grid over them contains cells that
## should not be run at all.
##
## ADMISSIBILITY IS BOUNDED ON BOTH SIDES, and both bounds are degeneracies of
## stage 2 rather than merely bad settings:
##
##   UNDER-SPLIT   no gap on the chromosome exceeds dcap, so the chromosome is a
##                 SINGLE RUN. Merging is unconstrained by distance and whatever
##                 structure stage 2 exists to preserve dissolves -- the failure
##                 mode the panel session identified for Eda.
##   OVER-SPLIT    every gap exceeds dcap, so every flagged cluster is its own
##                 run. No merge is ever possible and STAGE 2 IS A NO-OP: the
##                 output is stage 1 with extra steps.
##   DROPPED       fewer than two flagged clusters on the chromosome, so it
##                 leaves the analysis entirely rather than failing loudly.
##
## A cell is admissible for a chromosome when none of the three holds. The
## reported quantity is the fraction of chromosomes admissible, plus the merge
## opportunity -- the share of flagged clusters sitting in a run with at least
## one other cluster, which is what stage 2 can actually act on.
##
## THE TWO CRITERIA PULL IN OPPOSITE DIRECTIONS, which is why the region has an
## interior rather than an edge: admissibility wants SMALL dcap (large dcap
## under-splits) and merge opportunity wants LARGE dcap (small dcap isolates
## clusters). The canonical point sits near the corner -- 100% admissible with
## 98% merge opportunity -- which is a justification it did not previously have.
##
## RESULT: the region MOVES BETWEEN SIMULATION SETTINGS OF THE SAME SIMULATOR.
## The widest dcap holding 95% of chromosomes spans 50 to 1,000 kb across the
## four cells at ld_w = 0.0125, and 25 to 1,000 kb at ld_w = 0.05 where one cell
## has no admissible dcap at all. Only ld_w = 0.025 is 100% admissible in all
## four cells. ld_w = 0.10 is inadmissible everywhere (35-72% by cell), because
## 27.5% of chromosomes retain fewer than two flagged clusters and leave the
## analysis silently.
##
## UNDER-SPLITTING IS THE ONLY DEGENERACY THAT MATTERS IN PRACTICE: 21-46% of
## chromosome-cells against 0-6% over-split. Stage 2 becoming a no-op is a
## theoretical risk; the chromosome dissolving into one run is a real one.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/admissible_region.R
## Env: SIM_DATA, CELLS, TAGS, ENV, FILES, LDW, DCAP, CORES, OUT
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)
                  library(ggplot2)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/operating_points")
FIG   <- Sys.getenv("FIG", "module_sim_LDscnR/figures")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAGS  <- strsplit(Sys.getenv("TAGS", "nobgs,bgs"), ",")[[1]]
ENV   <- as.integer(Sys.getenv("ENV", "1"))
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
LDW   <- as.numeric(strsplit(Sys.getenv("LDW", "0.00625,0.0125,0.025,0.05,0.1"), ",")[[1]])
DCAP  <- as.numeric(strsplit(Sys.getenv("DCAP", "1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2.5e6,5e6"), ",")[[1]])
CORES <- as.integer(Sys.getenv("CORES", "8"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

jobs <- CJ(CELL = CELLS, TAG = TAGS, i = FILES, sorted = FALSE)
one <- function(k) {
  r <- jobs[k]
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, r$TAG, r$i, r$CELL, ENV)
  if (!file.exists(f)) return(NULL)
  x  <- readRDS(f)
  ms <- as.data.table(x$complexity_reduction$stage1$map_snp)
  ds <- as.data.table(x$LD_decay$decay_sum)
  rbindlist(lapply(LDW, function(t) {
    ## flagging matches ld_prune_and_eMLG.R:389 exactly -- any member, strict >
    ids <- ms[ld_w_095 > t, unique(CL_id)]
    ext <- ms[CL_id %in% ids, .(Chr = Chr[1], pmin = min(Pos), pmax = max(Pos)), by = CL_id]
    rbindlist(lapply(ds$Chr, function(ch) {
      e <- ext[Chr == ch][order(pmin)]
      nf <- nrow(e)
      if (nf < 2) return(data.table(CELL=r$CELL, TAG=r$TAG, i=r$i, Chr=ch, ldw=t,
                                    dcap=DCAP, n_flag=nf, n_runs=NA_integer_,
                                    merge_opp=NA_real_, state="dropped"))
      g <- e$pmin[-1] - e$pmax[-nf]
      rbindlist(lapply(DCAP, function(dc) {
        split_after <- g > dc
        run <- cumsum(c(TRUE, split_after))
        sz  <- tabulate(run)
        data.table(CELL=r$CELL, TAG=r$TAG, i=r$i, Chr=ch, ldw=t, dcap=dc,
                   n_flag=nf, n_runs=length(sz),
                   ## share of flagged clusters that can actually be merged
                   merge_opp = sum(sz[sz >= 2]) / nf,
                   state = if (length(sz) == 1L) "under_split"
                           else if (length(sz) == nf) "over_split" else "admissible")
      })) })) })) }

G <- rbindlist(Filter(Negate(is.null), mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES)))
fwrite(G, file.path(OUT, "admissible_region_by_chr.csv"))

S <- G[, .(n_chr = .N,
           pct_admissible = 100 * mean(state == "admissible"),
           pct_under      = 100 * mean(state == "under_split"),
           pct_over       = 100 * mean(state == "over_split"),
           pct_dropped    = 100 * mean(state == "dropped"),
           merge_opp      = mean(merge_opp, na.rm = TRUE),
           med_runs       = as.double(median(n_runs, na.rm = TRUE))),
       by = .(ldw, dcap)][order(ldw, dcap)]
fwrite(S, file.path(OUT, "admissible_region.csv"))

cat("== fraction of chromosomes ADMISSIBLE (neither degeneracy), % of 160\n")
print(dcast(S, ldw ~ dcap, value.var = "pct_admissible"))
cat("\n== MERGE OPPORTUNITY: share of flagged clusters in a run with >= 1 other\n")
print(dcast(S, ldw ~ dcap, value.var = "merge_opp"))
cat("\n== the canonical operating point\n")
print(S[ldw == 0.025 & dcap == 1e5])

## Does the region MOVE between simulation settings? This is the claim that
## makes it dataset-specific rather than a package constant.
B <- G[, .(pct = mean(state == "admissible")), by = .(CELL, ldw, dcap)]
Wd <- B[pct >= 0.95, .(max_ok_kb = max(dcap)/1e3), by = .(CELL, ldw)]
fwrite(dcast(Wd, ldw ~ CELL, value.var = "max_ok_kb"),
       file.path(OUT, "admissible_region_by_cell.csv"))
cat("\n== widest dcap holding >=95% of chromosomes, by simulation cell (kb)\n")
print(dcast(Wd, ldw ~ CELL, value.var = "max_ok_kb"))
cat("   NA = no tested dcap reaches 95% for that cell\n")
cat("\n== which degeneracy dominates (pooled over dcap)\n")
print(G[, .(under_split = 100*mean(state=="under_split"),
            over_split  = 100*mean(state=="over_split"),
            dropped     = 100*mean(state=="dropped")), by = ldw][order(ldw)])

P <- melt(S, id.vars = c("ldw","dcap"),
          measure.vars = c("pct_admissible","merge_opp"))
P[variable == "merge_opp", value := 100 * value]
P[, variable := factor(variable, c("pct_admissible","merge_opp"),
    c("chromosomes admissible (%)", "merge opportunity (% of flagged clusters)"))]
g <- ggplot(P, aes(factor(dcap/1e3), factor(ldw), fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f", value)), size = 2.6, colour = "grey15") +
  facet_wrap(~ variable, ncol = 1) +
  scale_fill_gradient(low = "#f7f7f7", high = "#2c7fb8", limits = c(0, 100), name = "%") +
  annotate("point", x = "100", y = "0.025", shape = 21, size = 7,
           stroke = 1.1, colour = "#c1272d", fill = NA) +
  labs(x = "distance_threshold (kb)", y = "ld_w_threshold",
       title = "Joint admissible region of the two coupled thresholds",
       subtitle = "red ring: the canonical operating point (ld_w 0.025, dcap 100 kb)") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(), panel.grid = element_blank())
ggsave(file.path(FIG, "admissible_region.png"), g, width = 8.2, height = 5.4, dpi = 200)
cat(sprintf("\nfigure: %s\n", file.path(FIG, "admissible_region.png")))
