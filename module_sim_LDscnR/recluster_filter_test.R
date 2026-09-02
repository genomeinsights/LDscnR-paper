## =============================================================================
## recluster_filter_test.R -- does the filter-then-test result survive
## re-clustering at the pilot decay settings?
##
## The main cluster-level result was computed on bundles whose LD_decay was fitted
## with n_win_decay = 10 (18 windows/chr, ~3 Mb spans, 55.8 kb half-decay). The
## newer pilot generation uses n_win_decay = 50 (95 windows/chr, ~0.56 Mb spans,
## 13.3 kb). At rho = 0.5 the half-decay IS the stage-2 clustering distance, so
## the test units themselves differ. That is the untested exposure this addresses.
##
## DESIGN. Both settings are run through the SAME recomputation path -- GDS,
## compute_LD_decay, ld_complexity_reduction, ld_prune_and_eMLG -- so any
## difference is attributable to n_win_decay and not to "recomputed vs stored".
## Running only the new setting and comparing against the stored bundles would
## confound the parameter with the recomputation.
##
## HELD FIXED: the association statistics. emx_p is taken from the bundle and is
## NOT recomputed, so the GRM and EMMAX fit are those of the original pipeline.
## This isolates the effect of decay settings on the UNITS. Rebuilding the GRM
## would change the p-values too and answer a different, larger question.
##
## NOT AVAILABLE HERE: the grm_markers identity check. Elsewhere in this module
## the stage-2 recomputation is required to reproduce the stored GRM marker set
## exactly; that check cannot apply when the clustering is deliberately different.
## Its absence is why this script reports the realised half-decay per bundle --
## a manipulation check that the parameter actually moved.
##
## Env: SIM_DATA, NWIN, CELLS, ENVS, FILES, KS, WINDOWS, ALPHA, CORES, OUT
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})

SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
NWIN  <- as.integer(Sys.getenv("NWIN", "50"))
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5"), ",")[[1]])
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
KS    <- as.integer(strsplit(Sys.getenv("KS", "500,1000,2000,5000,10000,20000"), ",")[[1]])
WIN   <- as.numeric(strsplit(Sys.getenv("WINDOWS", "10,50,100"), ",")[[1]]) * 1000
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
CORES <- as.integer(Sys.getenv("CORES", "1"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

RHO_GRID   <- c(seq(0.05, 0.95, by = 0.05), 0.99)
DECAY_ARGS <- list(min_maf_decay = 0.1, q = 0.95, n_sub_bg = 5000,
                   n_win_decay = NWIN, overlap = 0.5, max_SNPs_decay = Inf,
                   prob_robust = 0.95, max_pairs = 5000, ld_method = "corr",
                   n_strata = 20, keep_el = FALSE, slide = 1000,
                   rho_targets = c(0.99), cores = 1, ld_w_rho = RHO_GRID)
PRUNE_ARGS <- list(ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
                   score_threshold = 0.80, min_r2_rho = 0.5,
                   distance_threshold = 1e5, compute_unflagged_eMLG = FALSE)
CR_RHO <- 0.5

units_for <- function(cell, tag, env, i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, tag, i, cell, env)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); m <- as.data.table(x$map)
  if (!"emx_p" %in% names(m)) return(NULL)
  gp <- tempfile(fileext = ".gds"); on.exit(unlink(gp), add = TRUE)
  gds <- create_gds_from_geno(geno = x$GTs, map = m, gds_path = gp)
  D   <- do.call(compute_LD_decay, c(list(gds = gds, el_data_folder = NULL), DECAY_ARGS))
  s1  <- ld_complexity_reduction(map = m, LD_decay = D, rho = CR_RHO, gds = gds)
  pr  <- do.call(ld_prune_and_eMLG,
                 c(list(GTs = x$GTs, stage1 = s1, LD_decay = D, cores = 1), PRUNE_ARGS))
  ## manipulation check: did the parameter actually move the fit?
  ds  <- as.data.table(D$decay_sum)
  hk  <- median(1/ds$a[is.finite(ds$a) & ds$a > 0]) / 1000

  m[, ld_w_new := D$ld_ws[m$marker, "rho_0.95"]]
  g  <- as.data.table(pr$groups)
  ms <- data.table(marker = unlist(g$members, use.names = FALSE), CL_id = rep.int(g$group_id, lengths(g$members)))
  m <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL_id)]
  drv <- m[true_QTN %in% TRUE & MAF > 0.1 & p_Va > 0.05]
  m[, d_qtn := Inf]; m[, near_qtn := NA_character_]
  if (nrow(drv)) for (ch in unique(m$Chr)) {
    dd <- drv[Chr == ch]; if (!nrow(dd)) next
    ii <- which(m$Chr == ch)
    Dm <- abs(outer(m$Pos[ii], dd$Pos, "-")); j <- max.col(-Dm)
    m[ii, `:=`(d_qtn = Dm[cbind(seq_along(ii), j)],
               near_qtn = paste0(cell, "_", i, "_", ch, "_", dd$Pos[j]))]
  }
  m <- m[!(true_QTN %in% TRUE)]
  if (!nrow(m)) return(NULL)
  rep_mk <- intersect(pr$pruned, m$marker)
  u <- m[, .(ld_w = median(ld_w_new, na.rm = TRUE), MAF = median(MAF, na.rm = TRUE),
             n_loci = .N, d_qtn = min(d_qtn), near_qtn = near_qtn[which.min(d_qtn)]), by = CL_id]
  rp <- m[marker %in% rep_mk, .(CL_id, p = emx_p)][, .SD[1], by = CL_id]
  u  <- merge(u, rp, by = "CL_id")[is.finite(p)]
  if (!nrow(u)) return(NULL)
  u[, `:=`(CL_id = paste0(i, "_", CL_id), half_kb = hk)]
  u[]
}

score_sel <- function(u, idx, w) {
  s <- u[idx]
  if (!nrow(s)) return(data.table(n_sel = 0L, n_sig = 0L, n_tp = 0L, qtn_hit = 0L))
  q <- p.adjust(s$p, "BH"); sig <- which(q < ALPHA)
  if (!length(sig)) return(data.table(n_sel = nrow(s), n_sig = 0L, n_tp = 0L, qtn_hit = 0L))
  h <- s[sig]; tp <- h$d_qtn < w
  data.table(n_sel = nrow(s), n_sig = length(sig), n_tp = sum(tp),
             qtn_hit = uniqueN(h$near_qtn[tp]))
}

grid <- CJ(cell = CELLS, tag = c("nobgs","bgs"), env = ENVS, sorted = FALSE)
cat(sprintf("  n_win_decay = %d | %d panels x %d files, CORES=%d\n",
            NWIN, nrow(grid), length(FILES), CORES))

one_panel <- function(r) {
  cell <- grid$cell[r]; tag <- grid$tag[r]; env <- grid$env[r]
  u <- tryCatch(rbindlist(lapply(FILES, function(i) units_for(cell, tag, env, i)), fill = TRUE),
                error = function(e) { message("  FAIL ", cell," ",tag," env",env,": ",conditionMessage(e)); NULL })
  if (is.null(u) || !nrow(u)) return(NULL)
  set.seed(9000 + r)
  sels <- list(ld_w = order(-u$ld_w), MAF = order(-u$MAF),
               size = order(-u$n_loci), random = sample.int(nrow(u)))
  meta <- data.table(cell, tag, env, nwin = NWIN, n_units = nrow(u),
                     half_kb = round(median(u$half_kb), 1))
  out <- rbindlist(lapply(WIN, function(w) {
    base <- cbind(meta, method = "genome_wide", k = nrow(u), window_kb = w/1000,
                  score_sel(u, seq_len(nrow(u)), w))
    per <- rbindlist(lapply(KS, function(kk) {
      if (kk >= nrow(u)) return(NULL)
      rbindlist(lapply(names(sels), function(nm)
        cbind(meta, method = nm, k = kk, window_kb = w/1000, score_sel(u, head(sels[[nm]], kk), w))))
    }))
    rbind(base, per, fill = TRUE)
  }))
  cat(sprintf("    %-9s %-5s env%-3d  %6d units, half-decay %5.1f kb\n",
              cell, tag, env, nrow(u), meta$half_kb))
  out
}
res <- if (CORES > 1) {
  mclapply(seq_len(nrow(grid)), one_panel, mc.cores = CORES, mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(grid)), one_panel)
}
all <- rbindlist(Filter(Negate(is.null), res), fill = TRUE)
stopifnot(nrow(all) > 0)
fn <- file.path(OUT, sprintf("recluster_filter_test_nwin%d.csv", NWIN))
fwrite(all, fn)
cat(sprintf("\n  realised half-decay: median %.1f kb (range %.1f-%.1f)\n",
            median(all$half_kb), min(all$half_kb), max(all$half_kb)))
cat(sprintf("  median units/panel: %.0f\n", median(all$n_units)))
cat(sprintf("  written: %s (%d rows)\n", fn, nrow(all)))
