## =============================================================================
## module_sim_LDscnR / floor_sweep.R
##
## PK's open decision: 2c gets 1 discovery at minimum cluster size 2 on the
## stickleback panel and 14 at size 8, and the question is whether the higher
## floor DETECTS more or merely escapes multiplicity. BH over 2,631 units
## against 790,578 markers is a 300x smaller correction, so more discoveries at
## a higher floor is exactly what relief would produce -- with 95% of markers
## untestable as an unmeasured recall cost.
##
## THE SIMS CAN MEASURE THE COST THE PANEL CANNOT. QTN positions are known, so
## every floor gets precision AND recall, not just a discovery count. If the
## panel gains discoveries as the floor rises while the sims lose them, that
## asymmetry is what multiplicity relief predicts -- the relief scales with the
## marker count being escaped, and the panel's is 790k against ~30k per
## chromosome here.
##
## The floor is applied BEFORE BH, so the multiplicity actually changes; n_test
## is recorded so the relief is visible rather than inferred. Scoring follows
## snp_vs_cluster_dedup.R exactly -- one region per QTN, best tagger, satellites
## removed -- so the numbers are comparable to that table.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/floor_sweep.R
## Env: SIM_DATA, CELLS, TAGS, ENVS, FILES, FLOORS, ALPHA, CORES, OUT
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
SIM    <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT    <- Sys.getenv("OUT", "module_sim_LDscnR/results/floor_sweep")
CELLS  <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]
TAGS   <- strsplit(Sys.getenv("TAGS", "nobgs,bgs"), ",")[[1]]
ENVS   <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
FILES  <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
FLOORS <- as.integer(strsplit(Sys.getenv("FLOORS", "1,2,3,5,8,12"), ",")[[1]])
ALPHA  <- as.numeric(Sys.getenv("ALPHA", "0.05"))
CORES  <- as.integer(Sys.getenv("CORES", "8"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
simes <- function(p) { p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_); min(n * p / seq_len(n)) }

one_panel <- function(CELL, TAG, ENV) {
  clu <- list(); snp <- list(); lnk <- list()
  for (i in FILES) {
    f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV)
    if (!file.exists(f)) next
    x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
    pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
            LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
            score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
            compute_unflagged_eMLG = TRUE, cores = 1)
    g  <- as.data.table(pr$groups)
    ms <- rbindlist(lapply(seq_len(nrow(g)), function(k)
            data.table(marker = g$members[[k]], CL = paste0(i, "_", g$group_id[k]))))
    mm <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL)]
    th  <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
                            rho_r2 = 0.75, rho_d = 0.95, dmax_cap = 1e5)
    drv <- mm[true_pos_QTN %in% TRUE]
    if (nrow(drv)) lnk[[length(lnk)+1]] <- rbindlist(lapply(seq_len(nrow(drv)), function(j) {
      ch <- as.character(drv$Chr[j])
      near <- mm[as.character(Chr) == ch & abs(Pos - drv$Pos[j]) < th$dmax]
      if (!nrow(near)) return(NULL)
      r2 <- suppressWarnings(cor(x$GTs[, drv$marker[j]], x$GTs[, near$marker],
                                 use = "pairwise.complete.obs")^2)
      d <- data.table(CL = near$CL, r2 = as.numeric(r2))[is.finite(r2) & r2 >= th$r2min]
      if (!nrow(d)) return(NULL)
      d[, .(r2 = max(r2)), by = CL][, qtn := paste0(i, "_", drv$marker[j])][] }))
    clu[[length(clu)+1]] <- mm[, .(p_simes = simes(emx_p), n_loci = .N,
                                   has_qtn = any(true_pos_QTN %in% TRUE)), by = CL]
    snp[[length(snp)+1]] <- mm[, .(CL, p = emx_p)]
  }
  su <- rbindlist(clu, fill = TRUE); sn <- rbindlist(snp, fill = TRUE)
  lk <- rbindlist(lnk, fill = TRUE); nq <- sum(su$has_qtn)
  if (!nq) return(NULL)

  score <- function(fl, n_test, lab) {
    sub  <- lk[CL %in% fl]
    best <- if (nrow(sub)) sub[order(-r2)][, .SD[1], by = qtn] else sub
    kept <- setdiff(fl, setdiff(sub$CL, best$CL))
    tp_reg <- if (nrow(best)) uniqueN(best$CL)  else 0L
    tp_qtn <- if (nrow(best)) uniqueN(best$qtn) else 0L
    data.table(arm = lab, n_test = n_test, flagged = length(fl),
               dedup_regions = length(kept), dedup_tp = tp_reg,
               prec = tp_reg / max(length(kept), 1), rec = tp_qtn / nq)
  }
  out <- rbindlist(lapply(FLOORS, function(fo) {
    ok <- su$n_loci >= fo & is.finite(su$p_simes)
    q  <- rep(NA_real_, nrow(su)); q[ok] <- p.adjust(su$p_simes[ok], "BH")
    score(su$CL[which(!is.na(q) & q < ALPHA)], sum(ok), paste0("simes_floor", fo)) }))
  ## single-SNP reference: the multiplicity the floor is escaping
  qs <- p.adjust(sn$p, "BH")
  out <- rbind(out, score(unique(sn$CL[which(qs < ALPHA)]), nrow(sn), "snp"))
  out[, `:=`(cell = CELL, tag = TAG, env = ENV, n_qtn = nq, n_clusters = nrow(su))][]
}

grid <- CJ(cell = CELLS, tag = TAGS, env = ENVS, sorted = FALSE)
cat(sprintf("%d panels, CORES=%d\n", nrow(grid), CORES))
res <- rbindlist(Filter(Negate(is.null), mclapply(seq_len(nrow(grid)), function(z)
  tryCatch(one_panel(grid$cell[z], grid$tag[z], grid$env[z]),
           error = function(e) { cat("FAIL", z, conditionMessage(e), "\n"); NULL }),
  mc.cores = CORES)))
fwrite(res, file.path(OUT, "floor_sweep_raw.csv"))
cat(sprintf("\n%d panel-arms from %d panels\n", nrow(res), uniqueN(res[, .(cell,tag,env)])))
