## =============================================================================
## or_stage2_regions.R -- outlier regions that are the PRE-DEFINED stage-2
## clusters, scored by the regular OR convention.
##
## In the standard pipeline a region is CONSTRUCTED from the significant markers:
## ld_edges joins them by r2 and distance, ld_regions cuts the graph. That makes
## the region set a function of which markers happened to be significant, and it
## introduces two further thresholds (rho_ld, dcap).
##
## Here the candidate regions are fixed BEFORE any phenotype is seen -- they are
## the stage-2 LD clusters -- and a region is REPORTED when it contains at least
## l_min significant markers. Truth is then applied exactly as for constructed
## regions: evaluate_ors against qtn_ld_table, dedup so one region per QTN counts
## once and any further region on the same QTN is a false positive.
##
## Consequences worth being explicit about:
##   - region construction contributes no thresholds; rho_ld and dcap drop out
##     of the region step (dcap still enters the TP-matching window via
##     score_thresholds)
##   - the region set is IDENTICAL across arms, so arms differ only in which
##     regions they report -- a cleaner paired comparison than the standard
##     pipeline, where the arms build different regions
##   - a QTN in no stage-2 cluster is unreachable by construction; that ceiling
##     is reported as n_reachable
##
## k counts CLUSTERS selected, not SNPs.
## Env: SIM_DATA, CELLS, OUT, KS, LMINS, ALPHA, DCAP, ENGINE, CORES
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})

SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
CELLS <- strsplit(Sys.getenv("CELLS", ""), ",")[[1]]
KS    <- as.integer(strsplit(Sys.getenv("KS", "1000,5000,20000,50000"), ",")[[1]])
LMINS <- as.integer(strsplit(Sys.getenv("LMINS", "1,2,3"), ",")[[1]])
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
DCAP  <- as.numeric(Sys.getenv("DCAP", "1e5"))
ENG   <- Sys.getenv("ENGINE", "emmax")
CORES <- as.integer(Sys.getenv("CORES", "1"))
## MERGE = "ld" adds a SECOND clustering over the stage-2 clusters that were
## reported. Two stage-2 clusters tagging the same QTN are each within dmax of
## it, hence within 2*dmax of each other, so the merge cap defaults to twice the
## scoring cap. Without this step such clusters are separate regions and all but
## the best are false positives by dedup -- on one panel 115 clusters tag a QTN
## and only 21 survive, so 94 of the false positives are satellites of a locus
## that was already found.
MERGE      <- Sys.getenv("MERGE", "none")
MERGE_RHO  <- as.numeric(Sys.getenv("MERGE_RHO", "0.75"))
MERGE_DCAP <- as.numeric(Sys.getenv("MERGE_DCAP", as.character(2 * DCAP)))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
pcol  <- if (ENG == "emmax") "emx_p" else "lfmm_p"
if (!length(CELLS) || !nzchar(CELLS[1]))
  CELLS <- as.vector(outer(c("V0.5_c1","V0.5_c2","V1_c1.5","V2_c1"), paste0("_env", 1:10), paste0))

## pool 10 map sets AND their stage-2 partitions
pool2 <- function(tag, cell) {
  fs <- list.files(SIM, full.names = TRUE,
                   pattern = sprintf("^adapt_%s_chr[0-9]+_%s[.]rds$", tag, gsub("\\.", "[.]", cell)))
  if (!length(fs)) return(NULL)
  fs <- fs[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(fs))))]
  maps <- gts <- ldws <- decs <- grps <- prns <- vector("list", length(fs))
  for (i in seq_along(fs)) {
    d <- readRDS(fs[i]); m <- as.data.table(d$map)
    pr <- ld_prune_and_eMLG(GTs = d$GTs, stage1 = d$complexity_reduction$stage1,
            LD_decay = d$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
            score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
            compute_unflagged_eMLG = FALSE, cores = 1)
    if (!identical(sort(pr$pruned), sort(d$grm_markers)))
      stop("stage-2 does not reproduce grm_markers for ", basename(fs[i]))
    g <- as.data.table(pr$groups)
    pfx <- paste0("R", i, "_")
    grps[[i]] <- lapply(seq_len(nrow(g)), function(z) paste0(pfx, g$members[[z]]))
    prns[[i]] <- paste0(pfx, pr$pruned)
    m[, `:=`(Chr = paste0(pfx, Chr), marker = paste0(pfx, marker))]
    G <- d$GTs; colnames(G) <- m$marker
    lw <- d$ld_ws; rownames(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0(pfx, Chr)]
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
  }
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE))
  list(map = map, GTs = do.call(cbind, gts)[, map$marker],
       ld_ws = do.call(rbind, ldws)[map$marker, ],
       decay_sum = rbindlist(decs, fill = TRUE), groups = do.call(c, grps),
       pruned = unlist(prns))
}

one <- function(cell, tag) {
  P <- tryCatch(pool2(tag, cell), error = function(e) { message("  FAIL ", cell, " ", tag, ": ", conditionMessage(e)); NULL })
  if (is.null(P)) return(NULL)
  map <- P$map; p <- map[[pcol]]
  if (is.null(p) || all(is.na(p))) return(NULL)
  n_true <- sum(map$true_pos_QTN %in% TRUE); if (!n_true) return(NULL)
  th  <- score_thresholds(P$decay_sum, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = DCAP)
  grp <- P$groups
  prn <- P$pruned            # the real stage-2 representatives
  ldw <- P$ld_ws[, "rho_0.95"]; names(ldw) <- map$marker
  af  <- colMeans(P$GTs, na.rm = TRUE)/2; maf <- pmin(af, 1-af); names(maf) <- map$marker
  set.seed(7100 + nchar(cell) + nchar(tag))

  sig_of <- function(idx) { if (!length(idx)) return(character(0))
    q <- p.adjust(p[idx], "BH"); map$marker[idx][which(q < ALPHA)] }

  ## SELECTION IS ON CLUSTERS, NOT SNPS. ld_w is a local-LD statistic and is
  ## autocorrelated along the genome, so ranking SNPs by it returns a few large
  ## blocks many times over -- measured on one panel, the top 1,000 SNPs are 13
  ## clusters, while the top 1,000 clusters are 21,752 SNPs spread genome-wide.
  ## Ranking clusters gives one value per block, which is the level at which ld_w
  ## carries independent information. An earlier version ranked SNPs and lost 10
  ## QTN-bearing clusters at k=5,000 where cluster-ranking loses 1.
  ## marker -> group id, built VECTORISED. An earlier version created a named
  ## vector over all 294k markers and assigned into it once per group, which is
  ## a name lookup over the whole vector 100k times -- quadratic, and the process
  ## was killed before finishing. lengths()/rep.int does it in one pass.
  gmap <- data.table(marker = unlist(grp, use.names = FALSE),
                     g = rep.int(seq_along(grp), lengths(grp)))
  cl_stat <- merge(gmap, data.table(marker = map$marker, ldw = ldw, maf = maf),
                   by = "marker")
  cl_rank <- cl_stat[, .(ldw = median(ldw, na.rm = TRUE),
                         maf = median(maf, na.rm = TRUE), n = .N), by = g]
  mk_of  <- function(gs) cl_stat[g %in% gs]$marker
  idx_of <- function(mk) which(map$marker %in% mk)
  arms <- list(alpha = sig_of(seq_len(nrow(map))))
  for (kk in KS) { if (kk >= nrow(cl_rank)) next
    arms[[paste0("ld_w_",   kk)]] <- sig_of(idx_of(mk_of(head(cl_rank[order(-ldw)]$g, kk))))
    arms[[paste0("MAF_",    kk)]] <- sig_of(idx_of(mk_of(head(cl_rank[order(-maf)]$g, kk))))
    arms[[paste0("size_",   kk)]] <- sig_of(idx_of(mk_of(head(cl_rank[order(-n)]$g,   kk))))
    arms[[paste0("random_", kk)]] <- sig_of(idx_of(mk_of(cl_rank$g[sample.int(nrow(cl_rank), kk)]))) }

  ## truth table over every marker that any stage-2 cluster could report
  uni  <- unique(unlist(grp, use.names = FALSE))
  qtab <- qtn_ld_table(P$GTs, map, uni, 2e6, cores = 1)
  ## Ceiling: how many QTN ANY stage-2 cluster could reach. Derived from the
  ## truth table rather than by calling evaluate_ors on ~100k regions, which is
  ## the call that exhausted memory on the first attempt. qtn_ld_table returns
  ## columns marker / qtn_marker / r2 / dist_bp -- named explicitly, because an
  ## earlier version guessed them and would have missed dist_bp, silently
  ## dropping the distance constraint and overcounting.
  qt <- as.data.table(qtab)
  stopifnot(all(c("marker","qtn_marker","r2","dist_bp") %in% names(qt)))
  ## qtn_ld_table keys on type == "QTN" (ALL QTN), while n_true and evaluate_ors
  ## use flag_true_qtns' true_pos_QTN (MAF >= 0.1, within-chromosome p_Va >= 0.05).
  ## Counting reachability over the wider set made n_reachable exceed n_true,
  ## which is impossible. Restricted to the same set evaluate_ors scores against.
  drv_mk <- map[true_pos_QTN %in% TRUE]$marker
  n_reach <- uniqueN(qt[r2 >= th$r2min & abs(dist_bp) < th$dmax &
                        marker %in% uni & qtn_marker %in% drv_mk]$qtn_marker)
  stopifnot(n_reach <= n_true)

  res <- rbindlist(lapply(names(arms), function(nm) {
    mk <- arms[[nm]]
    parts <- strsplit(nm, "_(?=[0-9]+$)", perl = TRUE)[[1]]
    meth <- parts[1]; kv <- if (length(parts) > 1) as.integer(parts[2]) else NA_integer_
    hit <- if (length(mk)) vapply(grp, function(g) sum(g %in% mk), integer(1)) else integer(length(grp))
    rbindlist(lapply(LMINS, function(L) {
      r  <- grp[hit >= L]
      n_before <- length(r)
      if (MERGE == "ld" && length(r) > 1) {
        ## Join reported clusters through their pruned representatives, then
        ## expand each merged component back to every member marker of the
        ## clusters it absorbed.
        reps <- vapply(r, function(g) {
          z <- g[g %in% prn]; if (length(z)) z[1] else g[1]
        }, character(1))
        ok <- reps %in% map$marker
        if (sum(ok) > 1) {
          rk   <- reps[ok]
          e2   <- ld_edges(rk, P$GTs, map[, .(marker, Chr, Pos)], P$decay_sum,
                           rho_ld = MERGE_RHO, dcap = MERGE_DCAP)
          comp <- ld_regions(rk, e2)
          cid  <- rep(NA_integer_, length(rk))
          for (z in seq_along(comp)) cid[rk %in% comp[[z]]] <- z
          miss <- which(is.na(cid))
          if (length(miss)) cid[miss] <- length(comp) + seq_along(miss)
          merged <- unname(lapply(split(which(ok), cid), function(ii) unique(unlist(r[ii]))))
          r <- c(merged, r[!ok])
        }
      }
      ev <- if (length(r)) evaluate_ors(r, map, qtab, th$r2min, th$dmax)
            else list(Precision = NA_real_, Recall = 0)
      data.table(method = meth, k = kv, l_min = L, n_sig = length(mk),
                 n_regions = length(r), n_before_merge = n_before,
                 precision = ev$Precision, recall = ev$Recall,
                 PR = (if (is.na(ev$Precision)) 0 else ev$Precision) * ev$Recall)
    }))
  }))
  cat(sprintf("    %-16s %-5s  %6d clusters, %2d QTN (%d reachable), alpha calls %d\n",
              cell, tag, length(grp), n_true, n_reach, length(arms$alpha)))
  cbind(cell, tag, n_true, n_reachable = n_reach, n_clusters = length(grp), res)
}

grid <- CJ(cell = CELLS, tag = c("nobgs","bgs"), sorted = FALSE)
cat(sprintf("  %d panels, engine %s, dcap %g, CORES=%d\n", nrow(grid), ENG, DCAP, CORES))
run <- function(i) one(grid$cell[i], grid$tag[i])
res <- if (CORES > 1) {
  mclapply(seq_len(nrow(grid)), run, mc.cores = CORES, mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(grid)), run)
}
all <- rbindlist(Filter(Negate(is.null), res), fill = TRUE)
stopifnot(nrow(all) > 0)
fn <- file.path(OUT, sprintf("or_stage2_regions_%s.csv", ENG))
fwrite(all, fn)
cat(sprintf("\n  written: %s (%d rows, %d panels)\n", fn, nrow(all), uniqueN(all[, .(cell,tag)])))
