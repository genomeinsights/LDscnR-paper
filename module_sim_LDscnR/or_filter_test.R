## =============================================================================
## or_filter_test.R -- filter-then-test scored through the OUTLIER-REGION
## pipeline, so the numbers land on the same scale as the C-score benchmark.
##
## Everything else in this module scores filtering at the level of markers or
## clusters. The C-vs-alpha work scored REGIONS: significant markers are joined
## into outlier regions by ld_edges/ld_regions and evaluated by evaluate_ors
## against qtn_ld_table, with one region per QTN after dedup. Those two families
## of numbers are not comparable. This script runs the filter through the region
## machinery unchanged, so "filtered vs unfiltered" can be read directly against
## "C vs alpha".
##
## DESIGN. Per pooled panel, one candidate union, one ld_edges, one qtn_ld_table;
## each arm is then scored by subsetting. Arms differ ONLY in which markers are
## eligible and hence in the BH denominator:
##
##   alpha   BH over all markers                     (the conventional analysis)
##   ld_w    BH over the top k markers by ld_w       (phenotype-blind selection)
##   MAF     BH over the top k by MAF                (control: is it frequency?)
##   random  BH over k at random                     (control: is it just size?)
##
## Reported at the OPERATING POINT (alpha = 0.05), not as PR-AUC. Integrating
## over k would credit the method for operating points nobody deploys, which is
## how the C-score came to look good at a fitted tau.
##
## Env: SIM_DATA, CELLS, OUT, KS, LMINS, ALPHA, RHO_LD, DCAP, ENGINE, CORES
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})

SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
CELLS <- strsplit(Sys.getenv("CELLS", ""), ",")[[1]]
KS    <- as.integer(strsplit(Sys.getenv("KS", "1000,5000,20000,50000"), ",")[[1]])
LMINS <- as.integer(strsplit(Sys.getenv("LMINS", "1,3"), ",")[[1]])
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
RHO_LD<- as.numeric(Sys.getenv("RHO_LD", "0.75"))
DCAP  <- as.numeric(Sys.getenv("DCAP", "1e5"))
ENG   <- Sys.getenv("ENGINE", "emmax")
CORES <- as.integer(Sys.getenv("CORES", "1"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
pcol  <- if (ENG == "emmax") "emx_p" else "lfmm_p"

if (!length(CELLS) || !nzchar(CELLS[1]))
  CELLS <- as.vector(outer(c("V0.5_c1","V0.5_c2","V1_c1.5","V2_c1"),
                           paste0("_env", 1:10), paste0))

pool <- function(tag, cell) {
  fs <- list.files(SIM, full.names = TRUE,
                   pattern = sprintf("^adapt_%s_chr[0-9]+_%s[.]rds$", tag, gsub("\\.", "[.]", cell)))
  if (!length(fs)) return(NULL)
  fs <- fs[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(fs))))]
  maps <- gts <- ldws <- decs <- vector("list", length(fs))
  for (i in seq_along(fs)) {
    d <- readRDS(fs[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- m$marker
    lw <- d$ld_ws; rownames(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
  }
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE))
  list(map = map, GTs = do.call(cbind, gts)[, map$marker],
       ld_ws = do.call(rbind, ldws)[map$marker, ], decay_sum = rbindlist(decs, fill = TRUE))
}

one <- function(cell, tag) {
  P <- tryCatch(pool(tag, cell), error = function(e) NULL)
  if (is.null(P)) return(NULL)
  map <- P$map; p <- map[[pcol]]
  if (is.null(p) || all(is.na(p))) return(NULL)
  n_true <- sum(map$true_pos_QTN %in% TRUE); if (!n_true) return(NULL)
  th  <- score_thresholds(P$decay_sum, rho_r2 = RHO_LD, rho_d = 0.95, dmax_cap = DCAP)
  ldw <- P$ld_ws[, "rho_0.95"]
  af  <- colMeans(P$GTs, na.rm = TRUE)/2; maf <- pmin(af, 1-af)
  set.seed(7000 + nchar(cell) + nchar(tag))

  ## which markers each arm declares significant
  sig_of <- function(idx) {
    if (!length(idx)) return(character(0))
    q <- p.adjust(p[idx], "BH")
    map$marker[idx][which(q < ALPHA)]
  }
  arms <- list(alpha = sig_of(seq_len(nrow(map))))
  for (kk in KS) {
    if (kk >= nrow(map)) next
    arms[[paste0("ld_w_",   kk)]] <- sig_of(head(order(-ldw), kk))
    arms[[paste0("MAF_",    kk)]] <- sig_of(head(order(-maf), kk))
    arms[[paste0("random_", kk)]] <- sig_of(sample.int(nrow(map), kk))
  }
  uni <- unique(unlist(arms, use.names = FALSE))
  if (!length(uni)) return(NULL)

  ## ONE edge set and ONE truth table, shared by every arm
  edges <- ld_edges(uni, P$GTs, map[, .(marker, Chr, Pos)], P$decay_sum,
                    rho_ld = RHO_LD, dcap = DCAP)
  qtab  <- qtn_ld_table(P$GTs, map, uni, 2e6, cores = 1)

  res <- rbindlist(lapply(names(arms), function(nm) {
    mk <- arms[[nm]]
    parts <- strsplit(nm, "_(?=[0-9]+$)", perl = TRUE)[[1]]
    meth  <- parts[1]; kv <- if (length(parts) > 1) as.integer(parts[2]) else NA_integer_
    if (!length(mk)) {
      return(rbindlist(lapply(LMINS, function(L) data.table(method = meth, k = kv, l_min = L,
        n_sig = 0L, n_regions = 0L, precision = NA_real_, recall = 0, PR = 0))))
    }
    ra <- ld_regions(mk, edges)
    rbindlist(lapply(LMINS, function(L) {
      r  <- ra[lengths(ra) >= L]
      ev <- if (length(r)) evaluate_ors(r, map, qtab, th$r2min, th$dmax)
            else list(Precision = NA_real_, Recall = 0)
      data.table(method = meth, k = kv, l_min = L, n_sig = length(mk),
                 n_regions = length(r), precision = ev$Precision, recall = ev$Recall,
                 PR = (if (is.na(ev$Precision)) 0 else ev$Precision) * ev$Recall)
    }))
  }))
  cat(sprintf("    %-16s %-5s  %6d markers, %2d QTN, alpha calls %d\n",
              cell, tag, nrow(map), n_true, length(arms$alpha)))
  cbind(cell, tag, n_true, n_markers = nrow(map), res)
}

grid <- CJ(cell = CELLS, tag = c("nobgs","bgs"), sorted = FALSE)
cat(sprintf("  %d panels, engine %s, rho_ld %.2f, dcap %g, CORES=%d\n",
            nrow(grid), ENG, RHO_LD, DCAP, CORES))
run <- function(i) tryCatch(one(grid$cell[i], grid$tag[i]),
                            error = function(e) { message("  FAIL ", grid$cell[i], " ", grid$tag[i], ": ", conditionMessage(e)); NULL })
res <- if (CORES > 1) {
  mclapply(seq_len(nrow(grid)), run, mc.cores = CORES, mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(grid)), run)
}
all <- rbindlist(Filter(Negate(is.null), res), fill = TRUE)
stopifnot(nrow(all) > 0)
fn <- file.path(OUT, sprintf("or_filter_test_%s.csv", ENG))
fwrite(all, fn)
cat(sprintf("\n  written: %s (%d rows, %d panels)\n", fn, nrow(all), uniqueN(all[, .(cell,tag)])))
