## =============================================================================
## filter_then_test.R -- is ld_w useful as a phenotype-blind SELECTION device?
##
## This tests a DIFFERENT claim from everything else in this module. The C-score
## work asked whether ld_w improves a RANKING; it does not. This asks whether
## ld_w improves a SELECTION made BEFORE the phenotype is seen: choose k units on
## genotype information alone, run BH inside that set only, and see whether more
## true discoveries survive than a genome-wide BH makes.
##
## The mechanism is ordinary independent filtering: fewer tests means a less
## severe multiplicity correction. That is arithmetic and would work for ANY
## filter, so the informative comparison is not against the genome-wide scan --
## it is against other filters of the SAME SIZE:
##
##   random   -- is enrichment needed at all, or does any size reduction do it?
##   MAF      -- the sharp control. ld_w correlates with MAF, and MAF drives
##               power, so a MAF filter also concentrates power. If ld_w does not
##               beat it, "ld_w works" reduces to "high-MAF markers work".
##   rec_rate -- ld_w keys on low recombination; does raw recombination suffice?
##
## TRUTH is bp distance to the nearest DRIVING QTN with the causal markers
## themselves EXCLUDED. That truth shares no machinery with ld_w. Scoring an LD
## statistic against an LD-defined truth (r2-tagging) flatters it by construction
## -- that error is documented in this project and is not repeated here.
##
## WHAT THIS SCRIPT CANNOT DO. It measures POWER, not FDR VALIDITY. ld_w is
## phenotype-blind but NOT independent of the test statistic (both respond to
## MAF), so filtered BH is not guaranteed to control FDR. That requires
## permuting the phenotype and recomputing the filter per draw, which needs
## genotypes; these tracks have none. Treat a power gain here as necessary, not
## sufficient.
##
## Env: TRACKS, OUT, KS (selection sizes), WINDOWS (kb), ALPHA, ENGINE
## =============================================================================
suppressMessages({library(data.table)})

TR   <- Sys.getenv("TRACKS", "/Volumes/Nemo/Nemo_sim/pilot_v2_v05c1/tracks")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
KS   <- as.integer(strsplit(Sys.getenv("KS", "250,500,1000,2000,5000"), ",")[[1]])
WIN  <- as.numeric(strsplit(Sys.getenv("WINDOWS", "10,50,100"), ",")[[1]]) * 1000
ALPHA<- as.numeric(Sys.getenv("ALPHA", "0.05"))
ENG  <- Sys.getenv("ENGINE", "emx")
SEED <- as.integer(Sys.getenv("SEED", "20260901"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

pcol <- paste0(ENG, "_p")

## score one selection: BH within the selected set only
score_sel <- function(m, sel, w) {
  s <- m[sel]
  if (!nrow(s)) return(data.table(n_sel = 0L, n_sig = 0L, n_tp = 0L, n_fp = 0L,
                                  precision = NA_real_, qtn_hit = 0L))
  q <- p.adjust(s[[pcol]], "BH")
  sig <- which(q < ALPHA)
  if (!length(sig)) return(data.table(n_sel = nrow(s), n_sig = 0L, n_tp = 0L, n_fp = 0L,
                                      precision = NA_real_, qtn_hit = 0L))
  hit <- s[sig]
  tp  <- hit$d_qtn < w                       # d_qtn: bp to nearest driving QTN
  ## which distinct driving QTN were recovered
  nq <- uniqueN(hit$near_qtn[tp])
  data.table(n_sel = nrow(s), n_sig = length(sig), n_tp = sum(tp), n_fp = sum(!tp),
             precision = mean(tp), qtn_hit = nq)
}

one <- function(f) {
  t <- readRDS(f)
  m <- as.data.table(t$map)
  if (!pcol %in% names(m)) return(NULL)
  m <- m[is.finite(get(pcol))]
  if (!nrow(m)) return(NULL)

  ## driving QTN, and the distance from every marker to the nearest one.
  ## The causal markers themselves are REMOVED from the tested set -- otherwise
  ## the test rewards finding the variant rather than its neighbourhood.
  drv <- m[true_QTN %in% TRUE & MAF > 0.1 & p_Va > 0.05]
  if (!nrow(drv)) return(NULL)
  m[, `:=`(d_qtn = Inf, near_qtn = NA_character_)]
  for (ch in unique(m$Chr)) {
    dd <- drv[Chr == ch]; if (!nrow(dd)) next
    ii <- which(m$Chr == ch)
    D  <- abs(outer(m$Pos[ii], dd$Pos, "-"))
    j  <- max.col(-D)
    m[ii, `:=`(d_qtn = D[cbind(seq_along(ii), j)],
               near_qtn = paste0(ch, "_", dd$Pos[j]))]
  }
  m <- m[!(true_QTN %in% TRUE)]              # exclude causal variants
  if (!nrow(m)) return(NULL)

  meta <- data.table(variant = t$variant, arm = t$arm, set = t$set,
                     file = basename(f), n_mk = nrow(m), n_drv = nrow(drv))

  res <- rbindlist(lapply(WIN, function(w) {
    base <- cbind(meta, method = "genome_wide", k = nrow(m), window_kb = w/1000,
                  score_sel(m, seq_len(nrow(m)), w))
    per_k <- rbindlist(lapply(KS, function(k) {
      if (k >= nrow(m)) return(NULL)
      o_ld  <- head(order(-m$ld_w_095), k)
      o_maf <- head(order(-m$MAF), k)
      o_rec <- if ("rec_rate" %in% names(m)) head(order(m$rec_rate), k) else NULL   # LOW recomb
      o_rnd <- sample.int(nrow(m), k)
      sels <- list(ld_w = o_ld, MAF = o_maf, random = o_rnd)
      if (!is.null(o_rec)) sels$rec_low <- o_rec
      rbindlist(lapply(names(sels), function(nm)
        cbind(meta, method = nm, k = k, window_kb = w/1000, score_sel(m, sels[[nm]], w))))
    }))
    rbind(base, per_k, fill = TRUE)
  }))
  res
}

fs <- list.files(TR, pattern = "\\.rds$", full.names = TRUE)
cat(sprintf("  %d track files, engine %s, alpha %.2f\n", length(fs), ENG, ALPHA))
all <- rbindlist(lapply(seq_along(fs), function(i) {
  r <- tryCatch(one(fs[i]), error = function(e) {message("  skip ", basename(fs[i]), ": ", conditionMessage(e)); NULL})
  if (!is.null(r) && i %% 16 == 0) cat(sprintf("    %d/%d\n", i, length(fs)))
  r
}), fill = TRUE)
stopifnot(nrow(all) > 0)
fwrite(all, file.path(OUT, sprintf("filter_then_test_%s.csv", ENG)))

## ---- PAIRED reporting. Never a difference of means. ------------------------
sgn <- function(x, lab) {
  x <- x[is.finite(x)]; nz <- x[x != 0]
  if (!length(nz)) { cat(sprintf("    %-30s n=%3d  all ties\n", lab, length(x))); return(invisible(NULL)) }
  cat(sprintf("    %-30s n=%3d  median_d=%+6.1f  A>B %2d/%2d (%3.0f%%)  sign p=%.3g\n",
      lab, length(x), median(x), sum(nz > 0), length(nz), 100*sum(nz>0)/length(nz),
      binom.test(sum(nz>0), length(nz))$p.value))
}
key <- c("file","window_kb")
for (w in sort(unique(all$window_kb))) {
  cat(sprintf("\n=== window %.0f kb : TRUE DISCOVERIES (n_tp), paired within run ===\n", w))
  gw <- all[method == "genome_wide" & window_kb == w, c(key,"n_sig","n_tp","qtn_hit"), with = FALSE]
  setnames(gw, c("n_sig","n_tp","qtn_hit"), c("gw_sig","gw_tp","gw_qtn"))
  cat(sprintf("  genome-wide baseline: median %d significant, %d true, %d QTN recovered\n",
      median(gw$gw_sig), median(gw$gw_tp), median(gw$gw_qtn)))
  for (kk in KS) {
    sub <- all[window_kb == w & k == kk & method != "genome_wide"]
    if (!nrow(sub)) next
    mg <- merge(sub, gw, by = key)
    cat(sprintf("  k = %d\n", kk))
    for (mth in c("ld_w","MAF","rec_low","random")) {
      z <- mg[method == mth]
      if (nrow(z)) sgn(z$n_tp - z$gw_tp, sprintf("%s vs genome-wide", mth))
    }
    ldw <- mg[method == "ld_w"]
    for (mth in c("MAF","rec_low","random")) {
      z <- merge(ldw[, c(key,"n_tp"), with=FALSE], mg[method==mth, c(key,"n_tp"), with=FALSE],
                 by = key, suffixes = c("_ldw","_ctl"))
      if (nrow(z)) sgn(z$n_tp_ldw - z$n_tp_ctl, sprintf("ld_w vs %s (SAME SIZE)", mth))
    }
  }
}
cat(sprintf("\n  written: %s\n", file.path(OUT, sprintf("filter_then_test_%s.csv", ENG))))
