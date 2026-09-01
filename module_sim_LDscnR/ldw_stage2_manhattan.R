## =============================================================================
## ld_w manhattan coloured by STAGE-2 cluster, for one V0.5_c1 bundle.
##
## Stage-2 membership is NOT persisted in the bundles -- the parse stored only
## stage1 and the pruned marker vector -- so it is recomputed here with the same
## arguments as clusters_as_test_units.R, including the check that the
## recomputation reproduces the stored GRM marker set exactly.
##
## Colours are rotated along the genome so adjacent clusters differ; singleton
## clusters (n_loci == 1) are grey.
##
## Env: SIM_DATA, CELL, TAG, ENV, FILE, OUT, DCAP, Y (ld_w | q), ENGINE (emx | lfmm)
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(ggplot2)})

SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELL <- Sys.getenv("CELL", "V0.5_c1")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/figures")
DCAP <- as.numeric(Sys.getenv("DCAP", "1e5"))
YVAR <- Sys.getenv("Y", "ld_w")            ## ld_w or q
ENG  <- Sys.getenv("ENGINE", "emx")        ## emx or lfmm
PICK <- Sys.getenv("TAG", "")              ## non-empty => use TAG/ENV/FILE directly
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- build one bundle's stage-2 partition + detectability ------------------
build <- function(tag, env, i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, tag, i, CELL, env)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  m <- flag_true_qtns(as.data.table(x$map))
  pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
          LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
          score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
          compute_unflagged_eMLG = FALSE, cores = 1)
  if (!identical(sort(pr$pruned), sort(x$grm_markers)))
    stop(sprintf("stage-2 recomputation does not reproduce grm_markers for %s", basename(f)))
  g  <- as.data.table(pr$groups)
  ms <- rbindlist(lapply(seq_len(nrow(g)), function(k)
          data.table(marker = g$members[[k]], CL_id = g$group_id[k], n_loci = g$n_loci[k])))
  m  <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL_id)]
  m[, ld_w := x$ld_ws[marker, "rho_0.95"]]
  ## BH across BOTH chromosomes pooled, matching clusters_as_test_units.R
  pcol <- paste0(ENG, "_p")
  if (!pcol %in% names(m)) stop("no column ", pcol, " in map")
  m[, q := p.adjust(get(pcol), "BH")]

  ## detectable QTN: at least one retained marker at r2 >= r2min within dmax
  th  <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
                          rho_r2 = 0.75, rho_d = 0.95, dmax_cap = DCAP)
  qtn <- m[true_pos_QTN %in% TRUE]
  det <- rbindlist(lapply(seq_len(nrow(qtn)), function(j) {
    ch   <- as.character(qtn$Chr[j])
    near <- m[as.character(Chr) == ch & abs(Pos - qtn$Pos[j]) < th$dmax]
    if (!nrow(near)) return(NULL)
    r2 <- suppressWarnings(cor(x$GTs[, qtn$marker[j]], x$GTs[, near$marker],
                               use = "pairwise.complete.obs")^2)
    ok <- which(is.finite(r2) & r2 >= th$r2min)
    if (!length(ok)) return(NULL)
    data.table(marker = qtn$marker[j], Chr = qtn$Chr[j], Pos = qtn$Pos[j],
               n_tag = length(ok), best_r2 = max(r2[ok], na.rm = TRUE))
  }))
  list(m = m, det = det, n_qtn = nrow(qtn), file = i, tag = tag, env = env,
       r2min = th$r2min, dmax = th$dmax)
}

## ---- scan a few bundles, take the one with most detectable QTN -------------
cand <- list()
grid <- if (nzchar(PICK)) {
  list(list(PICK, as.integer(Sys.getenv("ENV", "1")), as.integer(Sys.getenv("FILE", "1"))))
} else {
  unlist(lapply(c("nobgs", "bgs"), function(tg)
    unlist(lapply(1:3, function(e) lapply(1:3, function(i) list(tg, e, i))),
           recursive = FALSE)), recursive = FALSE)
}
for (gg in grid) { tag <- gg[[1]]; env <- gg[[2]]; i <- gg[[3]]
  b <- tryCatch(build(tag, env, i), error = function(e) {message("  skip: ", conditionMessage(e)); NULL})
  if (is.null(b) || is.null(b$det) || !nrow(b$det)) next
  cat(sprintf("  %-5s env%-2d chr%-2d  QTN %d, detectable %d, markers %d, clusters %d (singleton %.0f%%)\n",
      tag, env, i, b$n_qtn, nrow(b$det), nrow(b$m), uniqueN(b$m$CL_id),
      100 * mean(b$m[!duplicated(CL_id)]$n_loci == 1)))
  cand[[length(cand) + 1]] <- b
}
stopifnot(length(cand) > 0)
best <- cand[[which.max(sapply(cand, function(z) nrow(z$det)))]]
cat(sprintf("\n  CHOSEN: %s env%d chr%d -- %d detectable QTN (r2min %.3f, dmax %.0f kb)\n",
    best$tag, best$env, best$file, nrow(best$det), best$r2min, best$dmax / 1000))

## ---- colour rotation: order clusters along the genome, cycle a palette -----
m <- copy(best$m)
m[, singleton := n_loci == 1]
ord <- m[singleton == FALSE, .(start = min(Pos)), by = .(Chr, CL_id)][order(Chr, start)]
pal <- c("#1F6F8B", "#C1622F", "#2E7156", "#8E5AA8", "#B0392B", "#3D7EA6", "#7A6A1F", "#456A8C")
ord[, col := pal[(seq_len(.N) - 1L) %% length(pal) + 1L]]
m <- merge(m, ord[, .(CL_id, col)], by = "CL_id", all.x = TRUE)
m[singleton == TRUE, col := "grey72"]

det <- best$det
cat(sprintf("  non-singleton clusters: %d ; singleton: %d\n",
    nrow(ord), uniqueN(m[singleton == TRUE]$CL_id)))

m[, yval := if (YVAR == "q") -log10(pmax(q, .Machine$double.xmin)) else ld_w]
ylab <- if (YVAR == "q") {
  bquote(-log[10](q)~"("*.(ENG)*", BH pooled over both chromosomes)")
} else {
  expression(ld[w]~"("*rho*" = 0.95)")
}
sig  <- if (YVAR == "q") -log10(0.05) else NA_real_

p <- ggplot(m[order(singleton, decreasing = TRUE)],
            aes(Pos / 1e6, yval, colour = col)) +
  geom_vline(data = det, aes(xintercept = Pos / 1e6),
             colour = "grey35", linetype = "22", linewidth = .4) +
  {if (YVAR == "q") geom_hline(yintercept = sig, colour = "grey35", linewidth = .35) else NULL} +
  geom_point(size = .55, alpha = .85) +
  scale_colour_identity() +
  facet_wrap(~ Chr, ncol = 1, scales = "free_x") +
  labs(x = "Position (Mb)", y = ylab,
       title = sprintf("%s %s env%d chr%d - %s by stage-2 cluster",
                       CELL, best$tag, best$env, best$file,
                       if (YVAR == "q") sprintf("-log10(q) [%s]", ENG) else "ld_w"),
       subtitle = sprintf(paste0("%d markers in %d stage-2 clusters; singletons grey. ",
                          "Dashed lines: %d detectable QTN (r2 >= %.2f within %.0f kb)."),
                          nrow(m), uniqueN(m$CL_id), nrow(det), best$r2min, best$dmax / 1000)) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", plot.subtitle = element_text(colour = "grey30", size = 8))

ytag <- if (YVAR == "q") paste0("q", ENG) else "ldw"
fn <- file.path(OUT, sprintf("%s_stage2clusters_%s_%s_env%d_chr%d.png",
                             ytag, CELL, best$tag, best$env, best$file))
ggsave(fn, p, width = 10, height = 6.2, dpi = 190)
cat(sprintf("\n  written: %s\n", fn))
