## module_sticklebacks/14_poster_manhattan.R
## REVERSE-ENGINEER the poster Manhattan (R_080726/ld_w_filtering_3sp.R). Poster recipe:
##   1. C-score = summarise_stability over run_one_grid at l_min=10, r2_th=0.6 (this is
##      the POSTER's cluster-based C: consistency of >=l_min cluster membership across
##      rho x th_ldw x r2_th at alpha=0.05).
##   2. select markers with C_<method> > 0.05.
##   3. re-cluster them at r2_th=0.4 (precompute_LD_edges r2_min=0.1/max_bp=1e6 +
##      LD_igraph_components), keep n_loci >= 10.
##   4. Manhattan: -log10(<method>_q) faceted by Chr, coloured by cluster; q=0.05 line.
## Uses emx_p_GC / lfmm_P / lfmm_q (poster engine) + ld_ws_3sp_MAF01 (grid 0.25-0.99;
## poster used 0-0.95 + a permuted rho=0 control -- minor deviation, noted).
## Run: Rscript module_sticklebacks/14_poster_manhattan.R
## Output: module_sticklebacks/fig_poster_manhattan_{LFMM,EMX}.png

suppressMessages({ library(data.table); library(igraph); library(ggplot2); library(parallel) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"; CORES <- 4L

## ---------- poster functions (verbatim from ld_w_filtering_3sp.R) ----------
precompute_LD_edges <- function(GTs, map, r2_min = 0.1, max_bp = Inf, cores = 1) {
  markers <- intersect(colnames(GTs), map$marker); map_sub <- copy(map[marker %in% markers]); setkey(map_sub, Chr, Pos)
  out <- mclapply(unique(map_sub$Chr), function(ch) {
    chr_map <- map_sub[Chr == ch]; chr_markers <- chr_map$marker
    if (length(chr_markers) < 2) return(data.table(Chr = ch, marker1 = chr_markers, marker2 = chr_markers, r2 = 1, dist_bp = 0))
    gts <- as.matrix(GTs[, chr_markers, drop = FALSE]); storage.mode(gts) <- "double"
    R2 <- cor(gts, use = "pairwise.complete.obs")^2; R2[is.na(R2)] <- 0; diag(R2) <- 0
    idx <- which(R2 >= r2_min, arr.ind = TRUE); idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
    if (nrow(idx) == 0) return(data.table(Chr = ch, marker1 = chr_markers, marker2 = chr_markers, r2 = 1, dist_bp = 0))
    pos <- chr_map$Pos
    dt <- data.table(Chr = ch, marker1 = chr_markers[idx[, 1]], marker2 = chr_markers[idx[, 2]],
                     r2 = R2[idx], dist_bp = abs(pos[idx[, 1]] - pos[idx[, 2]]))
    if (is.finite(max_bp)) dt <- dt[dist_bp <= max_bp]; dt }, mc.cores = cores)
  out <- rbindlist(out, use.names = TRUE, fill = TRUE); setkey(out, Chr, marker1, marker2); out }

LD_igraph_components <- function(el, markers, r2_th = 0.8, bp_th = Inf) {
  markers <- unique(markers)
  if (length(markers) == 0) return(data.table(marker = character(), CL_id = integer(), n_loci = integer()))
  if (length(markers) == 1) return(data.table(marker = markers, CL_id = 1L, n_loci = 1L))
  edges <- el[marker1 %in% markers & marker2 %in% markers & r2 >= r2_th]
  if (is.finite(bp_th)) edges <- edges[dist_bp <= bp_th]
  if (nrow(edges) == 0) return(data.table(marker = markers, CL_id = seq_along(markers), n_loci = 1L))
  g <- graph_from_data_frame(edges[, .(from = marker1, to = marker2)], directed = FALSE, vertices = data.table(name = markers))
  comp <- components(g); clusters <- data.table(marker = names(comp$membership), CL_id = as.integer(comp$membership))
  clusters[, n_loci := comp$csize[CL_id]]; clusters }

empty_result <- function(rho, th_ldw, p_names, r2_grid, lmin_grid) {
  out <- CJ(r2_th = r2_grid, l_min = lmin_grid); out[, `:=`(th_ldw = th_ldw, rho = rho)]
  for (nm in p_names) out[, (nm) := list(list(character()))]; out[] }

run_one_grid <- function(map, el = NULL, ld_ws, rho, th_ldw, p_cols, p_names = names(p_cols),
                         alpha = 0.05, r2_grid, lmin_grid, bp_th = Inf, cores) {
  common <- intersect(map$marker, rownames(ld_ws)); map_sub <- copy(map[marker %in% common]); ld_sub <- ld_ws[map_sub$marker, , drop = FALSE]
  keep <- if (th_ldw > 0) { k <- ld_sub[, rho] > quantile(ld_sub[, rho], th_ldw, na.rm = TRUE); k[is.na(k)] <- FALSE; k } else rep(TRUE, nrow(ld_sub))
  if (!any(keep)) return(cbind(empty_result(rho, th_ldw, p_names, r2_grid, lmin_grid), n_loci = 0))
  markers_keep <- map_sub[keep, marker]; outliers <- setNames(vector("list", length(p_cols)), p_names)
  for (i in seq_along(p_cols)) { pc <- p_cols[[i]]; q <- p.adjust(unlist(map_sub[keep, ..pc]), method = "fdr"); outliers[[p_names[i]]] <- markers_keep[q < alpha] }
  if (length(unique(unlist(outliers))) == 0) return(cbind(empty_result(rho, th_ldw, p_names, r2_grid, lmin_grid), n_loci = length(which(keep))))
  out <- rbindlist(mclapply(r2_grid, function(r2_th) {
    clusters <- lapply(outliers, function(mk) LD_igraph_components(el = el, markers = mk, r2_th = r2_th, bp_th = bp_th))
    rbindlist(lapply(lmin_grid, function(l_min) { row <- data.table(r2_th = r2_th, l_min = l_min, th_ldw = th_ldw, rho = rho)
      for (nm in p_names) { tmp <- clusters[[nm]][n_loci >= l_min]; cls <- split(tmp$marker, tmp$CL_id)
        row[, (nm) := if (length(cls) > 1) list(cls) else list(list(cls))] }; row }), fill = TRUE) }, mc.cores = cores), fill = TRUE)
  out[, n_loci := length(which(keep))]; out }

get_potential_outliers <- function(map, ld_ws, th_ldw_grid, p_cols, alpha = 0.05) {
  potential <- character()
  for (rho in colnames(ld_ws)) { ld_vec <- ld_ws[, rho]
    for (th_ldw in th_ldw_grid) { keep <- ld_vec > quantile(ld_vec, th_ldw, na.rm = TRUE); keep[is.na(keep)] <- FALSE
      if (!any(keep)) next
      for (p_col in p_cols) { q <- p.adjust(unlist(map[keep, ..p_col]), method = "fdr"); potential <- c(potential, map[keep, marker][q < alpha]) } } }
  unique(potential) }

summarise_stability <- function(outliers, map, p_names) { map_C <- copy(map)
  for (nm in p_names) { C <- outliers[, table(unlist(get(nm))) / .N]
    if (length(C) > 0) { C <- data.table(C); setnames(C, c("V1", "N"), c("marker", paste0("C_", nm)))
      map_C <- C[map_C, on = "marker"]; map_C[is.na(get(paste0("C_", nm))), (paste0("C_", nm)) := 0] } }
  for (nm in p_names) if (!paste0("C_", nm) %in% names(map_C)) map_C[, (paste0("C_", nm)) := 0]   # guard: no clusters -> C=0
  map_C[] }

## ---------- data ----------
sr <- as.data.table(readRDS("/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/3sp/SNP_res_3sp.rds"))
map_3sp <- sr[, .(Chr, Pos, marker, emx_p_GC, lfmm_P, lfmm_q, emx_q)]
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs_3sp <- e$GTs_3sp; colnames(GTs_3sp) <- e$map_3sp$marker
ld_ws <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/ld_ws_3sp_MAF01.rds"); ld_ws[is.na(ld_ws)] <- 0
map_3sp <- map_3sp[marker %in% rownames(ld_ws)]; ld_ws <- ld_ws[map_3sp$marker, ]
th_ldw_grid <- 1 - 10^seq(log10(1), log10(0.01), length.out = 20)
p_cols <- c(EMX = "emx_p_GC", LFMM = "lfmm_P")
r2_grid <- seq(0.6, 0.9, by = 0.05); lmin_grid <- c(1, 5, 10, 20, 40, 80, 160)
cat(sprintf("map %d SNPs; ld_ws grid %s\n", nrow(map_3sp), paste(colnames(ld_ws), collapse=",")))

potential <- get_potential_outliers(map_3sp, ld_ws, th_ldw_grid, p_cols, alpha = 0.05)
cat(sprintf("potential outliers: %d\n", length(potential)))
el_pot <- precompute_LD_edges(GTs_3sp[, potential, drop = FALSE], map_3sp[marker %in% potential], r2_min = 0.1, max_bp = 1e6, cores = CORES)
param_grid <- CJ(th_ldw = th_ldw_grid, rho = colnames(ld_ws))
outliers_grid <- rbindlist(mclapply(seq_len(nrow(param_grid)), function(i)
  run_one_grid(map_3sp, el_pot, ld_ws, param_grid$rho[i], param_grid$th_ldw[i], p_cols, alpha = 0.05,
               r2_grid = r2_grid, lmin_grid = lmin_grid, bp_th = Inf, cores = 1), mc.cores = CORES), fill = TRUE)
map_C <- summarise_stability(outliers_grid[l_min == 10 & r2_th == 0.6], map_3sp, names(p_cols))
saveRDS(list(map_C = map_C, outliers_grid = outliers_grid), file.path(mod, "poster_cscore.rds"))
cat(sprintf("C_EMX>0.05: %d markers | C_LFMM>0.05: %d markers\n", map_C[C_EMX > 0.05, .N], map_C[C_LFMM > 0.05, .N]))

col_vector <- rep(c("#B2DF8A","#FFD92F","firebrick","#33A02C","#7FC97F","#CAB2D6","#FB8072","#E6AB02","steelblue",
  "#FB9A99","#1B9E77","#BC80BD","#E31A1C","#7570B3","#A6761D","#A6CEE3","salmon","forestgreen","#BF5B17"), 20)
manhattan <- function(Ccol, qcol, tag) {
  outl <- map_C[get(Ccol) > 0.05, marker]; if (!length(outl)) { cat("no outliers for", tag, "\n"); return(invisible()) }
  el <- precompute_LD_edges(GTs_3sp[, outl, drop = FALSE], map_3sp[marker %in% outl], r2_min = 0.1, max_bp = 1e6, cores = CORES)
  cls <- LD_igraph_components(el, outl, r2_th = 0.4, bp_th = Inf)
  mh <- cls[n_loci >= 10, .(marker, CL_id)][map_3sp, on = "marker"]
  mh[, CL_id := as.character(CL_id)][is.na(CL_id), CL_id := "ns"]
  cat(sprintf("[%s] %d >=10-SNP clusters across %d chr\n", tag, mh[CL_id!="ns", uniqueN(CL_id)], mh[CL_id!="ns", uniqueN(Chr)]))
  p <- ggplot(mh, aes(Pos, -log10(get(qcol)))) +
    geom_point(data = mh[CL_id == "ns"], color = "grey35", size = 0.5, alpha = 0.3) +
    geom_point(data = mh[CL_id != "ns"], aes(color = CL_id), size = 1.8) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.6) +
    facet_wrap(~ factor(Chr, levels = paste0("Chr", c(1:18,20,21))), nrow = 1, scales = "free_x") +
    scale_color_manual(values = rep(col_vector, 10)) +
    labs(x = "Genomic position", y = bquote(-log[10](italic(q))~"("*.(tag)*")"),
         title = sprintf("Poster Manhattan (reverse-engineered): %s, C_%s>0.05, >=10-SNP clusters (r2=0.4)", tag, tag)) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", panel.grid = element_blank(), strip.text = element_text(face = "bold", size = 7),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(0.05, "lines"))
  ggsave(file.path(mod, sprintf("fig_poster_manhattan_%s.png", tag)), p, width = 18, height = 4.5, dpi = 200)
  cat("wrote fig_poster_manhattan_", tag, ".png\n", sep = "") }
manhattan("C_LFMM", "lfmm_q", "LFMM")
manhattan("C_EMX",  "emx_q",  "EMX")
