## =====================================================================
## module_sticklebacks_LDscnR / manhattan_emp_pvals.R
##
## Genome-wide C-score Manhattan for 3sp with the outlier regions coloured by
## their LOCATION-MATCHED empirical p-value from a structure null (default the
## among-population ecotype permutation). Non-region markers are grey; region
## markers are coloured by -log10(emp_p), so the peaks that the structure-
## preserving null seldom rebuilds at their locus (the robust, ecotype-specific
## signals) light up and the structure-prone peaks stay dark.
##
## Same reference plotting style as manhattan_regions.R (single-row chromosome
## facets, grey background, no x-axis, bold strips). Recomputes regions + emp_p
## in-script (fast) so it is self-contained.
##
## Run from the LDscnR-paper root:
##   Rscript module_sticklebacks_LDscnR/manhattan_emp_pvals.R [tau_C] [l_min] [rho_ld] [null.rds]
##   defaults: 0.05 3 0.60 null_popperm_3sp.rds
## =====================================================================

suppressMessages({ library(data.table); library(LDscnR); library(ggplot2) })

a      <- commandArgs(trailingOnly = TRUE)
TAU    <- if (length(a) >= 1) as.numeric(a[1]) else 0.05
LMIN   <- if (length(a) >= 2) as.integer(a[2]) else 3L
RHO_LD <- if (length(a) >= 3) as.numeric(a[3]) else 0.60
NULLF  <- if (length(a) >= 4) a[4] else "module_sticklebacks_LDscnR/results/null_popperm_3sp.rds"
DCAP   <- 5e5
BUNDLE <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
OUTFIG <- "module_sticklebacks_LDscnR/figures"; if (!dir.exists(OUTFIG)) dir.create(OUTFIG, recursive = TRUE)

## ---- 1. data, null, regions + empirical p ----------------------------
d    <- readRDS(BUNDLE); map <- as.data.table(d$map)
null <- readRDS(NULLF); BASIS <- if (!is.null(null$basis)) null$basis else "null"; B <- length(null$C_surr)
edges <- ld_edges(null$universe, d$GTs, map[, .(marker, Chr, Pos)],
                  as.data.table(d$LD_decay$decay_sum), rho_ld = RHO_LD, dcap = DCAP)
mpos <- stats::setNames(map$Pos, map$marker); mchr <- stats::setNames(as.character(map$Chr), map$marker)

reglist <- function(Cvec) { mk <- names(Cvec)[Cvec >= TAU]
  if (!length(mk)) return(list())
  r <- ld_regions(mk, edges); r[lengths(r) >= LMIN] }
score  <- function(Cvec, r) sum(Cvec[r])
span   <- function(r) list(Chr = unname(mchr[r[1]]), lo = min(mpos[r]), hi = max(mpos[r]))

obs_regs <- reglist(null$C_obs)
stopifnot(length(obs_regs) > 0)
null_regs <- lapply(null$C_surr, function(Cs) { rr <- reglist(Cs)
  if (!length(rr)) return(data.table(Chr = character(), lo = numeric(), hi = numeric(), s = numeric()))
  rbindlist(lapply(rr, function(r) { sp <- span(r); data.table(Chr = sp$Chr, lo = sp$lo, hi = sp$hi, s = score(Cs, r)) })) })

emp_p <- vapply(obs_regs, function(r) { sp <- span(r); s_obs <- score(null$C_obs, r)
  best <- vapply(null_regs, function(nr) { if (!nrow(nr)) return(0)
    h <- nr[Chr == sp$Chr & lo <= sp$hi & hi >= sp$lo]; if (!nrow(h)) 0 else max(h$s) }, numeric(1))
  (1 + sum(best >= s_obs)) / (1 + B) }, numeric(1))
cat(sprintf("[1] %s null B=%d ; %d regions ; emp_p<0.05: %d\n", BASIS, B, length(obs_regs), sum(emp_p < 0.05)))

## ---- 2. per-marker frame: region id, emp_p; ns = grey ----------------
mh <- copy(map[, .(marker, Chr, Pos, C = null$C_obs[map$marker])]); mh[is.na(C), C := 0]
mh[, `:=`(reg = NA_integer_, emp_p = NA_real_)]
for (i in seq_along(obs_regs)) mh[marker %in% obs_regs[[i]], `:=`(reg = i, emp_p = emp_p[i])]
mh[, nlp := -log10(emp_p)]
chr_lev <- paste0("Chr", sort(as.integer(gsub("Chr", "", unique(map$Chr)))))
mh[, Chr := factor(Chr, levels = chr_lev)]
## region label anchors (top marker of each region), significant ones only
lab <- mh[!is.na(reg), .SD[which.max(C)], by = reg][emp_p < 0.05]
lab[, txt := sprintf("%s p=%.3f", sub("Chr", "chr", Chr), emp_p)]

## ---- 3. plot ---------------------------------------------------------
g <- ggplot(mh, aes(Pos, C)) +
  geom_point(data = mh[is.na(reg)], color = "grey78", size = 0.4, alpha = 0.5) +
  geom_point(data = mh[!is.na(reg)], aes(color = nlp), size = 1.7) +
  geom_text(data = lab, aes(label = txt), size = 2.4, vjust = -0.6, hjust = 0.5, color = "grey15") +
  facet_wrap(~ Chr, nrow = 1, scales = "free_x") +
  scale_color_viridis_c(option = "inferno", end = 0.9, direction = 1,
    name = expression(-log[10]~italic(p)[emp]),
    breaks = -log10(c(0.5, 0.1, 0.05, 0.02)), labels = c("0.5", "0.1", "0.05", "0.02")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Genomic position",
       y = "C-score",
       title = sprintf("3sp outlier regions coloured by %s empirical p (tau=%.2f, l_min=%d, rho_ld=%.2f, B=%d)",
                       BASIS, TAU, LMIN, RHO_LD, B)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0.05, "lines"),
        legend.position = "right")
outfn <- file.path(OUTFIG, sprintf("manhattan_emppval_%s_tau%.2f_lmin%d_rho%.2f.png", BASIS, TAU, LMIN, RHO_LD))
ggsave(outfn, g, width = 18, height = 4.6, dpi = 180)
cat(sprintf("[3] wrote %s\n", outfn))
