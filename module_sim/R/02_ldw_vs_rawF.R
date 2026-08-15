## module_sim/02_ldw_vs_rawF.R
## THE mechanism figure: raw association F vs ld_w, per chromosome, for a pooled
## 20-chromosome genome (top row = QTN-bearing, bottom = neutral). Use RAW F, not
## a percentile -- percentiles normalise each method's distribution and HIDE the
## effect. EMMAX F is well-calibrated (median ~0.42 = null F(1,158)) and flat in
## ld_w on every chromosome; LFMM F is baseline-inflated AND rises with ld_w on
## every chromosome INCLUDING the neutral ones -> the ld_w filter concentrates
## LFMM's spurious LD-correlated inflation into high-LD regions -> neutral FP.
## Run from LDscnR-paper/:  Rscript module_sim/02_ldw_vs_rawF.R [V c env]
## Output (git-ignored): module_sim/ldw_vs_rawF_V{V}_c{c}_env{env}.png

suppressMessages({ library(data.table); library(ggplot2) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sim"
dir_fig <- file.path(mod, "figures"); if (!dir.exists(dir_fig)) dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
SIM_DATA <- Sys.getenv("LDSCNR_SIM_DATA",
                       unset = "/Volumes/Nemo/Nemo_sim/parsed_sim_data2")
if (!dir.exists(SIM_DATA)) SIM_DATA <- "parsed_sim_data"
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"

files <- list.files(SIM_DATA, full.names = TRUE,
                    pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
reps  <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
pm <- rbindlist(lapply(seq_along(files), function(i) {
  d <- readRDS(files[i]); m <- as.data.table(d$map)
  lwc <- if ("rho_0.95" %in% colnames(d$ld_ws)) "rho_0.95" else "0.95"
  data.table(Chr = paste0("R", i, "_", m$Chr), ld_w = d$ld_ws[, lwc],
             EMMAX = m$emx_F, LFMM = m$lfmm_F)
}))
pm[, ldw_q := ecdf(ld_w)(ld_w)]
long <- melt(pm[, .(Chr, ldw_q, EMMAX, LFMM)], id.vars = c("Chr", "ldw_q"),
             variable.name = "method", value.name = "F")
setorder(long, method, Chr, ldw_q)
long[, Froll := frollapply(F, 201, median, align = "center"), by = .(method, Chr)]
ord <- c(paste0("R", 1:10, "_Chr1"), paste0("R", 1:10, "_Chr2"))
long[, Chr := factor(Chr, levels = ord, labels = c(paste0("Q", 1:10), paste0("N", 1:10)))]

p <- ggplot(long, aes(ldw_q, colour = method)) +
  geom_point(aes(y = F), size = 0.15, alpha = 0.05) +
  geom_line(aes(y = Froll), linewidth = 0.4, na.rm = TRUE) +
  scale_colour_manual(values = c(EMMAX = "#F21A00", LFMM = "#3B9AB2")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  facet_wrap(~Chr, nrow = 2) + coord_cartesian(ylim = c(0, 8)) +
  labs(x = "ld_w quantile (genome-wide)",
       y = "raw F (points = raw, line = windowed median w=201; y clipped at 8)",
       title = sprintf("ld_w vs raw F per chromosome - V%s_c%s_env%s (Q = QTN-bearing top, N = neutral bottom)", V, cc, env)) +
  theme_bw(base_size = 9) + theme(legend.position = "top", strip.background = element_blank())
ggsave(file.path(dir_fig, sprintf("ldw_vs_rawF_V%s_c%s_env%s.png", V, cc, env)), p,
       width = 15, height = 6, dpi = 140)
cat("wrote", sprintf("ldw_vs_rawF_V%s_c%s_env%s.png", V, cc, env), "\n")
