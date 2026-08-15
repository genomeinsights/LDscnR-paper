## module_sim/03_ldw_qtn_manhattan.R
## ld_w Manhattan coloured by max r^2 with a true-positive QTN. This shows WHY
## the caller fails in the structure-limited (c2) regime: if high ld_w tracked
## QTNs the peaks would be red; blue peaks are structure LD with no QTN linkage.
## Quantifies the fraction of high-ld_w candidates that are actually QTN-linked.
## Reads the meta_* written by 01_score_pooled.R (map carries max_LD_with_QTN,
## ld_w_095, true_pos_QTN).
## Run from LDscnR-paper/:  Rscript module_sim/03_ldw_qtn_manhattan.R [V c env]
## Output (git-ignored): module_sim/ldw_qtn_manhattan_V{V}_c{c}_env{env}.png

suppressMessages({ library(data.table); library(ggplot2) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sim"
dir_data <- file.path(mod, "data"); dir_fig <- file.path(mod, "figures")
for (d in c(dir_data, dir_fig)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "1"
cc  <- if (length(a) >= 2) a[2] else "2"
env <- if (length(a) >= 3) a[3] else "3"
m   <- as.data.table(readRDS(file.path(dir_data, sprintf("meta_V%s_c%s_env%s.rds", V, cc, env)))$map)

## genome coordinate transform (moduleB style)
m[, chr_num := as.integer(factor(Chr, levels = unique(Chr)))]
clen <- m[, .(len = max(Pos)), by = .(Chr, chr_num)][order(chr_num)]
spg  <- 0.005 * sum(clen$len); clen[, start := c(0, head(cumsum(len + spg), -1))]
m    <- clen[, .(Chr, start)][m, on = "Chr"][, x := Pos + start]
cmid <- clen[, .(Chr, mid = start + len / 2)]
shade <- clen[chr_num %% 2 == 0, .(xmin = start, xmax = start + len)]
qs   <- quantile(m$ld_w_095, 0.855)          # illustrative candidate threshold
qtn  <- m[true_pos_QTN == TRUE]
setorder(m, max_LD_with_QTN)                 # plot QTN-linked points on top

p <- ggplot() +
  geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93") +
  geom_point(data = m, aes(x, ld_w_095, color = max_LD_with_QTN), size = 0.35, alpha = 0.7) +
  geom_point(data = qtn, aes(x, ld_w_095), shape = 25, fill = "black", color = "black", size = 1.8) +
  geom_hline(yintercept = qs, linetype = 2, color = "red", linewidth = 0.4) +
  scale_color_gradientn(colors = c("#3A3A98", "#5DC8CD", "#F2E205", "#F29724", "#D62828"),
                        limits = c(0, 1), name = expression(max ~ r^2 ~ with ~ QTN)) +
  scale_x_continuous(breaks = cmid$mid, labels = gsub("R([0-9]+)_Chr", "\\1.", cmid$Chr),
                     expand = c(0.01, 0.01)) +
  labs(x = "pooled chromosome (Ri.Chrj; 10 QTN-bearing + 10 neutral)",
       y = expression(ld[w] ~ (rho == 0.95)),
       title = sprintf("V%s_c%s: ld_w Manhattan coloured by max r^2 with a true-positive QTN", V, cc),
       subtitle = "down-triangle = true-positive QTN. Red peaks = QTN linkage; blue peaks = structure LD, not QTN linkage") +
  theme_bw(base_size = 9) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 6), legend.position = "top")
ggsave(file.path(dir_fig, sprintf("ldw_qtn_manhattan_V%s_c%s_env%s.png", V, cc, env)), p,
       width = 16, height = 5.5, dpi = 150)

## quantify: among high-ld_w candidates, how many are actually QTN-linked?
cand <- m[ld_w_095 >= qs]
cat(sprintf("candidates (ld_w>=%.3f): n=%d\n", qs, nrow(cand)))
cat(sprintf("  max_LD_with_QTN among candidates: median=%.3f  frac r2>0.5=%.3f  frac r2>0.8=%.3f\n",
    median(cand$max_LD_with_QTN), mean(cand$max_LD_with_QTN > 0.5), mean(cand$max_LD_with_QTN > 0.8)))
cat(sprintf("  cor(ld_w, max_LD_with_QTN) genome-wide = %.3f\n",
    cor(m$ld_w_095, m$max_LD_with_QTN, method = "spearman")))
cat("wrote", sprintf("ldw_qtn_manhattan_V%s_c%s_env%s.png", V, cc, env), "\n")
