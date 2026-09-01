## The precision-recall trade-off, and the fact that the metric is a claim.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
R <- "module_sim_LDscnR/results/filter_then_test"; OUT <- "module_sim_LDscnR/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
d <- fread(file.path(R, "snp_vs_cluster_dedup_allpanels.csv"))
lab <- c(snp = "single SNP", rep = "representative", simes = "Simes", emlg = "eMLG")
d <- d[analysis %in% names(lab)]; d[, an := factor(lab[analysis], levels = lab)]
d[, `:=`(P = dedup_prec, R = dedup_rec)]
pal <- c("single SNP" = "#B0392B", "representative" = "#7A6A1F",
         "Simes" = "#1F6F8B", "eMLG" = "#2E7156")
fb <- function(P, R, b) { v <- (1+b^2)*P*R/(b^2*P + R); v[!is.finite(v)] <- 0; v }
BET <- c(0.25, 0.35, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4)

## A -- precision vs recall, per arm
a <- d[, .(P = median(P, na.rm = TRUE), R = median(R)), by = .(an, tag)]
pA <- ggplot(a, aes(R, P, colour = an, shape = tag)) +
  geom_point(size = 3.2) +
  geom_line(aes(group = an), linewidth = .5, alpha = .5) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_shape_manual(values = c(nobgs = 16, bgs = 17), name = NULL) +
  labs(x = "recall", y = "precision",
       title = "A  Clustering trades recall for precision",
       subtitle = "Medians over 40 panels per arm. Circles = nobgs, triangles = bgs. Lines join the same analysis across arms.") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.subtitle = element_text(colour = "grey30", size = 7.4))

## B -- F_beta curves
cur <- rbindlist(lapply(BET, function(b) { d[, f := fb(P, R, b)]
  d[, .(beta = b, F = median(f, na.rm = TRUE)), by = an] }))
pB <- ggplot(cur, aes(beta, F, colour = an)) +
  geom_line(linewidth = .8) + geom_point(size = 1.5) +
  scale_colour_manual(values = pal, name = NULL) + scale_x_log10(breaks = BET) +
  labs(x = expression(beta~"  (< 1 weights precision, > 1 weights recall; log scale)"),
       y = expression("median "*F[beta]),
       title = expression("B  Which analysis is better depends entirely on "*beta),
       subtitle = "F_beta = (1+b^2)PR/(b^2 P + R). The curves cross: no single number settles the comparison.") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.subtitle = element_text(colour = "grey30", size = 7.4))

## C -- paired win rate against beta, with the crossover band
wr <- rbindlist(lapply(BET, function(b) { d[, f := fb(P, R, b)]
  v <- dcast(d, cell + tag + env ~ an, value.var = "f")
  x <- v$Simes - v$`single SNP`; x <- x[is.finite(x)]; nz <- x[x != 0]
  data.table(beta = b, win = 100*sum(nz>0)/length(nz),
             p = binom.test(sum(nz>0), length(nz))$p.value) }))
pC <- ggplot(wr, aes(beta, win)) +
  annotate("rect", xmin = 1.5, xmax = 2.0, ymin = -Inf, ymax = Inf, fill = "grey88") +
  geom_hline(yintercept = 50, linetype = "22", colour = "grey30") +
  geom_line(colour = "#1F6F8B", linewidth = .8) +
  geom_point(aes(shape = p < 0.05), colour = "#1F6F8B", size = 2.2) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("n.s.", "p < 0.05"), name = NULL) +
  scale_x_log10(breaks = BET) + ylim(0, 100) +
  labs(x = expression(beta), y = "% of 80 panels where Simes beats the single-SNP scan",
       title = expression("C  The crossover sits between "*beta*" = 1.5 and 2"),
       subtitle = "Shaded: the region where the two are indistinguishable. Left of it clustering wins, right of it the single-SNP scan does.") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.subtitle = element_text(colour = "grey30", size = 7.4))

ggsave(file.path(OUT, "fbeta_tradeoff.png"), pA / pB / pC, width = 9, height = 10, dpi = 190)
cat(sprintf("  written: %s\n", file.path(OUT, "fbeta_tradeoff.png")))
