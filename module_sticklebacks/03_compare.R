## module_sticklebacks/03_compare.R
## Compare the three cluster-level analyses (EMMAX+perm, EMMAX+bg, LFMM+bg):
## RMSC curves, background rates, cluster concordance, top regions.
## Run from LDscnR-paper/:  Rscript module_sticklebacks/03_compare.R

suppressMessages({ library(data.table); library(ggplot2) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
R <- readRDS(file.path(mod, "results_three_nulls.rds"))
nm <- c(emx_perm = "EMMAX+perm", emx_bg = "EMMAX+bg", lfmm_bg = "LFMM+bg")

## ---- per-method headline + diagnostics ------------------------------------
cat("=== headline ===\n")
for (k in names(nm)) {
  d <- R[[k]]$diagnostics
  cat(sprintf("%-9s q*=%.3f (ld_w>=%.3f) | cand=%d clusters=%d sig=%d | bg_rate=%.4f | no_elbow=%s all_chr=%s\n",
              nm[k], R[[k]]$ld_w_threshold, R[[k]]$ld_w_value, nrow(R[[k]]$candidates),
              nrow(R[[k]]$clusters), sum(R[[k]]$clusters$significant),
              ifelse(is.na(d$background_rate), NA, d$background_rate),
              d$no_elbow, d$all_chr_flagged))
}

## ---- RMSC curves (check the q* boundary) ----------------------------------
rc <- rbindlist(lapply(names(nm), function(k) data.table(method = nm[k], R[[k]]$rmsc)))
pk <- rc[, .SD[which.max(rejections)], by = method]
p1 <- ggplot(rc, aes(q, rejections, colour = method)) +
  geom_line() + geom_point(data = pk, size = 2.5) +
  geom_vline(data = pk, aes(xintercept = q, colour = method), linetype = 2, alpha = .5) +
  labs(x = "ld_w quantile threshold", y = "FDR discoveries",
       title = "RMSC selection curves (peak = q*)") +
  theme_bw(base_size = 9)
ggsave(file.path(mod, "fig_rmsc_curves.png"), p1, width = 150, height = 85, units = "mm", dpi = 130)
cat("\nRMSC peak location:\n"); print(pk)

## ---- cluster concordance (clusters identical across methods) --------------
base <- R$emx_perm$clusters[, .(cluster, Chr, start, end, n)]
mk_sig <- function(r) r$clusters[, .(cluster, sig = significant, q = emp_q, n_sig)]
m <- Reduce(function(a, b) merge(a, b, by = "cluster", all = TRUE),
            list(base,
                 setnames(mk_sig(R$emx_perm), c("cluster","sig_ep","q_ep","nsig_emx")),
                 setnames(mk_sig(R$emx_bg),   c("cluster","sig_eb","q_eb","nsig_emx2")),
                 setnames(mk_sig(R$lfmm_bg),  c("cluster","sig_lb","q_lb","nsig_lfmm"))))
m[, `:=`(nsig_emx2 = NULL)]
for (c in c("sig_ep","sig_eb","sig_lb")) m[is.na(get(c)), (c) := FALSE]

cat("\n=== concordance (counts of significant clusters) ===\n")
cat(sprintf("EMMAX perm : %d\nEMMAX bg   : %d\nLFMM  bg   : %d\n",
            sum(m$sig_ep), sum(m$sig_eb), sum(m$sig_lb)))
cat(sprintf("EMMAX perm & bg (both)      : %d\n", sum(m$sig_ep & m$sig_eb)))
cat(sprintf("all three                   : %d\n", sum(m$sig_ep & m$sig_eb & m$sig_lb)))
cat(sprintf("EMMAX(either) & LFMM        : %d\n", sum((m$sig_ep|m$sig_eb) & m$sig_lb)))
cat(sprintf("EMMAX-only (not LFMM)       : %d\n", sum((m$sig_ep|m$sig_eb) & !m$sig_lb)))
cat(sprintf("LFMM-only (not EMMAX)       : %d\n", sum(m$sig_lb & !(m$sig_ep|m$sig_eb))))

## ---- top regions: union of significant, ranked by best emp_q --------------
m[, nsig_any := pmax(nsig_emx, nsig_lfmm, na.rm = TRUE)]
m[, best_q := pmin(q_ep, q_eb, q_lb, na.rm = TRUE)]
top <- m[sig_ep | sig_eb | sig_lb][order(best_q, -n)]
top[, Mb := sprintf("%.2f-%.2f", start/1e6, end/1e6)]
cat("\n=== significant clusters (union), top 25 by best emp_q ===\n")
print(top[1:min(25,.N), .(cluster, Chr, Mb, n, nsig_emx, nsig_lfmm,
                          q_ep = round(q_ep,3), q_eb = round(q_eb,3), q_lb = round(q_lb,3))])

saveRDS(m, file.path(mod, "cluster_concordance.rds"))
cat("\nwrote fig_rmsc_curves.png + cluster_concordance.rds\n")
