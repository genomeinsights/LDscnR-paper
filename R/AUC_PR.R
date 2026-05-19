auc_trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

compute_auc_prstar <- function(PR, n_boot = 500, k_max = 500) {
  PR <- PR[is.finite(PR)]

  curves <- replicate(n_boot, {
    pr_sub <- sample(PR, k_max, replace = FALSE)
    cummax(pr_sub)
  })

  mean_curve <- rowMeans(curves)

  x <- seq_len(k_max)

  auc <- auc_trapz(x, mean_curve)

  max_pr <- max(PR, na.rm = TRUE)

  data.table(
    AUC_PRstar = auc,
    max_PR = max_pr,
    R = auc / (max_pr * (k_max - 1))
  )
}

R = AUC_PRstar / (max_PR * (k_max - 1))

x <- seq(0, 1, length.out = k_max)
auc <- auc_trapz(x, mean_curve)
R <- auc / max_pr
