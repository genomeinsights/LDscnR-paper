

res[, param_id := paste(alpha, qt, rho, rho_ld, l_min, sep = "_")]
res_split <- split(res, res$param_id)

true_positives <- map[true_pos == TRUE, marker]
n_true <- uniqueN(map[true_QTN == TRUE, marker])

perf <- rbindlist(mclapply(res_split, function(par_comb) {

  OR <- par_comb$OR

  true_pos_OR <- sapply(OR, function(mrk) {
    any(mrk %in% true_positives)
  })

  QTN_f <- sapply(OR, function(mrk) {
    tmp <- map[marker %in% mrk & !is.na(focal_QTN)]
    if (nrow(tmp) == 0) return(NA_character_)
    tmp[which.max(max_LD_with_QTN), focal_QTN]
  })

  max_ld <- sapply(OR, function(mrk) {
    tmp <- map[marker %in% mrk]
    max(tmp$max_LD_with_QTN, na.rm = TRUE)
  })

  eval_dt <- data.table(
    QTN_f = QTN_f,
    true_pos_OR = true_pos_OR,
    max_ld = max_ld
  )

  # keep only the best OR per focal QTN
  setorder(eval_dt, QTN_f, -max_ld)
  eval_dt[, keep_TP := !duplicated(QTN_f) & !is.na(QTN_f)]

  TP <- eval_dt[, sum(true_pos_OR & keep_TP, na.rm = TRUE)]
  FP <- eval_dt[, sum(!true_pos_OR | !keep_TP, na.rm = TRUE)]

  captured_QTNs <- unique(eval_dt[true_pos_OR & keep_TP, QTN_f])
  FN <- n_true - length(captured_QTNs)

  out <- data.table(
    par_comb[1, .(qt, alpha, rho, rho_ld, l_min)],
    TP = TP,
    FP = FP,
    FN = FN,
    precision = TP / (TP + FP),
    recall = TP / (TP + FN)
  )

  out[, PR := precision * recall]
  out

}, mc.cores = cores))

baseline <- perf[
  rho == "none",
  .(PR0 = mean(PR)),
  by = .(alpha, rho_ld, l_min)
]

tmp <- merge(
  perf[rho != "none"],
  baseline,
  by = c("alpha","rho_ld","l_min")
)

tmp[, delta_PR := PR - PR0]

tmp2 <- merge(
  tmp,
  r2_F_ld,
  by = c("rho","qt")
)


tmp2[,E:=R2_F_ld*(n_loci/nrow(map))]
summary(lm(delta_PR ~ E, data = tmp2))

ggplot(tmp2,aes(E,delta_PR)) +
  geom_point()
#alpha<2 & l_min==0
hm <- tmp2[,
           .(mean_delta_PR = mean(delta_PR, na.rm = TRUE),
             mean_n_loci = mean(n_loci, na.rm = TRUE),
             mean_cor = mean(R2_F_ld, na.rm = TRUE),
             mean_precision = mean(precision, na.rm = TRUE),
             mean_recall = mean(recall, na.rm = TRUE),
             mean_E = mean(E, na.rm = TRUE)),

           by = .(rho, qt)
]

hm[, rho := factor(rho, levels = sort(unique(rho)))]
hm[, qt := factor(qt, levels = sort(unique(qt)))]
hm[, better_than_full := mean_delta_PR > 0]


p1 <- ggplot(hm, aes(rho, qt, fill = mean_precision)) +
  geom_tile(color="white") +
  geom_tile(
    data = hm[mean_delta_PR > 0],
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_viridis_c(option="turbo") +
  theme_bw() +
  ggtitle("All parameter combinations")

p2 <- ggplot(hm, aes(rho, qt, fill = mean_recall)) +
  geom_tile(color="white") +
  geom_tile(
    data = hm[mean_delta_PR > 0],
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_viridis_c(option="turbo") +
  theme_bw() +
  ggtitle("All parameter combinations")


p3 <- ggplot(hm, aes(rho, qt, fill = mean_delta_PR)) +
  geom_tile(color="white") +
  geom_tile(
    data = hm[mean_delta_PR > 0],
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_viridis_c(option="turbo") +
  theme_bw() +
  ggtitle("All parameter combinations")

p_all <- p1 | p2 | p3


hm <- tmp2[alpha<2 & l_min==0 & rho_ld>0.5,
           .(mean_delta_PR = mean(delta_PR, na.rm = TRUE),
             mean_n_loci = mean(n_loci, na.rm = TRUE),
             mean_cor = mean(R2_F_ld, na.rm = TRUE),
             mean_precision = mean(precision, na.rm = TRUE),
             mean_recall = mean(recall, na.rm = TRUE),
             mean_E = mean(E, na.rm = TRUE)),
           by = .(rho, qt)
]

hm[, rho := factor(rho, levels = sort(unique(rho)))]
hm[, qt := factor(qt, levels = sort(unique(qt)))]
hm[, better_than_full := mean_delta_PR > 0]

p1 <- ggplot(hm, aes(rho, qt, fill = mean_precision)) +
  geom_tile(color="white") +
  geom_tile(
    data = hm[mean_delta_PR > 0],
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_viridis_c(option="turbo") +
  theme_bw() +
  ggtitle("alpha<2 & l_min==0")

p2 <- ggplot(hm, aes(rho, qt, fill = mean_recall)) +
  geom_tile(color="white") +
  geom_tile(
    data = hm[mean_delta_PR > 0],
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_viridis_c(option="turbo") +
  theme_bw() +
  ggtitle("alpha<2 & l_min==0")


p3 <- ggplot(hm, aes(rho, qt, fill = mean_delta_PR)) +
  geom_tile(color="white") +
  geom_tile(
    data = hm[mean_delta_PR > 0],
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_viridis_c(option="turbo") +
  theme_bw() +
  ggtitle("alpha<2 & l_min==0")


p_sub <- p1 | p2 | p3

p_all / p_sub

# library(patchwork)
# summary(lm(delta_PR ~ E, data = tmp2)) # r2
# summary(lm(delta_PR ~ E, data = tmp2[alpha<2 & l_min==0]))
#
# tmp2[,cor(delta_PR,E)^2]
# tmp2[alpha<2 & l_min==0,cor(delta_PR,E)^2]
#
#
# p1 | p2 | p3
#
# perf[alpha<2 & l_min==0 & rho=="none",mean(PR),by=.(rho_ld)]
# perf[alpha<2 & l_min==0 & rho!="none",mean(PR),by=.(rho_ld)]
#
#
# perf[,cor.test(rho_ld,PR)]
