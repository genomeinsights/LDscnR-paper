library(data.table)
library(igraph)
library(parallel)
library(ggplot2)
rho_ld_grid <- c(seq(0.1, 0.95, by = 0.1), 0.99)

qt_grid <- c(seq(0.1, 0.95, by = 0.05), 0.99)

alpha_grid <- c(
  0.05,
  0.01,
  0.005,
  0.001,
  0.0005,
  0.0001
)

cores = 8
folder <- "./single_SNP_results/"
out_folder <- "./OR_performance/"
files <- list.files(folder,full.names = TRUE)
done <- list.files(out_folder,full.names = TRUE)

#files <- files[!basename(files) %in% basename(done)]
#file <- files[1]


for(file in done){
  message("Processing ",basename(file))
  data <- readRDS(paste0(out_folder,basename(file)))

  param_grid <- data$param_grid
  r2_F_ld <- data$r2_F_ld
  meth <- data$emx
  perf <- meth$perf

  perf_full <- merge(
    param_grid,
    perf,
    by = c("qt", "alpha", "rho", "rho_ld"),
    all.x = TRUE
  )

  FN_max <- perf_full[TP==0,max(FN,na.rm=TRUE)]
  perf_full[is.na(FN),FN:=FN_max]
  perf_full[is.na(TP),TP:=0]
  perf_full[is.na(FP),FP:=0]
  perf_full[,precision := TP / (TP + FP)]
  perf_full[,recall := TP / (TP + FN)]
  perf_full[,PR := precision * recall]
  perf_full[is.nan(precision),precision := 0]
  perf_full[is.nan(PR),PR := 0]

  baseline <- perf_full[alpha==0.05 & rho == "none",
    .(PR0 = mean(PR),
      TP0 = mean(TP),
      FP0 = mean(FP)),
    by = .(alpha, rho_ld)
  ]

  dt <- merge(
    perf_full[rho != "none" & alpha==0.05],
    baseline,
    by = c("alpha","rho_ld")
  )

  # C_scores <- tmp[, .(
  #   C = uniqueN(param_id) / nrow(param_grid)
  # ), by = marker]


  dt[, delta_PR := PR - PR0]
  dt[, delta_TP := TP - TP0]
  dt[, delta_FP := FP - FP0]


  dt2 <- merge(
    dt,
    r2_F_ld,
    by = c("rho","qt")
  )

  dt2[is.na(delta_PR)] <- 0
  #tmp2[,hist(delta_PR)]

  dt2[,E:=R2_F_ld*(n_loci/nrow(map))^2]
  summary(lm(delta_PR ~ E, data = dt2))

  ggplot(dt2,aes(R2_F_ld,delta_PR)) +
    geom_point() +
    geom_vline(xintercept = median(dt2$R2_F_ld,na.rm=TRUE))
  #E > median(E,na.rm=TRUE) &


  hm_red <- dt2[ ,
                 .(mean_delta_PR = mean(delta_PR, na.rm = TRUE),
                   mean_delta_TP = mean(delta_TP, na.rm = TRUE),
                   mean_delta_FP = mean(delta_FP, na.rm = TRUE),
                   mean_n_loci = mean(n_loci, na.rm = TRUE),
                   mean_cor = mean(R2_F_ld, na.rm = TRUE),
                   mean_precision = mean(precision, na.rm = TRUE),
                   mean_recall = mean(recall, na.rm = TRUE),
                   mean_TP = mean(TP, na.rm = TRUE),
                   mean_FP = mean(FP, na.rm = TRUE),
                   mean_E = mean(E, na.rm = TRUE)),

                 by = .(rho, qt)
  ]

  hm_red[, rho := factor(rho, levels = sort(unique(rho)))]
  hm_red[, qt := factor(qt, levels = sort(unique(qt)))]
  hm_red[, better_than_full := mean_delta_PR > 0]

  hm_full <- dt2[ ,
                  .(mean_delta_PR = mean(delta_PR, na.rm = TRUE),
                    mean_delta_TP = mean(delta_TP, na.rm = TRUE),
                    mean_delta_FP = mean(delta_FP, na.rm = TRUE),
                    mean_n_loci = mean(n_loci, na.rm = TRUE),
                    mean_cor = mean(R2_F_ld, na.rm = TRUE),
                    mean_precision = mean(precision, na.rm = TRUE),
                    mean_recall = mean(recall, na.rm = TRUE),
                    mean_TP = mean(TP, na.rm = TRUE),
                    mean_FP = mean(FP, na.rm = TRUE),
                    mean_E = mean(E, na.rm = TRUE)),

                  by = .(rho, qt)
  ]

  hm_full[, rho := factor(rho, levels = sort(unique(rho)))]
  hm_full[, qt := factor(qt, levels = sort(unique(qt)))]
  hm_full[, better_than_full := mean_delta_PR > 0]

  p1 <- ggplot(hm_full, aes(rho, qt, fill = mean_delta_TP)) +
    geom_tile(color="white",fill="white") +
    geom_tile(data=hm_red,color="white") +
    geom_tile(
      data = hm_red[mean_delta_PR > 0],
      fill = NA,
      color = "black",
      linewidth = 0.8
    ) +
    scale_fill_viridis_c(option="turbo") +
    theme_bw() +
    ggtitle("E>median(E)")

  p1

}
