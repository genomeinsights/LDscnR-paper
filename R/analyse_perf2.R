make_perf <- function(data, method, param_grid) {

  x <- data[[method]]

  if (is.null(x) || is.null(x$perf) || nrow(x$perf) == 0) {
    out <- copy(param_grid)

    out[, `:=`(
      method = method,
      TP = 0L,
      FP = 0L,
      FN = 0L,
      precision = 0,
      recall = 0,
      PR = 0
    )]

    return(out[])
  }

  out <- merge(
    param_grid,
    x$perf,
    by = c("qt", "alpha", "rho", "rho_ld"),
    all.x = TRUE
  )

  out[, method := method]

  out[is.na(TP), TP := 0L]
  out[is.na(FP), FP := 0L]
  out[is.na(FN), FN := 0L]
  out[is.na(precision), precision := 0]
  out[is.na(recall), recall := 0]
  out[is.na(PR), PR := 0]

  out[]
}

file <- files[20]

performances <- rbindlist(lapply(files, function(file) {

  message("Processing ", basename(file))

  data <- readRDS(paste0(out_folder,basename(file)))

  sim_par <- strsplit(basename(file), "_", fixed = TRUE)[[1]]

  param_grid <- data$param_grid
  r2_F_ld <- data$r2_F_ld

  perf_emx <- make_perf(data, "emx", param_grid)
  perf_lfmm <- make_perf(data, "lfmm", param_grid)

  perf <- rbindlist(list(perf_emx, perf_lfmm), fill = TRUE)

  perf[, `:=`(
    sel_int = sim_par[1],
    disp_kern = sim_par[2],
    rep = gsub(".rds", "", sim_par[3], fixed = TRUE)
  )]

  merge(
    perf,
    r2_F_ld,
    by = c("rho", "qt"),
    all.x = TRUE
  )

}), fill = TRUE)

perf_full <- copy(performances)

n_true_tab <- perf_full[!is.na(FN), .(
  n_true = max(TP + FN, na.rm = TRUE)
), by = .(sel_int, disp_kern, rep, method)]

perf_full[n_true_tab, on = .(sel_int, disp_kern, rep, method),
          n_true := i.n_true]

perf_full[is.na(TP), TP := 0L]
perf_full[is.na(FP), FP := 0L]
perf_full[is.na(FN), FN := n_true]


perf_full[, precision := fifelse(TP + FP > 0, TP / (TP + FP), 0)]
perf_full[, recall := fifelse(TP + FN > 0, TP / (TP + FN), 0)]
perf_full[, PR := precision * recall]


perf_full[,baseline := ifelse(rho=="none",TRUE,FALSE)]

ggplot(perf_full[qt <=0.9 & qt >=0.2 , .(mean_PR = mean(PR)), by=.(rho_ld,disp_kern,baseline,method,alpha)],aes(rho_ld,mean_PR,linetype=baseline,col=disp_kern))+
  geom_line() +
  facet_grid(method~alpha)




baseline <- perf_full[alpha==0.05 & rho == "none",
                      .(PR0 = mean(PR),
                        TP0 = mean(TP),
                        FP0 = mean(FP)),
                      by = .(alpha, rho_ld,sel_int,disp_kern,rep,method)
]


dt <- merge(
  perf_full[rho != "none" & alpha==0.05],
  baseline,
  by = c("alpha","rho_ld","sel_int","disp_kern","rep","method")
)


dt[, delta_PR := PR - PR0]
dt[, delta_TP := TP - TP0]
dt[, delta_FP := FP0 - FP ]


#dt[,p_improved := mean(delta_PR > 0)]
ggplot(dt[,mean(PR),by=.(rho_ld,disp_kern)],aes(rho_ld,V1,col=disp_kern))+
  geom_line()


hm_full <- dt[rho_ld==0.9,
                .(mean_delta_PR = mean(delta_PR, na.rm = TRUE),
                  mean_delta_TP = mean(delta_TP, na.rm = TRUE),
                  mean_delta_FP = -mean(delta_FP, na.rm = TRUE),
                  p_improved = mean(delta_PR > 0),
                  p_improved_or_equal = mean(delta_PR >= 0),
                  #p_improved_FP = mean(delta_FP z= 0),
                  mean_n_loci = mean(n_loci, na.rm = TRUE),
                  mean_cor = mean(R2_F_ld, na.rm = TRUE),
                  mean_precision = mean(precision, na.rm = TRUE),
                  mean_recall = mean(recall, na.rm = TRUE),
                  mean_TP = mean(TP, na.rm = TRUE),
                  mean_FP = mean(FP, na.rm = TRUE),
                  mean_E = mean(E, na.rm = TRUE)),

                by = .(rho, qt,sel_int,disp_kern,method)
]

rho_levels <- hm_full[, sort(unique(as.numeric(as.character(rho))))]
qt_levels  <- hm_full[, sort(unique(as.numeric(as.character(qt))))]

hm_full[, rho := factor(rho, levels = as.character(rho_levels))]
hm_full[, qt := factor(qt, levels = as.character(qt_levels))]


rho_breaks <- levels(hm_full$rho)[seq(1, length(levels(hm_full$rho)), 2)]
qt_breaks  <- levels(hm_full$qt)[seq(1, length(levels(hm_full$qt)), 2)]

hm_full[,sel_int := c("high","medium","low")[match(sel_int,c("V0.5","V1","V2"))]]
hm_full[,gene_flow := c("high","medium","low")[match(disp_kern,c("c1","c1.5","c2"))]]
hm_full[,sim_scenario := factor(paste(sel_int,gene_flow,sep=" | "))]
hm_full[,sim_scenario := factor(sim_scenario,
                                  levels=c("high | high", "high | medium", "high | low",
                                           "medium | high", "medium | medium","medium | low",
                                           "low | high", "low | medium","low | low"))]
tile_template <- function(hm_full,fill){
  #fill <- "mean_delta_PR"
  p1 <- ggplot(hm_full, aes(rho, qt, fill = !!rlang::sym(fill))) +
    geom_tile(data=hm_full) +
    facet_grid(method~sim_scenario)+
    scale_x_discrete(breaks = rho_breaks) +
    scale_y_discrete(breaks = qt_breaks) +
    # geom_tile(
    #   data = hm_full[p_improved_or_equal > 0.75],
    #   fill = NA,
    #   color = "grey40",
    #   linewidth = 0.25
    # ) +
    geom_tile(
      data = hm_full[p_improved > 0.5],
      fill = NA,
      color = "black",
      linewidth = 0.5
    ) +
    #scale_fill_viridis_c(option="turbo") +
    # geom_contour(
    #   aes(z = p_improved),
    #   breaks = 0.8,
    #   color = "black",
    #   linewidth = 1
    # )+
    theme_bw() +
    theme(axis.text.x = element_text(angle=90),
          strip.background = element_blank())
}



plots <- lapply(c("mean_delta_PR","mean_delta_TP","mean_delta_FP"),function(fl){

  tile_template(hm_full,fl)

})

(plots[[1]] + scale_fill_gradientn(
  colours = c("blue3", "white", "orange", "red3"),
  values = scales::rescale(c(-0.15, 0, 0.05, 0.15)),
  limits = c(-0.15, 0.15),
  oob = scales::squish, name=expression(Delta*PR))) /

  ( plots[[2]]+ scale_fill_gradientn(
    colours = c("blue3", "white", "orange", "red3"),
    values = scales::rescale(c(-3, 0, 1, 3)),
    limits = c(-3, 3),
    oob = scales::squish, name=expression(Delta*TP))) /

  (plots[[3]] + scale_fill_gradientn(
    colours = c("blue3", "white", "orange", "red3"),
    values = scales::rescale(c(-3, 0, 1, 3)),
    limits = c(-3, 3),
    oob = scales::squish, name=expression(-Delta*FP)))


