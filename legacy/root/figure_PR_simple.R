performances

keep <- perf_full[,alpha==0.05]
dt <- perf_full[keep
  ,
  .(mean_PR = mean(PR, na.rm = TRUE)),
  by = .(method, rho, qt,alpha,sim_scenario)
]

# baseline at qt = 0 from rho == "none"

baseline <- dt[
  rho == "none" ,
  .(method, qt, mean_PR,alpha,rho,sim_scenario)
]

baseline0 <- dt[
  rho == "none" & qt == 0 ,
  .(method, alpha, mean_PR,sim_scenario)
]

rho_levels <- sort(unique(dt[rho != "none", rho]))

baseline0_expanded <- CJ(
  method = unique(dt$method),
  alpha = unique(dt$alpha),
  sim_scenario = unique(dt$sim_scenario),
  rho = rho_levels
)[
  baseline0,
  on = .(method, alpha,sim_scenario)
]

baseline0_expanded[, qt := 0]


# remove original qt=0 rows for rho != none
dt_plot <- dt[!(rho != "none" & qt == 0)]

# add corrected common starting point
dt_plot <- rbind(
  dt_plot,
  baseline0_expanded,
  fill = TRUE
)

ggplot(
  dt_plot[method %in% c("lfmm", "emx") & rho!="none" & rho %in% c(0.1,0.3,0.5,0.7,0.9)],
  aes(
    x = qt,
    y = mean_PR,
    colour = rho,
    group = rho
  )
) +
  geom_line(size = 1) +

  # baseline
  geom_line(
    data = baseline,
    aes(qt, mean_PR, linetype = "Random subsampling"),
    colour = "black",
    linewidth = 1.2,
    linetype=2
  ) +
  #facet_grid(rho~.)+
  scale_color_viridis_d(option="turbo")+

  facet_grid(sim_scenario~method, scales="free_y") +

  labs(
    x = "qt",
    y = "Mean PR",
    colour = expression(rho),
    title = "LD-informed filtering vs random subsampling"
  ) +

  theme_bw()
