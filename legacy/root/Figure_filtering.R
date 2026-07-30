
perf_full

# Filter to comparison range
pf <- copy(perf_full)[alpha == 0.05 & rho_ld >= 0.8]

# E threshold within each replicate / parameter combination
perf_full[, keep_E := E2>quantile(E2,0.25, na.rm = TRUE) & E2<quantile(E,0.75, na.rm = TRUE),
          by = .(rep, rho_ld, sim_scenario, method, rho,rep)]

#pf[, keep_E := E > E_median]
pf[rho == "none", keep_E := TRUE]


# Function to summarize baseline vs filtered
make_diff <- function(x, label) {

  baseline <- x[
    rho == "none",
    .(mean_PR_baseline = mean(PR, na.rm = TRUE)),
    by = .(rep, rho_ld, sim_scenario, method)
  ]

  filtered <- x[
    rho != "none",
    .(mean_PR_filtered = mean(PR, na.rm = TRUE)),
    by = .(rep, rho_ld, sim_scenario, method)
  ]

  out <- merge(
    baseline,
    filtered,
    by = c("rep", "rho_ld", "sim_scenario", "method")
  )

  out[, delta_PR := mean_PR_filtered - mean_PR_baseline]
  out[, "Pre filtering" := label]
  out
}

diff_all <- make_diff(pf, "all loci")
diff_highE <- make_diff(pf[keep_E == TRUE], "filtered")

dt_diff <- rbind(diff_all, diff_highE)

ggplot(dt_diff, aes(x = `Pre filtering`, y = delta_PR, fill = `Pre filtering`)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_grid(method ~ sim_scenario) +
  labs(
    x = NULL,
    y = "Difference in mean PR: filtered - baseline",

  ) +
  scale_x_discrete(
    labels = c(
      "all loci" = "unfiltered",
      #"filtered" = expression(E %in% IQR(E))
      "filtered" = expression(Q[0.25]~"<"~E~"<"~Q[0.25])
    ))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_blank())+
  theme(legend.position = "none")

ggplot(dt_diff, aes(x = method, y = delta_PR, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_grid(. ~ sim_scenario) +
  labs(
    x = NULL,
    y = "Difference in mean PR: filtered - baseline",

  ) +
  scale_x_discrete(
    labels = c(
      "all loci" = "unfiltered",
      #"filtered" = expression(E %in% IQR(E))
      "filtered" = expression(Q[0.25]~"<"~E~"<"~Q[0.25])
    ))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_blank())+
  theme(legend.position = "none")

# ggplot(dt_diff, aes(x = filtering, y = delta_PR, colour = filtering)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   stat_summary(fun = mean, geom = "point", size = 2) +
#   stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.15) +
#   facet_grid(method ~ sim_scenario) +
#   labs(
#     x = NULL,
#     y = "Mean ΔPR: filtered - baseline"
#   ) +
#   theme_bw() +
#   theme(legend.position = "none")

p_improved <- dt_diff[
  ,
  .(
    p_gt0 = mean(delta_PR > 0, na.rm = TRUE),
    p_ge0 = mean(delta_PR >= 0, na.rm = TRUE),
    mean_delta = mean(delta_PR, na.rm = TRUE)
  ),
  by = .(method, sim_scenario, `Pre filtering`)
]

ggplot(
  p_improved,
  aes(`Pre filtering`, sim_scenario, fill = p_ge0)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(p_ge0, 2)), size = 4) +
  facet_grid(method ~ .) +
  scale_fill_gradient2(
    low = "#d73027",
    mid = "white",
    high = "#1a9850",
    midpoint = 0.5,
    limits = c(0, 1),
    name = "P(ΔPR ≥ 0)"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

p_improved_wide <- dcast(
  p_improved,
  method + sim_scenario ~ `Pre filtering`,
  value.var = "p_ge0"
)

p_improved_wide[
  ,
  delta_p := `E > median(E)` - `all loci`
]

ggplot(
  p_improved_wide,
  aes(sim_scenario, method, fill = delta_p)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(delta_p, 2)), size = 5) +
  scale_fill_gradient2(
    low = "#d73027",
    mid = "white",
    high = "#1a9850",
    midpoint = 0,
    limits = c(-1,1),
    name = expression(Delta * "P")
  ) +
  theme_bw()
