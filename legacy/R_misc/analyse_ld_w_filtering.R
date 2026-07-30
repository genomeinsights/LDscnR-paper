library(ggplot2)
library(patchwork)
library(reshape2)
library(data.table)
??melt
# How number of outliers change with increaased LD-filtering -------------------
emp_dat <- outliers_100_perm_regional[,.(rho,th_ldw, r2_th,l_min,EMX,LFMM)]

perm_data_long <- melt(
  outliers_100_perm[],
  id.vars = c("r2_th", "l_min", "th_ldw", "rho"),
  measure.vars = patterns("perm"),
  variable.name = "perm_rep",
  value.name = "emx_perm")
perm_data_long[,"Number of outliers":=lengths(emx_perm)]
dt_perm <- perm_data_long[,.("Number of outliers"=mean(`Number of outliers`)),by=.(q_t=th_ldw,rho=as.numeric(rho))]


emp_dat[,"Number of outliers":=lengths(EMX)]
dt_EMX <- emp_dat[,.("Number of outliers"=mean(`Number of outliers`)),by=.(q_t=th_ldw,rho=as.numeric(rho))]
emp_dat[,"Number of outliers":=lengths(LFMM)]
dt_LFMM <- emp_dat[,.("Number of outliers"=mean(`Number of outliers`)),by=.(q_t=th_ldw,rho=as.numeric(rho))]

#dt_perm[,hist(`Number of outliers`)]

x_breaks <- seq(0.25, 1.00, by = 0.10)

base_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=90)
  )

base_scale <- scale_colour_viridis_c(
  option = "turbo",
  name = expression("Quantile\nthreshold"~(q[t]))
)

make_plot <- function(dat) {
  ggplot(dat[q_t >= 0.9],
         aes(rho, `Number of outliers`,
             col = q_t, group = factor(q_t))) +
    geom_line() +
    base_scale +
    scale_x_continuous(breaks = x_breaks) +
    xlab(expression("Window size for"~ld[w]~(rho))) +
    base_theme
}

p1 <- make_plot(dt_perm) + ggtitle("A | Simulated")
p2 <- make_plot(dt_EMX)  + ggtitle("B | EMMAX")
p3 <- make_plot(dt_LFMM) + ggtitle("C | LFMM")

(p1 | p2 | p3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


map_C_perm <- summarise_stability(
  outliers = outliers_100_perm[l_min>5 & n_loci>1 & th_ldw>0.8 & rho>0.6],
  map = map_3sp_with_perm,
  p_names = names(p_cols_perm)
)



#outliers_100_perm_regional <- copy(outliers_100_perm)
outliers_100_perm <- outliers_100_perm_grm_H2_05
map_C_perm <- summarise_stability(
  outliers = outliers_100_perm[l_min>5 & n_loci>1 & th_ldw>0.8 & rho>0.6],
  map = map_3sp_with_perm,
  p_names = names(p_cols_perm)
)


DT_long_C <- melt(
  map_C_perm,
  id.vars = c("marker","Chr","Pos"),
  measure.vars = patterns("C_"),
  variable.name = "perm_rep",
  value.name = "C_perm"
)


#DT_long_C[,sum(C_perm),by=perm_rep]


perm_C <- DT_long_C[,.(C_mean=mean(C_perm),C_max=max(C_perm),sd=sd(C_perm)),by=.(marker,Chr,Pos)]
perm_C[,plot(C_max)]


## Use outliers_100_perm_regional for this
map_C_EMX <- summarise_stability(
  outliers = outliers_100_perm_regional[l_min>5 & n_loci>1 & th_ldw>0.8 & rho>0.6],
  map = map_3sp,
  p_names = "EMX"
)


DT_long <- melt(
  outliers_100_perm[l_min>5 & n_loci>1 & th_ldw>0.8 & rho>0.6],
  id.vars = c("r2_th", "l_min", "th_ldw", "rho"),
  measure.vars = patterns("perm"),
  variable.name = "perm_rep",
  value.name = "emx_perm")


outliers_perm <- unique(unlist(DT_long$emx_perm))
outliers_emx <- unique(unlist(outliers_100_perm_regional[l_min>5 & n_loci>1 & th_ldw>0.8 & rho>0.6]$EMX))
joint_outliers <- unique(c(outliers_perm,outliers_emx))

#joint_outliers <- map_C[C_PC2>0 | C_PC1>0,marker]

#length(joint_outliers)

el_joint_outl <- precompute_LD_edges(
  GTs = GTs_3sp[, , drop = FALSE],
  map = map_3sp[marker %in% joint_outliers],
  r2_min = 0.1,
  max_bp = 1e6,
  cores = cores
)


joint_ORs <- LD_igraph_components(bp_th = 1e6,r2_th = 0.25,markers = joint_outliers,el = el_joint_outl)[n_loci>=5]


joint_ORs[, in_EMX := marker %in% outliers_emx]
joint_ORs[, in_perm := marker %in% outliers_perm]

joint_ORs <- map_C_EMX[,.(marker,C_EMX)][joint_ORs,on="marker"]
joint_ORs <- perm_C[,.(marker,C_max)][joint_ORs,on="marker"]


OR_class <- joint_ORs[
  ,
  .(
    n_loci = .N,
    any_EMX  = any(in_EMX),
    any_perm = any(in_perm)
  ),
  by = CL_id
]

OR_class[
  ,
  outlier_status := fifelse(
    any_EMX & any_perm, "Both",
    fifelse(any_EMX, "EMX_only", "perm_only")
  )
]

joint_ORs <- merge(
  joint_ORs,
  OR_class[, .(CL_id, outlier_status)],
  by = "CL_id"
)


OR_class <- joint_ORs[
  ,
  .(
    n_loci = .N,
    C_EMX  = max(C_EMX, na.rm=TRUE),
    C_perm = max(C_max, na.rm=TRUE),
    prop_EMX  = mean(in_EMX),
    prop_perm = mean(in_perm),
    any_EMX  = any(in_EMX),
    any_perm = any(in_perm)
  ),
  by = CL_id
]
OR_class[
  ,
  outlier_status := fifelse(
    any_EMX & any_perm, "Both",
    fifelse(any_EMX, "EMX_only", "perm_only")
  )
]
OR_class[, delta_C := C_EMX - C_perm]
OR_class[, delta_prop := prop_EMX - prop_perm]


DT_perm_long <- melt(
  outliers_100_perm[l_min>5 & n_loci>1 & th_ldw>0.8 & rho>0.6],
  id.vars = c("r2_th", "l_min", "th_ldw", "rho"),
  measure.vars = patterns("perm"),
  variable.name = "perm_rep",
  value.name = "marker"
)[
  ,
  .(marker = unlist(marker)),
  by = .(r2_th, l_min, th_ldw, rho, perm_rep)
]

perm_cols <- grep("perm", names(outliers_100_perm), value = TRUE)
n_perm_total <- length(perm_cols)

DT_perm_long <- melt(
  outliers_100_perm[l_min > 5 & n_loci > 1 & th_ldw > 0.8 & rho > 0.6],
  id.vars = c("r2_th", "l_min", "th_ldw", "rho"),
  measure.vars = perm_cols,
  variable.name = "perm_rep",
  value.name = "marker"
)[
  ,
  .(marker = unlist(marker)),
  by = .(r2_th, l_min, th_ldw, rho, perm_rep)
]

perm_region_hits <- merge(
  joint_ORs[, .(marker, CL_id)],
  DT_perm_long,
  by = "marker",
  allow.cartesian = TRUE
)

perm_recurrence <- perm_region_hits[
  ,
  .(
    n_perm_reps = uniqueN(perm_rep)
  ),
  by = CL_id
][
  ,
  prop_perm_reps := n_perm_reps / n_perm_total
]

OR_class <- merge(
  OR_class,
  perm_recurrence,
  by = "CL_id",
  all.x = TRUE
)

OR_class[
  is.na(n_perm_reps),
  `:=`(
    n_perm_reps = 0L,
    prop_perm_reps = 0
  )
]
OR_class[
  ,
  .(
    n_regions = .N,
    median_loci = median(n_loci),
    median_C_EMX = median(C_EMX),
    median_C_perm = median(C_perm),
    median_deltaC = median(delta_C),
    median_perm_recurrence = median(prop_perm_reps)
  ),
  by = outlier_status
]

## ----------------------------------------
## Parse marker coordinates
## ----------------------------------------

joint_ORs[
  ,
  c("Chr", "BP") := tstrsplit(marker, ":", fixed = TRUE)
]

joint_ORs[, BP := as.integer(BP)]

## ----------------------------------------
## Region summaries
## ----------------------------------------

region_coords <- joint_ORs[
  ,
  .(
    Chr = unique(Chr)[1],
    start = min(BP, na.rm = TRUE),
    end = max(BP, na.rm = TRUE),
    range_bp = max(BP, na.rm = TRUE) - min(BP, na.rm = TRUE),
    midpoint = round(mean(range(BP, na.rm = TRUE)))
  ),
  by = CL_id
]

## ----------------------------------------
## Merge into OR_class
## ----------------------------------------

OR_class <- merge(
  region_coords,
  OR_class,
  by = "CL_id",
  all.x = TRUE
)

setorder(OR_class,-C_EMX)
OR_class
#--------------------------------------------
# Summarize region classes
#--------------------------------------------

plot_dt <- OR_class[
  ,
  .(
    n_regions = .N,
    median_loci = median(n_loci),
    mean_loci = mean(n_loci),
    median_deltaC = median(delta_C),
    median_perm_recurrence = median(prop_perm_reps, na.rm = TRUE)
  ),
  by = outlier_status
]

# Order classes manually
plot_dt[, outlier_status := factor(
  outlier_status,
  levels = c("EMX_only", "Both", "perm_only")
)]

#--------------------------------------------
# Main barplot
#--------------------------------------------

ggplot(plot_dt,
       aes(x = outlier_status,
           y = n_regions,
           fill = outlier_status)) +

  geom_col(width = 0.75,
           color = "black") +

  # Add region counts
  geom_text(aes(label = n_regions),
            vjust = -0.8,
            size = 5) +

  # Add summary statistics inside bars
  geom_text(
    aes(
      label = paste0(
        "median loci = ", round(median_loci), "\n",
        "median ΔC = ", round(median_deltaC, 2), "\n",
        "perm recurrence = ", round(median_perm_recurrence, 2)
      )
    ),
    color = "white",
    size = 3.5,
    lineheight = 1.0,
    vjust = 1.1
  ) +

  scale_fill_manual(
    values = c(
      "EMX_only" = "#1f78b4",
      "Both" = "#ff7f00",
      "perm_only" = "#e31a1c"
    )
  ) +

  labs(
    x = NULL,
    y = "Number of outlier regions",
    fill = "Region class"
  ) +

  theme_bw(base_size = 14) +

  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


OR_class[,plot(delta_C)]
abline(h=0)

joint_ORs[,plot(C_max,C_EMX)]

map_joint_anal <- map_3sp[]


dt <- melt(map_joint_anal,measure.vars = c("C_EMX","C_max"))


dt[is.na(outlier_status),outlier_status:="non-outlier"]
library(ggplot2)

dt[outlier_status=="perm_only",plot(value)]

png("manh_emx_perm.png",height = 3,width = 15,res = 100,units = "in")
ggplot(
  ,
  aes(
    Pos,
    value,
    col = outlier_status,
    size = outlier_status
  )
) +
  geom_point(data = dt[outlier_status == "non-outlier"]) +
  geom_point(data = dt[outlier_status != "non-outlier"]) +
  facet_grid(variable ~ Chr, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    #legend.position = "none",
    strip.background = element_blank()

  )+
  scale_color_manual(
    values = c(
      "Both" = "firebrick4",
      "EMX_only" = "steelblue",
      "non-outlier" = "grey40",
      "perm_only" = "black"
    ),name="Outlier status"
  ) +
  scale_size_manual(
    values = c(
      "Both" = 0.75,
      "EMX_only" = 0.75,
      "non-outlier" = 0.25,
      "perm_only" = 0.75
    ),guide=NULL
  )
dev.off()
