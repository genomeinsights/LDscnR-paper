library(patchwork)
library(ggplot2)
library(data.table)
library(igraph)
library(parallel)
library(LDscnR)

# ----------------------------
# LD edge precomputation
# ----------------------------
precompute_LD_edges <- function(
    GTs,
    map,
    r2_min = 0.1,
    max_bp = Inf,
    cores = 1
) {
  markers <- intersect(colnames(GTs), map$marker)
  if (length(markers) == 0) stop("No overlapping markers between GTs and map.")

  map_sub <- copy(map[marker %in% markers])
  setkey(map_sub, Chr, Pos)

  chr_levels <- unique(map_sub$Chr)

  out <- mclapply(chr_levels, function(ch) {
    chr_map <- map_sub[Chr == ch]
    chr_markers <- chr_map$marker

    if (length(chr_markers) < 2) {
      return(data.table(
        Chr = ch,
        marker1 = chr_markers,
        marker2 = chr_markers,
        r2 = 1,
        dist_bp = 0
      ))
    }

    gts <- as.matrix(GTs[, chr_markers, drop = FALSE])
    storage.mode(gts) <- "double"

    R2 <- cor(gts, use = "pairwise.complete.obs")^2
    R2[is.na(R2)] <- 0
    diag(R2) <- 0

    idx <- which(R2 >= r2_min, arr.ind = TRUE)
    idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]

    if (nrow(idx) == 0) {
      return(data.table(
        Chr = ch,
        marker1 = chr_markers,
        marker2 = chr_markers,
        r2 = 1,
        dist_bp = 0
      ))
    }

    pos <- chr_map$Pos

    dt <- data.table(
      Chr = ch,
      marker1 = chr_markers[idx[, 1]],
      marker2 = chr_markers[idx[, 2]],
      r2 = R2[idx],
      dist_bp = abs(pos[idx[, 1]] - pos[idx[, 2]])
    )

    if (is.finite(max_bp)) dt <- dt[dist_bp <= max_bp]

    dt
  }, mc.cores = cores)

  out <- rbindlist(out, use.names = TRUE, fill = TRUE)
  setkey(out, Chr, marker1, marker2)

  out
}


# ----------------------------
# LD clusters from edge list
# ----------------------------
LD_igraph_components <- function(
    el,
    markers,
    r2_th = 0.8,
    bp_th = Inf
) {
  markers <- unique(markers)

  if (length(markers) == 0) {
    return(data.table(marker = character(), CL_id = integer(), n_loci = integer()))
  }

  if (length(markers) == 1) {
    return(data.table(marker = markers, CL_id = 1L, n_loci = 1L))
  }

  edges <- el[
    marker1 %in% markers &
      marker2 %in% markers &
      r2 >= r2_th
  ]

  if (is.finite(bp_th)) edges <- edges[dist_bp <= bp_th]

  if (nrow(edges) == 0) {
    return(data.table(
      marker = markers,
      CL_id = seq_along(markers),
      n_loci = 1L
    ))
  }

  g <- graph_from_data_frame(
    edges[, .(from = marker1, to = marker2)],
    directed = FALSE,
    vertices = data.table(name = markers)
  )

  comp <- components(g)

  clusters <- data.table(
    marker = names(comp$membership),
    CL_id = as.integer(comp$membership)
  )

  clusters[, n_loci := comp$csize[CL_id]]
  clusters
}


# ----------------------------
# Empty result, generalized
# ----------------------------
empty_result <- function(
    rho,
    th_ldw,
    p_names,
    r2_grid,
    lmin_grid
) {
  out <- CJ(r2_th = r2_grid, l_min = lmin_grid)

  out[, `:=`(
    th_ldw = th_ldw,
    rho = rho
  )]

  for (nm in p_names) {
    out[, (nm) := list(list(character()))]
  }

  out[]
}


# ----------------------------
# One grid point, generalized
# ----------------------------
run_one_grid <- function(
    map,
    el=NULL,
    ld_ws,
    rho,
    th_ldw,
    p_cols,
    p_names = names(p_cols),
    alpha = 0.05,
    r2_grid,
    lmin_grid,
    bp_th = Inf,
    cores
) {
  stopifnot(length(p_cols) == length(p_names))

  # Important: align ld_ws and map by marker before filtering
  common_markers <- intersect(map$marker, rownames(ld_ws))


  map_sub <- copy(map[marker %in% common_markers])
  ld_sub <- ld_ws[map_sub$marker, , drop = FALSE]

  keep_ld_w <- ld_sub[, rho] > quantile(ld_sub[, rho], th_ldw, na.rm = TRUE)

  keep <- ld_sub[, rho] > quantile(ld_sub[, rho], th_ldw, na.rm = TRUE)
  keep[is.na(keep)] <- FALSE
  #table(keep)

  if (!any(keep)) {
    return(cbind(empty_result(rho, th_ldw, p_names, r2_grid, lmin_grid),n_loci=length(which(keep))))
  }

  markers_keep <- map_sub[keep, marker]

  outliers <- setNames(vector("list", length(p_cols)), p_names)

  #i <- 1
  for (i in seq_along(p_cols)) {
    p_col <- p_cols[i]
    nm <- p_names[i]

    q <- p.adjust(unlist(map_sub[keep, ..p_col]), method = "fdr")
    outliers[[nm]] <- markers_keep[q < alpha]
  }

  if (length(unique(unlist(outliers))) == 0) {
    return(cbind(empty_result(rho, th_ldw, p_names, r2_grid, lmin_grid),n_loci=length(which(keep))))
  }


  if(is.null(el)){
    all_outliers <- unique(unlist(outliers))
    el <- precompute_LD_edges(
      GTs = GTs[, all_outliers, drop = FALSE],
      map = map_sub[marker %in% all_outliers],
      r2_min = 0.1,
      max_bp = 1e6,
      cores = 1
    )
  }
  #r2_th = 0.8
  out <- rbindlist(mclapply(r2_grid, function(r2_th) {
    clusters <- lapply(outliers, function(markers) {
      LD_igraph_components(
        el = el,
        markers = markers,
        r2_th = r2_th,
        bp_th = bp_th
      )
    })

    out <- rbindlist(lapply(lmin_grid, function(l_min) {
      row <- data.table(
        r2_th = r2_th,
        l_min = l_min,
        th_ldw = th_ldw,
        rho = rho
      )
  #l_min = 1
      #nm = p_names[1]
      for (nm in p_names) {
        tmp <- clusters[[nm]][n_loci >= l_min, ]
        cls <- split(tmp$marker,tmp$CL_id)

        if(length(cls)>1){
          row[, (nm) := list(cls)]
        }else{
          row[, (nm) := list(list(cls))]
        }
      }
      row

    }), fill = TRUE)
  },mc.cores=cores), fill = TRUE)

  out[,n_loci:=length(which(keep))]

}
# ----------------------------
# Potential outliers, generalized
# ----------------------------
get_potential_outliers <- function(
    map,
    ld_ws,
    th_ldw_grid,
    p_cols,
    alpha = 0.05
) {
  common_markers <- intersect(map$marker, rownames(ld_ws))

  #map <- copy(map[marker %in% common_markers])
  #ld <- ld_ws[map$marker, , drop = FALSE]

  potential <- character()

  for (rho in colnames(ld_ws)) {
    message("processing rho = ",rho)
    ld_vec <- ld_ws[, rho]

    potential <- unique(potential)
    for (th_ldw in th_ldw_grid) {
      keep <- ld_vec > quantile(ld_vec, th_ldw, na.rm = TRUE)

        keep[is.na(keep)] <- FALSE

        if (!any(keep)) next
        ##p_col <- p_cols[1]
        for (p_col in p_cols) {
          q <- p.adjust(unlist(map[keep, ..p_col]), method = "fdr")
          potential <- c(potential, map[keep, marker][q < alpha])
        }
      }
    }

  return(unique(potential))

}

# ----------------------------
# Summerize
# ----------------------------
summarise_stability <- function(outliers, map, p_names) {
  map_C <- copy(map)
  #nm <- "PC1"
  for (nm in p_names) {
    C <- outliers[, table(unlist(get(nm))) / .N]
    if(length(C)>0){
      C <- data.table(C)
      setnames(C, c("V1", "N"), c("marker", paste0("C_", nm)))
      map_C <- C[map_C, on = "marker"]
      map_C[is.na(get(paste0("C_", nm))), (paste0("C_", nm)) := 0]
    }
  }

  map_C[]
}
# ----------------------------
# Example use
# ----------------------------
# read in raw data


ld_ws_3sp <- precalculate_ld_w(seq(0,0.95,by=0.05),ld_decay_3sp)
rownames(ld_ws_3sp) <- map_3sp$marker
colnames(ld_ws_3sp)
ld_ws_3sp[is.na(ld_ws_3sp)] <- 0

th_ldw_grid <- 1 - 10^seq(log10(1), log10(0.01), length.out = 20)

p_cols <- c(
  EMX = "emx_p_GC",
  EMX_perm = "emx_p_perm",
  LFMM = "lfmm_P"
)

p_cols <- c(
  #EMX = "emx_p_GC",
  EMX_perm = "emx_p_perm"
  #LFMM = "lfmm_P"
)

# Regional permutation ----------------------------
emmax_perm_reginal <- readRDS("./3sp_data/emmax_perm_reginal.rds")
emmax_perm_reginal <- do.call(cbind,emmax_perm_reginal)

colnames(emmax_perm_reginal) <- paste0("emx_p_perm",1:ncol(emmax_perm_reginal))

map_3sp_with_perm <- cbind(map_3sp[,.(Chr,Pos,marker,emx_p_GC,lfmm_P)],emmax_perm_reginal)


p_cols <- c(
  EMX = "emx_p_GC",
  LFMM = "lfmm_P"
)
p_cols_perm <- colnames(emmax_perm_reginal)
names(p_cols_perm) <- colnames(emmax_perm_reginal)
p_cols_all <- p_cols#c(p_cols,p_cols_perm)


potential_outliers <- get_potential_outliers(
  map = map_3sp_with_perm,
  ld_ws = ld_ws_3sp,
  th_ldw_grid = th_ldw_grid,
  p_cols = p_cols_all,
  alpha = 0.05
)

length(potential_outliers)

el_potential <- precompute_LD_edges(
  GTs = GTs_3sp[, potential_outliers, drop = FALSE],
  map = map_3sp[marker %in% potential_outliers],
  r2_min = 0.1,
  max_bp = 1e6,
  cores = cores
)

r2_grid   <- seq(0.6,0.9,by=0.05)
lmin_grid <- c(1,5,10,20,40,80,160)


param_grid <- CJ(
  rho = colnames(ld_ws_3sp),
  th_ldw = th_ldw_grid
)

#i <- 500
outliers_100_perm <- rbindlist(
  mclapply(seq_len(nrow(param_grid)), function(i) {
    pars <- param_grid[i]
    cat(i,"..")
    out <- run_one_grid(
      map = map_3sp_with_perm,
      el = el_potential,
      ld_ws = ld_ws_3sp,
      rho = pars$rho,
      th_ldw = pars$th_ldw,
      p_cols = p_cols_all,
      alpha = 0.05,
      r2_grid = r2_grid,
      lmin_grid = lmin_grid,
      bp_th = Inf,
      cores = 1
    )
  },mc.cores=1),
  fill = TRUE
)
#saveRDS(outliers_100_perm,"outliers_3sp_100_perm.rds")
#q("no")

#save.image()


#tmp1 <- outliers_100_perm[r2_th==0.85,.("Number of outliers"=mean(lengths(EMX))),by=.(rho=factor(rho),th_ldw)]


tmp1 <- outliers_100_perm[r2_th==0.85 & l_min==5,.("Number of outliers"=mean(lengths(EMX))),by=.(rho=factor(rho),th_ldw)]
#tmp1[rho==0.99,rho:=1]
#tmp1[th_ldw==0.99,th_ldw:=1]
p1 <- ggplot(tmp1, aes(
  x = as.numeric(as.character(rho)),
  y = th_ldw,
  fill = `Number of outliers`
)) +
  scale_y_log10(
    breaks = th_ldw_grid,
    labels = round(th_ldw_grid,3)
  ) +
  geom_tile() +
  scale_fill_gradientn(
    colors = wes_palette("Zissou1", 100, type = "continuous")
  ) +
  labs(
    x = expression("LD window size relative to LD-decay ("*rho*")"),
    y = expression("Quantile trehshold for local LD"),
    fill = "No. outliers"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_blank())
p1

#sapply(outliers_100_perm[rho=="0.95",LFMM],length)

#outliers_100_perm[rho=="0.99"][21,]
#outliers_100_perm[r2_th==0.85]

#outliers_100_perm[l_min==25 & lengths(LFMM)>100][23,LFMM]
sapply(outliers_100_perm$EMX,function(x) is(x)[1]=="character")

#outliers_100_perm[sapply(outliers_100_perm$EMX,function(x) is(x)[1]=="character"),EMX]
tmp1.2 <- outliers_100_perm[r2_th==0.7 & l_min>=5,.("Number of outliers"=mean(lengths(LFMM))),by=.(rho=factor(rho),th_ldw)]
#tmp1.2[which.max()]
#tmp1.2[rho==0.99,rho:=1]
#tmp1.2[th_ldw==0.99,th_ldw:=1]

dt[, th_ldw_f := factor(
  th_ldw,
  levels = th_ldw_grid,
  ordered = TRUE
)]

#dt[,th_ldw]

tmp1.2[,data:="LFMM"]
tmp1[,data:="EMMAX"]

dt <- rbind(tmp1,tmp1.2)

dt[, th_ldw_f := factor(
  round(th_ldw,3),
  levels = round(th_ldw_grid,3),
  ordered = TRUE
)]

p1 <- ggplot(dt, aes(
  x = as.numeric(as.character(rho)),
  y = th_ldw_f,
  fill = `Number of outliers`
)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(. ~ data) +
  scale_fill_gradientn(
    colors = wes_palette("Zissou1", 100, type = "continuous"),
    name = "No. outlier\nregions"
  ) +
  # scale_y_continuous(
  #   labels = scales::percent_format(accuracy = 1)
  # ) +
  labs(
    x = expression("LD window size relative to LD decay (" * rho * ")"),
    y = "Local LD filtering threshold"#,
    #title = "Number of outlier regions across LD-filtering parameters"
  ) +
  # scale_x_continuous(
  #   breaks = seq(0.2, 2, by = 0.1)
  # ) +
  # scale_y_discrete(
  #   breaks = th_ldw_grid[c(1, 6, 11, 16, 20)],
  #   labels = scales::number(th_ldw_grid[c(1, 6, 11, 16, 20)], accuracy = 0.01)
  # ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(3, 3, 3, 3)
  )
p1
ggsave(p1,filename="the_heatmap.png",height = 5,width = 12,units = "in",dpi=300)
# #library(wesanderson)
# p1 <- ggplot(tmp, aes(rho,`Number of outliers`,col=th_ldw,group=th_ldw)) +
#   geom_line(linewidth=1.5) +
#   #scale_color_viridis_c(option="turbo",name=expression("Quantile \nthreshold"~q[t])) +
#   scale_color_gradientn(
#     colors = wes_palette("Zissou1", 100, type = "continuous"),name=expression("Quantile \nthreshold"~q[t])
#   ) +
#   theme_bw() +
#   scale_fill_manual(values = wes_palette("Zissou1")) +
#   theme(axis.text.x=element_text(angle=90)) +
#   xlab(expression("Window size for local LD"~rho[ld[w]]))
# p1

DT_long <- melt(
  outliers_100_perm,
  id.vars = c("r2_th", "l_min", "th_ldw", "rho"),
  measure.vars = patterns("perm"),
  variable.name = "perm_rep",
  value.name = "emx_perm"
)

tmp2 <- DT_long[th_ldw>=0.7,.("Number of outliers"=mean(lengths(emx_perm))),by=.(rho=factor(rho),th_ldw)]

tmp2[rho=="0.99",rho:="1"]
p2 <-ggplot(tmp2, aes(
  x = as.numeric(as.character(rho)),
  y = th_ldw,
  fill = `Number of outliers`
)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = wes_palette("Zissou1", 100, type = "continuous")
  ) +
  labs(
    x = expression("LD window size relative to LD-decay ("*rho*")"),
    y = expression("Quantile trehshold for local LD"),
    fill = "No. outliers"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_blank())

# #p1 | p2
#
# p2 <- ggplot(tmp, aes(rho,`Number of outliers`,col=th_ldw,group=th_ldw)) +
#   geom_line(linewidth=1.5) +
#   scale_color_viridis_c(name=expression("Quantile \nthreshold"~q[t])) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90)) +
#   ylab("Mean number of outliers")+
#   xlab(expression("Window size for local LD"~rho[ld[w]]))
# (p1 + ggtitle("(A) Observed data")) | (p2 + ggtitle("(B) Mean from 100 permuted data set"))
#
#
# p1 <- ggplot(tmp, aes(r2_th,`No. outl.`,col=l_min,group=l_min)) +
#   geom_line()+
#   scale_color_gradientn(colours = c("steelblue", "forestgreen", "orange", "firebrick4")) +
#   theme_bw()
#
# tmp <- outliers_100_perm[th_ldw>0.7,.("Number of outliers"=mean(lengths(LFMM))),by=.(rho=factor(rho),th_ldw)]
#
# p2 <- ggplot(tmp, aes(rho,`Number of outliers`,col=th_ldw,group=th_ldw)) +
#   geom_line(linewidth=2) +
#   scale_color_viridis_c(option="turbo",name=expression("Quantile \nthreshold"~q[t])) +
#   theme_bw() +
#   xlab(expression("Window size for local LD"~rho[ld[w]]))
# p2
#
# p1 | p2

DT_long <- melt(
  outliers_100_perm,
  id.vars = c("r2_th", "l_min", "th_ldw", "rho"),
  measure.vars = patterns("perm"),
  variable.name = "perm_rep",
  value.name = "emx_perm"
)

tmp2 <- DT_long[,.("No. outl."=mean(lengths(emx_perm))),by=.(r2_th=factor(r2_th),l_min)]

p3 <- ggplot(tmp2, aes(r2_th,`No. outl.`,col=l_min,group=l_min)) +
  geom_line()+
  scale_color_gradientn(colours = c("steelblue", "forestgreen", "orange", "firebrick4")) +
  theme_bw()

tmp <- DT_long[,.("No. outl."=mean(lengths(emx_perm))),by=.(rho=factor(rho),th_ldw)]

p4 <- ggplot(tmp, aes(rho,`No. outl.`,col=th_ldw,group=th_ldw)) +
  geom_line() +

  scale_color_gradientn(colours = c("steelblue", "forestgreen", "orange", "firebrick4")) +
  theme_bw()

(p1+ggtitle("EMMAX") | p2+ggtitle("EMMAX")) / (p3+ggtitle("EMMX perm.") |  p4+ggtitle("EMMX perm."))


map_C <- summarise_stability(
  outliers = outliers[l_min>5 & n_loci>1 ],
  map = map_3sp,
  p_names = names(p_cols)
)

par(mfcol=c(2,1))
map_C[,plot(C_EMX_perm,pch=20,ylim=c(0,0.2))]
map_C[,plot(C_EMX,pch=20,ylim=c(0,0.2))]
#map_C[,plot(C_LFMM,pch=20,ylim=c(0,0.5))]

# rests ------------------------

prep_outlier_summary <- function(outliers) {
  dt <- copy(outliers)
  setDT(dt)

  long <- rbindlist(list(
    dt[, .(
      method = "LFMM",
      r2_th,
      l_min,
      th_ldw,
      rho = as.numeric(rho),
      n_loci,
      n_outliers = lengths(LFMM)
    )],
    dt[, .(
      method = "EMMAX",
      r2_th,
      l_min,
      th_ldw,
      rho = as.numeric(rho),
      n_loci,
      n_outliers = lengths(EMX)
    )],
    dt[, .(
      method = "EMMAX permuted",
      r2_th,
      l_min,
      th_ldw,
      rho = as.numeric(rho),
      n_loci,
      n_outliers = lengths(EMX_perm)
    )]
  ), use.names = TRUE)

  long[, prop_outliers := n_outliers / n_loci]

  long[]
}

summarise_outliers <- function(long) {
  long[, .(
    mean_outliers = mean(n_outliers, na.rm = TRUE),
    mean_prop = mean(prop_outliers, na.rm = TRUE),
    mean_n_loci = mean(n_loci, na.rm = TRUE)
  ), by = .(method, r2_th, l_min, th_ldw, rho)]
}

plot_by_r2 <- function(sum_dt, y = "mean_outliers") {
  ggplot(
    sum_dt,
    aes(
      x = r2_th,
      y = .data[[y]],
      colour = l_min,
      group = l_min
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~ method, scales = "free_y") +
    scale_colour_gradientn(
      colours = c("steelblue", "forestgreen", "orange", "firebrick4")
    ) +
    labs(
      x = expression(r^2~"clustering threshold"),
      y = ifelse(y == "mean_outliers", "No. outliers", "Outlier proportion"),
      colour = expression(l[min])
    ) +
    theme_bw()
}

plot_by_ldw <- function(sum_dt, y = "mean_outliers") {
  ggplot(
    sum_dt,
    aes(
      x = rho,
      y = .data[[y]],
      colour = th_ldw,
      group = th_ldw
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~ method, scales = "free_y") +
    scale_colour_gradientn(
      colours = c("steelblue", "forestgreen", "orange", "firebrick4")
    ) +
    labs(
      x = expression(rho~"(LD decay distance parameter)"),
      y = ifelse(y == "mean_outliers", "No. outliers", "Outlier proportion"),
      colour = "LD-window\nquantile"
    ) +
    theme_bw()
}
long_out <- prep_outlier_summary(outliers_144)
######
sum_r2 <- long_out[
  rho == 0.25 & th_ldw == 0,
  .(
    mean_outliers = mean(n_outliers),
    mean_prop = mean(prop_outliers)
  ),
  by = .(method, r2_th, l_min)
]

sum_rho <- long_out[
  r2_th == 0.8 & l_min == 5,
  .(
    mean_outliers = mean(n_outliers),
    mean_prop = mean(prop_outliers)
  ),
  by = .(method, rho, th_ldw)
]

ggplot(sum_r2, aes(r2_th, mean_outliers, colour = l_min, group = l_min)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ method, scales = "free_y") +
  theme_bw()

sum_out <- summarise_outliers(long_out)

ggplot(sum_out,aes())

p_counts_r2  <- plot_by_r2(sum_out, y = "mean_outliers")
p_counts_ldw <- plot_by_ldw(sum_out, y = "mean_outliers")

p_props_r2   <- plot_by_r2(sum_out, y = "mean_prop")
p_props_ldw  <- plot_by_ldw(sum_out, y = "mean_prop")

(p_counts_r2 / p_counts_ldw) | (p_props_r2 / p_props_ldw)

emx_enrichment <- dcast(
  sum_out[method %in% c("EMMAX", "EMMAX permuted")],
  r2_th + l_min + th_ldw + rho ~ method,
  value.var = "mean_outliers"
)

emx_enrichment[, enrichment := (EMMAX + 1) / (`EMMAX permuted` + 1)]

ggplot(
  emx_enrichment,
  aes(
    x = rho,
    y = enrichment,
    colour = th_ldw,
    group = th_ldw
  )
) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_colour_gradientn(
    colours = c("steelblue", "forestgreen", "orange", "firebrick4")
  ) +
  labs(
    x = expression(rho),
    y = "EMMAX / permuted enrichment",
    colour = "LD-window\nquantile"
  ) +
  theme_bw()

outliers_144
tmp <- outliers_144[,.("No. outl."=mean(lengths(LFMM))),by=.(r2_th=factor(r2_th),l_min)]

p1 <- ggplot(tmp, aes(r2_th,`No. outl.`,col=l_min,group=l_min)) +
  geom_line()+
  scale_color_gradientn(colours = c("steelblue", "forestgreen", "orange", "firebrick4")) +
  theme_bw()

tmp <- outliers_144[,.("No. outl."=mean(lengths(LFMM))),by=.(rho=factor(rho),th_ldw)]

p2 <- ggplot(tmp, aes(rho,`No. outl.`,col=th_ldw,group=th_ldw)) +
  geom_line() +

  scale_color_gradientn(colours = c("steelblue", "forestgreen", "orange", "firebrick4")) +
  theme_bw()


tmp <- outliers_144[,.("No. outl."=mean(lengths(EMX))),by=.(r2_th=factor(r2_th),l_min)]

p3 <- ggplot(tmp, aes(r2_th,`No. outl.`,col=l_min,group=l_min)) +
  geom_line()+
  scale_color_gradientn(colours = c("steelblue", "forestgreen", "orange", "firebrick4")) +
  theme_bw()

tmp <- outliers_144[,.("No. outl."=mean(lengths(EMX))),by=.(rho=factor(rho),th_ldw)]

p4 <- ggplot(tmp, aes(rho,`No. outl.`,col=th_ldw,group=th_ldw)) +
  geom_line() +
  scale_color_gradientn(colours = c("steelblue", "forestgreen", "orange", "firebrick4")) +
  theme_bw()


tmp <- outliers_144[,.("No. outl."=mean(lengths(EMX_perm))),by=.(r2_th=factor(r2_th),l_min)]

p5 <- ggplot(tmp, aes(r2_th,`No. outl.`,col=l_min,group=l_min)) +
  geom_line()+
  scale_color_gradientn(colours = c("steelblue", "forestgreen", "orange", "firebrick4")) +
  theme_bw()

tmp <- outliers_144[,.("No. outl."=mean(lengths(EMX_perm))),by=.(rho=factor(rho),th_ldw)]

p5 <- ggplot(tmp, aes(rho,`No. outl.`,col=th_ldw,group=th_ldw)) +
  geom_line() +
  scale_color_gradientn(colours = c("steelblue", "forestgreen", "orange", "firebrick4")) +
  theme_bw()

(p1+ggtitle("LFMM") | p3+ggtitle("EMMAX")) / (p2+ggtitle("LFMM") |  p4+ggtitle("EMMAX")) / (p5+ggtitle("EMMAX permuted") |  p5+ggtitle("EMMAX permuted"))



map_C <- summarise_stability(
  outliers = outliers_54[l_min>5 & n_loci>1 & r2_th<=0.7 & th_ldw > 0.7 ],
  map = map_3sp,
  p_names = names(p_cols)
)


par(mfcol=c(2,1))
map_C[,plot(C_LFMM,pch=20)]
map_C[,plot(C_EMX,pch=20)]
# global permutation ------------------------------
emmax_perm_grm <- readRDS("./3sp_data/out_emx_grm_05.rds")
emmax_perm_grm <- do.call(cbind,emmax_perm_grm)

colnames(emmax_perm_grm) <- paste0("emx_p_perm",1:ncol(emmax_perm_grm))

map_3sp_with_perm_grm <- cbind(map_3sp[,.(Chr,Pos,marker,lfmm_P,emx_p_GC)],emmax_perm_grm)


p_cols_perm_grm <- colnames(emmax_perm_grm)
names(p_cols_perm_grm) <- colnames(emmax_perm_grm)

potential_outliers <- get_potential_outliers(
  map = map_3sp_with_perm_grm,
  ld_ws = ld_ws_3sp,
  th_ldw_grid = th_ldw_grid,
  p_cols = p_cols_perm_grm,
  alpha = 0.05
)


el_potential <- precompute_LD_edges(
  GTs = GTs_3sp[, potential_outliers, drop = FALSE],
  map = map_3sp[marker %in% potential_outliers],
  r2_min = 0.1,
  max_bp = 1e6,
  cores = cores
)

r2_grid   <- seq(0.6,0.9,by=0.1)
lmin_grid <- 5


param_grid <- CJ(
  rho = colnames(ld_ws_3sp)[-(1:5)],
  th_ldw = th_ldw_grid
)


outliers_100_perm_grm_H2_05 <- rbindlist(
  mclapply(seq_len(nrow(param_grid)), function(i) {
    pars <- param_grid[i]
    cat(i,"..")
    out <- run_one_grid(
      map = map_3sp_with_perm_grm,
      el = el_potential,
      ld_ws = ld_ws_3sp,
      rho = pars$rho,
      th_ldw = pars$th_ldw,
      p_cols = p_cols_all,
      alpha = 0.05,
      r2_grid = r2_grid,
      lmin_grid = lmin_grid,
      bp_th = Inf,
      cores = 1
    )
  },mc.cores=4),
  fill = TRUE
)
save.image()


DT_long <- melt(
  outliers_100_perm_grm_H2_05,
  id.vars = c("r2_th", "l_min", "th_ldw", "rho"),
  measure.vars = patterns("perm"),
  variable.name = "perm_rep",
  value.name = "emx_perm"
)

tmp3 <- DT_long[th_ldw>=0.7,.("Number of outliers"=mean(lengths(emx_perm))),by=.(rho=factor(rho),th_ldw)]

tmp3[rho=="0.99",rho:="1"]
p3 <- ggplot(tmp, aes(
  x = as.numeric(as.character(rho)),
  y = th_ldw,
  fill = `Number of outliers`
)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = wes_palette("Zissou1", 100, type = "continuous")
  ) +
  labs(
    x = expression("LD window size relative to LD-decay ("*rho*")"),
    y = expression("Quantile trehshold for local LD"),
    fill = "No. outliers"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.margin = margin(0, 0, 0, 0),
        panel.grid = element_blank(),
        panel.grid.major = element_blank())

(p1 | p2 | p3)

dt <- rbind(cbind(data="Observed data",tmp1),
      cbind(data="Permuted (within geographic region)", tmp2),
      cbind(data="Simulated (based on GRM)",tmp3))

ggplot(dt, aes(
  x = as.numeric(as.character(rho)),
  y = th_ldw,
  fill = `Number of outliers`
)) +
  facet_grid(.~data)+
  geom_tile() +
  scale_fill_gradientn(
    colors = wes_palette("Zissou1", 100, type = "continuous")
  ) +
  labs(
    x = expression("LD window size relative to LD-decay ("*rho*")"),
    y = expression("Quantile trehshold for local LD"),
    fill = "No. outliers"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.margin = margin(0, 0, 0, 0),
        panel.grid = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))


make_p <- function(d) {
  ggplot(d, aes(
    x = as.numeric(as.character(rho)),
    y = th_ldw,
    fill = `Number of outliers`
  )) +
    geom_tile() +
    scale_fill_gradientn(
      colors = wes_palette("Zissou1", 100, type = "continuous")
    ) +
    labs(
      x = NULL,#expression("LD window size relative to LD-decay ("*rho*")"),
      y = NULL,#expression("Quantile threshold for local LD"),
      fill = "No. outliers"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.margin = margin(1, 4, 1, ),
      panel.grid = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle(d[1,data])
}

p_list <- split(dt, dt$data) |> lapply(make_p)

wrap_plots(p_list, nrow = 1)
