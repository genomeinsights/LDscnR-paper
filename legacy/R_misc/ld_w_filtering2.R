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

  if(th_ldw>0){
    keep <- ld_sub[, rho] > quantile(ld_sub[, rho], th_ldw, na.rm = TRUE)
    keep[is.na(keep)] <- FALSE
  }else{
    keep <- rep(TRUE,nrow(ld_sub))
  }
  #table(keep)

  if (!any(keep)) {
    return(cbind(empty_result(rho, th_ldw, p_names, r2_grid, lmin_grid),n_loci=0))
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
ld_decay_3sp$by_chr$Chr1$decay[,agg]

dt <- rbindlist(lapply(ld_decay_3sp$by_chr, function(ch){
  rbindlist(ch$decay[1/a<500,agg])[d_mid<10000,.(d_mid,r2_q)]
}))

ch <- ld_decay_3sp$by_chr$Chr1
rbindlist(ch$decay[1/a<500,agg])[d_mid<10000,plot(d_mid,r2_q)]
b <- ch$decay$b[1]
c <- ch$decay$c[1]
a <- ch$decay$a[1]

dt <- ch$decay[3,raw][[1]][d<50000,.(b,c,a,d,r2)]
dt[,r2_pred := b + (c - b)/(c + a * d)]

p_decay <- ggplot(dt[d<10000], aes(d,r2)) +
  geom_point(fill="black",shape=20,alpha=0.5,size=5) +
  theme_bw(base_size = 18) +
  geom_line(aes(d,r2_pred),col="#E31A1C",linewidth = 2)+
  labs(x=expression(bold("Distance (bps)")),
       y=expression(bold(r^2)))+
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  geom_vline(xintercept = dt[1,d_from_rho(a,rho = 0.95)],col="#1F78B4",linetype=2,linewidth=1.5) +
  annotate("text",x=6000,y=0.85,label=expression(bold(w[rho==0.95]==2900~bps)),col="#1F78B4",size=12)
p_decay

#b + (1 - b)/(1 + a * d)
# +
  #geom_hline(yintercept = dt[1,b],col="#E31A1C",linetype=2,linewidth=2)#+
  #geom_hline(yintercept = dt[1,ld_from_rho(b = b,c = c,rho = 0.95)],col="steelblue",linetype=2,linewidth=1.5)

ggsave(p_decay,filename="p_decay.png",width = 10,height = 6,units = "in",dpi = 300)


rbindlist(ld_decay_3sp$by_chr$Chr5$decay[1/a<500,agg])[d_mid<10000,plot(d_mid,r2_q)]
rbindlist(ld_decay_3sp$by_chr$Chr5$decay[1/a<500,agg])[d_mid<10000,plot(d_mid,r2_q)]
ld_decay_3sp$by_chr$Chr1$decay[,agg][[1]][,plot(d_mid,r2_q)]

ld_ws_3sp <- precalculate_ld_w(seq(0,0.95,by=0.05),ld_decay_3sp)
rownames(ld_ws_3sp) <- map_3sp$marker
colnames(ld_ws_3sp)
ld_ws_3sp[is.na(ld_ws_3sp)] <- 0
ld_ws_3sp[,1]<-sample(ld_ws_3sp[,"0.95"])

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
  th_ldw = th_ldw_grid,
  rho = colnames(ld_ws_3sp)
)

#i <- 10
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
## heatmap --------------------------------------------------
tmp1 <- outliers_100_perm[l_min==10 & r2_th==0.6,.("Number of outliers"=mean(lengths(EMX))),by=.(rho=factor(rho),th_ldw)]
tmp1.2 <- outliers_100_perm[l_min==10 & r2_th==0.6,.("Number of outliers"=mean(lengths(LFMM))),by=.(rho=factor(rho),th_ldw)]

tmp1.2[,data:="LFMM"]
tmp1[,data:="EMMAX"]

dt <- rbind(tmp1,tmp1.2)
dt <- rbind(tmp1.2)

dt[, th_ldw_f := factor(
  round(th_ldw,3),
  levels = round(th_ldw_grid,3),
  ordered = TRUE
)]

dt[, rho_f := factor(
  rho,
  levels = unique(rho),
  ordered = TRUE
)]

p1 <- ggplot(dt, aes(
  x = rho_f,
  y = th_ldw_f,
  fill = `Number of outliers`
)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(. ~ data) +
  scale_fill_gradientn(
    colors = wes_palette("Zissou1", 100, type = "continuous"),
    name = "Outlier \nregions"
  ) +
  scale_x_discrete(
    breaks = factor(seq(0,0.9,by=0.1)),
    labels = as.character(seq(0,0.9,by=0.1))
  ) +
  scale_y_discrete(
    breaks = round(th_ldw_grid,3)[c(1,4,8,12,16,20)],
    labels = paste(as.character(round(th_ldw_grid,2)[c(1,4,8,12,16,20)]*100),"%")
  ) +
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
    plot.margin = margin(1, 1, 1, 1)
  )
p1
ggsave(p1,filename="the_heatmap.png",height = 6,width = 7,units = "in",dpi=300)

## Manhattan plots -------------------------------
col_vector <- c("#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
                "#FB8072", "#E6AB02", "steelblue",
                "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
                "#A6CEE3", "salmon",  "forestgreen", "#BF5B17",
                "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#BEBADA",
                 "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
                "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C")

map_C_perm <- summarise_stability(
  outliers = outliers_100_perm[l_min==10 & r2_th==0.6],
  map = map_3sp,
  p_names = names(p_cols)
)
outl_markers <- map_C_perm[C_LFMM > 0.05,marker]
el <- precompute_LD_edges(
  GTs = GTs_3sp[, outl_markers, drop = FALSE],
  map = map_3sp[marker %in% outl_markers],
  r2_min = 0.1,
  max_bp = 1e6,
  cores = cores
)
outl_cls <- LD_igraph_components(
  el = el,
  markers = outl_markers,
  r2_th = 0.4,
  bp_th = Inf
)

map_3sp_manh <- outl_cls[n_loci>=10,.(marker,CL_id)][map_3sp,on="marker"]
map_3sp_manh[,CL_id:=as.character(CL_id)]
map_3sp_manh[is.na(CL_id),CL_id:="ns"]
length(map_3sp_manh[,table(CL_id)])



p_manh <- ggplot(map_3sp_manh, aes(Pos, -log10(lfmm_q))) +
  geom_point(
    data = map_3sp_manh[CL_id == "ns"],
    color = "grey35", size = 0.7, alpha = 0.35
  ) +
  geom_point(
    data = map_3sp_manh[CL_id != "ns"],
    aes(color = CL_id),
    size = 2.2, alpha = 1
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    linewidth = 1
  ) +
  facet_wrap(~ Chr, nrow = 1, scales = "free_x") +
  scale_color_manual(values = rep(col_vector, 10)) +
  labs(
    x = "Genomic position",
    y = expression(-log[10](italic(q))~"(LFMM)")
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(3, 3, 3, 3),
    panel.spacing.x = unit(0.05, "lines")
  )
ggsave(p_manh,filename="the_manhattan.png",height = 5,width = 18,units = "in",dpi=600)

map_C_perm[,plot(LFMM)]
map_C_perm[,plot(C_EMX)]
map_C_perm[C_LFMM > 0.1 | C_EMX>0.1,marker]

## -----------
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
