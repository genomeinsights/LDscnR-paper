## Figure for the c151 novel-vs-corroborated test (R/19). Both classifications shown,
## because neither is clean: c151 is nested inside c155, so classifying by c155 is
## deflationary for the novel class and classifying by c150 alone is inflationary.
suppressMessages({ library(data.table); library(ggplot2) })
P <- path.expand("~/gitlab/LDscnR-paper/kingman2021")
d <- rbindlist(lapply(c("both","c150"), function(k)
  fread(file.path(P,"data",sprintf("novel_regions_c151_%s.csv",k)))), fill=TRUE)
d[, engine := ifelse(grepl("EMMAX", group), "EMMAX", "LFMM")]
d[, cls := ifelse(grepl("NOVEL", group), "novel (no EcoPeak)", "corroborated")]
d[, panel := factor(ifelse(peakset=="both", "classified by c155 + c150\n(deflationary for novel)",
                                            "classified by c150 only\n(inflationary for novel)"))]
d[, lab := sprintf("%.2fx\np=%.3f\n%d reg / %.1f Mb", fold, p, n_regions, span_Mb)]
g <- ggplot(d, aes(cls, fold, fill = engine)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "grey30", linewidth = 0.25) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey30") +
  geom_text(aes(label = lab), position = position_dodge(0.8), vjust = -0.15, size = 2.5, lineheight = 0.95) +
  facet_wrap(~ panel) +
  scale_y_log10(expand = expansion(mult = c(0, 0.30))) +
  scale_fill_manual(values = c(EMMAX = "#E1AF00", LFMM = "#3B9AB2"), name = NULL) +
  labs(x = NULL, y = "fold enrichment of the top 1% of c151 p\n(vs within-chromosome rotation null, log scale)",
       title = "Do the regions that overlap no EcoPeak carry signal in a geography-matched cohort?",
       subtitle = "Corroborated regions are the internal positive control: c151 has the power. LFMM's novel regions sit exactly on the null under both classifications.") +
  theme_minimal(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
        panel.grid.major.x = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 7.5))
ggsave(file.path(P,"figures","fig9_novel_c151.png"), g, width = 9.5, height = 5, dpi = 200)
cat("wrote figures/fig9_novel_c151.png\n"); print(d[, .(peakset, group, n_regions, span_Mb, fold, p)])
