# Comparison of STAR and kallisto for quantifying total RNA data

ks.scatter <- wrap_plots(
  ggplot(
    vis.merged@meta.data[sample(rownames(vis.merged@meta.data)),],
    aes(
      x=nCount_STARsolo,
      y=nCount_kallisto,
      color=rnase_inhib
    )
  )+
    geom_abline()+
    geom_point(
      alpha=0.6,size=1
    )+
    scale_color_manual(values=mckolors$txg[c(1,4,2)])+
    scale_x_continuous(labels=scales::scientific)+
    scale_y_continuous(labels=scales::scientific)+
    scTheme$scatter,
  ggplot(
    vis.merged@meta.data[sample(rownames(vis.merged@meta.data)),],
    aes(
      x=nFeature_STARsolo,
      y=nFeature_kallisto,
      color=rnase_inhib
    )
  )+
    geom_abline()+
    geom_point(
      alpha=0.6,size=1
    )+
    scale_color_manual(values=mckolors$txg[c(1,4,2)])+
    scale_x_continuous(labels=scales::scientific)+
    scale_y_continuous(labels=scales::scientific)+
    scTheme$scatter,
  guides="collect"
)

# Maps of feature counts & UMI counts compared between star & kallisto ----
wrap_plots(
  visListPlot(
    skm.list,
    sample.titles = stringr::str_remove_all(meta_skm$sample,pattern = "Vis_") %>%
      stringr::str_remove_all(pattern ="_SkM"),
    reduction="space",
    pt.size=0.4,
    legend.position = "bottom",
    font.size = 6,
    axis.title.angle.y=0,
    nrow = 1,
    combine = T,
    verbose=F,
    colormap = mckolors$Spectral %>% rev(),
    colormap.same.scale = T,
    features = c(
      "nFeature_STARsolo_collapsed","nFeature_kallisto_collapsed"
    )
  )&theme(
    legend.text = element_text(size=small.font)
  )&coord_fixed(ratio=1/1.6),
  
  visListPlot(
    skm.list,
    sample.titles = stringr::str_remove_all(meta_skm$sample,pattern = "Vis_") %>%
      stringr::str_remove_all(pattern ="_SkM"),
    reduction="space",
    pt.size=0.4,
    legend.position = "bottom",
    font.size = 6,
    axis.title.angle.y=0,
    nrow = 1,
    combine = T,
    verbose=F,
    colormap = mckolors$Spectral %>% rev(),
    colormap.same.scale = T,
    features = c(
      "nCount_STARsolo_collapsed","nCount_kallisto_collapsed"
    )
  )&theme(
    legend.text = element_text(size=small.font),
    axis.title.y = element_blank()
  )&coord_fixed(ratio=1/1.6),
  nrow=1
)

# wrap and save plots ----
wrap_plots(
  ks.scatter
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_kallisto_star.pdf",
  device="pdf",
  units="cm",
  width = 24*2,
  height = 10*2
)
