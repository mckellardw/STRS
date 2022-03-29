# Comparison of STAR and kallisto for quantifying total RNA data


ggplot(
  vis.merged@meta.data,
  aes(
    x=kal.protein_coding,
    y=Ss_protein_coding
  )
)+
  geom_point()+
  scTheme$scatter


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

