
# Gene biotype spatial maps

## Skeletal muscle
i = c(
  1:3, # SkM ctrl
  # 4,6,7,5,8,9, # SkM SUPERase
  10:13 # # SkM Protector
)
visListPlot(
  vis.list[i],
  sample.titles = stringr::str_remove_all(meta_vis$sample[i],pattern = "Vis_") %>%
    stringr::str_remove_all(pattern ="yPAP_")%>%
    stringr::str_remove_all(pattern ="ctrl_")%>%
    stringr::str_remove_all(pattern ="_Heart")%>%
    stringr::str_remove_all(pattern ="_SkM"),
  # assay="STARsolo_collapsed",
  reduction="space",
  # slot = 'counts',
  pt.size=0.3,
  legend.position = "bottom",
  font.size = 8,
  axis.title.angle.y=0,
  nrow = 1,
  combine = T,
  verbose=F,
  colormap="plasma",
  # colormap = mckolors$RdYlBu %>% rev(),
  features = c(
    "pct.protein_coding",
    "pct.miRNA",
    "pct.rRNA",
    "pct.Mt_rRNA",
    "pct.lncRNA",
    "pct.misc_RNA",     
    "pct.Mt_tRNA",
    "pct.snoRNA",
    "pct.snRNA",   
    "pct.pseudogene",
    "pct.ribozyme",                
    "pct.scaRNA"
  )
)&coord_fixed(ratio=1/1.6)&scTheme$space&theme(legend.text = element_text(size=6))

## Heart
i = c(
  14:17
)
visListPlot(
  vis.list[i],
  sample.titles = stringr::str_remove_all(meta_vis$sample[i],pattern = "Vis_") %>%
    stringr::str_remove_all(pattern ="yPAP_")%>%
    stringr::str_remove_all(pattern ="ctrl_")%>%
    stringr::str_remove_all(pattern ="_Heart")%>%
    stringr::str_remove_all(pattern ="_SkM"),
  # assay="STARsolo_collapsed",
  reduction="space",
  # slot = 'counts',
  pt.size=0.3,
  legend.position = "bottom",
  font.size = 8,
  axis.title.angle.y=0,
  nrow = 1,
  combine = T,
  verbose=F,
  colormap="plasma",
  # colormap = mckolors$RdYlBu %>% rev(),
  features = c(
    "pct.protein_coding",
    "pct.miRNA",
    "pct.rRNA",
    "pct.Mt_rRNA",
    "pct.lncRNA",
    "pct.misc_RNA",     
    "pct.Mt_tRNA",
    "pct.snoRNA",
    "pct.snRNA",   
    "pct.pseudogene",
    "pct.ribozyme",                
    "pct.scaRNA"
  )
)&coord_fixed(ratio=1.6)&scTheme$space&theme(legend.text = element_text(size=6))


# xGene and spTotal efficiency
wrap_plots(
  ggplot( #std vs xGen
    heart.merged@meta.data[sample(Cells(heart.merged)),],
    aes(
      x=reo.log2p1,
      y=nCount_xGen.kallisto.log2p1,
      # x=reo.counts,
      # y=nCount_xGen.kallisto,
      color=polyA
    )
  )+
    geom_abline()+
    geom_point(
      alpha=0.7
    )+
    scTheme$scatter+
    scale_color_manual(
      values=mckolors$txg[c(1,4)]
    )+
    labs(
      x="Standard library prep\nlog2(Reovirus UMIs+1)",
      y="xGen Enrichment\nlog2(Reovirus UMIs+1)"
    ),
  plot_spacer(),
  ggplot( # sense vs. antisense
    heart.merged@meta.data[sample(Cells(heart.merged)),],
    aes(
      x=nCount_xGen.kallisto.sense.log2p1,
      y=nCount_xGen.kallisto.as.log2p1,
      color=polyA
    )
  )+
    geom_abline()+
    geom_point(
      alpha=0.7
    )+
    scTheme$scatter+
    scale_color_manual(
      values=mckolors$txg[c(1,4)]
    )+
    labs(
      x="Sense strand, all Reovirus segments\nlog2(UMIs+1)",
      y="Antisense strand, all Reovirus segments\nlog2(UMIs+1)"
    ),
  
  guides="collect",ncol = 1,
  heights=c(1,0.1,1)
)

