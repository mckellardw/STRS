
# Gene biotype spatial maps

## Skeletal muscle
i = c(
  1:3, # SkM ctrl
  # 4,6,7,5,8,9, # SkM SUPERase
  10:13 # # SkM Protector
)
suppressMessages(
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
    pt.size=0.4,
    legend.position = "bottom",
    font.size = small.font,
    axis.title.angle.y=0,
    nrow = 1,
    combine = T,
    verbose=F,
    colormap="plasma",
    # colormap = mckolors$RdYlBu %>% rev(),
    features = c(
      "kal.protein_coding",
      "kal.rRNA",
      "kal.Mt_rRNA",
      "kal.miRNA",
      "kal.lncRNA",
      "kal.Mt_tRNA",
      "kal.snoRNA",
      "kal.snRNA",   
      "kal.ribozyme",  
      "kal.misc_RNA",              
      "kal.scaRNA"
    )
  )&coord_fixed(ratio=1/1.6)&
    scTheme$space&
    theme(
      legend.text = element_text(size=6*2)
    )
    # scale_color_manual(labels=scales::percent)
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_biotypes_SkM.pdf",
  device="pdf",
  units="cm",
  width = 24*2,
  height = 12*2
)

## Heart ----
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
  pt.size=0.35,
  legend.position = "bottom",
  font.size = small.font,
  axis.title.angle.y=0,
  nrow = 1,
  combine = T,
  verbose=F,
  # colormap.same.scale = T,
  colormap="plasma",
  # colormap = mckolors$RdYlBu %>% rev(),
  features = c(
    "kal.protein_coding",
    "kal.rRNA",
    "kal.Mt_rRNA",
    "kal.miRNA",
    "kal.lncRNA",
    "kal.Mt_tRNA",
    "kal.snoRNA",
    "kal.snRNA",   
    "kal.ribozyme",  
    "kal.misc_RNA",              
    "kal.scaRNA"
  )
)&coord_fixed(ratio=1.6)&
  scTheme$space&
  theme(
    legend.text = element_text(size=6*2)
  )

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_biotypes_heart.pdf",
  device="pdf",
  units="cm",
  width = 24*2,
  height = 11*2
)


# xGene and spTotal efficiency ----
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

