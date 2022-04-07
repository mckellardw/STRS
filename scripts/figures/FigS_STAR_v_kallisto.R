# Comparison of STAR and kallisto for quantifying total RNA data


# UMI & feature count comparison scatter plots ----
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


# intergenic counts ----
# prep meta_data
meta_vis <- read.csv("/workdir/dwm269/totalRNA/spTotal/resources/metadata_sheets/meta_sheet_visium.csv")

tmp.df <- lapply(
  meta_vis$data.dir.STARsolo,
  FUN=function(DIR) read.csv(
    file = paste0(DIR,"/Solo.out/Gene/Summary.csv"),
    header = F,
    row.names = 1
  ) %>% t()
) %>%
  do.call(what=rbind)
colnames(tmp.df) <- stringr::str_replace_all(colnames(tmp.df),pattern=" ",replacement = "_")
colnames(tmp.df) <- stringr::str_replace_all(colnames(tmp.df),pattern="\\+",replacement = "_")
colnames(tmp.df) <- stringr::str_remove_all(colnames(tmp.df),pattern="\\:")
meta_vis <- cbind(meta_vis, tmp.df)

# plot!
unique.scatter <- ggplot(
  meta_vis,
  aes_string(
    x="Reads_Mapped_to_Genome_Unique",
    y="Reads_Mapped_to_Gene_Unique_Gene",
    color="rnase_inhib",
    shape="tissue"
  )
)+
  geom_point(
    alpha=0.8,
    size=4
  )+
  labs(
    color="RNase Inhibitor/\nChemistry",
    shape="Tissue",
    x="Reads Mapped to Genome\n(Unique Only)",
    y="Reads Mapped to Annotated Genes\n(Unique Only)"
  )+
  scale_color_manual(values=mckolors$txg[c(1,4,2)])+
  scale_x_continuous(
    labels=scales::percent
    # limits = c(0.9,1)
  )+
  scale_y_continuous(
    labels=scales::percent
    # limits=c(NA,1)
  )+
  scTheme$scatter
unique.scatter

# Include multimappers
multi.scatter <- ggplot(
  meta_vis,
  aes_string(
    x="Reads_Mapped_to_Genome_Unique_Multiple",
    y="Reads_Mapped_to_Gene_Unique_Multipe_Gene",
    color="rnase_inhib",
    shape="tissue"
  )
)+
  geom_point(
    alpha=0.8,
    size=4
  )+
  labs(
    color="RNase Inhibitor/\nChemistry",
    shape="Tissue",
    x="Reads Mapped to Genome\n(Unique + Multimappers)",
    y="Reads Mapped to Annotated Genes\n(Unique + Multimappers)"
  )+
  scale_color_manual(values=mckolors$txg[c(1,4,2)])+
  scale_x_continuous(
    labels=scales::percent,
    limits = c(0.9,1)
  )+
  scale_y_continuous(
    labels=scales::percent,
    limits=c(NA,1)
  )+
  scTheme$scatter


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
