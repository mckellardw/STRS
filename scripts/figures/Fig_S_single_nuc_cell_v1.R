

# plotting settings ----

# fonts, sizes, etc.
small.font = 6*2
big.font = 8*2
line.width = 0.5*2
pt.size=0.01*2
pt.stroke=0.3*2
label.size=2*2

scTheme <- scThemes(
  small.font = small.font,
  big.font = big.font,
  line.width = line.width,
  pt.size=pt.size,
  pt.stroke=pt.stroke,
  label.size=label.size
)

# Feature violins ----
tmp.features <- c(
  # poly(A)+
  "Malat1",
  "Neat1",
  "Srsf11",
  "Paxbp1",
  
  # poly(A)-
  "Gm15564",
  "AY036118",
  "H1f3",
  "Gm37357",
  "Gm19220",
  "H1f2",
  "Gm42826",
  "Snord118",
  "Rpph1",
  "Scarna2",
  "Snord17",
  "7SK"
)

c2c.vln <- VlnPlot(
  nuc.merged,
  group.by="polyA",
  features=tmp.features,
  pt.size=0,
  combine=F
) %>% lapply(
  FUN = function(X) X +
    scTheme$vln +
    theme(
      plot.title=element_text(face="bold.italic", hjust=0.5,vjust=0),
      plot.margin = unit(rep(0,4),"cm"),
      legend.position = "none",
      axis.line.y = element_line(color="black"),
      axis.ticks.y = element_line(color="black"),
      axis.title.y=element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    scale_fill_manual(values=mckolors$txg[c(1,4)]) +
    scale_color_manual(values=mckolors$txg[c(1,4)])
) %>% wrap_plots(
  nrow=4,
  guides="collect"
)

# uTAR ----
nuc.merged$kallisto_vs_TAR_UMI <- nuc.merged$nCount_kallisto_collapsed / nuc.merged$nCount_TAR
nuc.merged$kallisto_vs_TAR_Feature <- nuc.merged$nFeature_kallisto_collapsed / nuc.merged$nFeature_TAR

tmp.feat <- c(
  # "nFeature_kallisto_collapsed",
  # "nCount_kallisto",
  # "nFeature_TAR",
  "nCount_TAR",
  "percent.uTAR"
  # "kal.miRNA"
)

utar.vln <- VlnPlot(
    nuc.merged,
    group.by="polyA",
    features=tmp.feat,
    cols=mckolors$txg[c(1,4)],
    # ncol = 1,
    pt.size = 0,
    combine = F
) %>%
  lapply(
    FUN = function(X) X +
      scTheme$vln +
      theme(
        legend.position="none",
        axis.line.y=element_line(color="black"),
        axis.title.y=element_blank(),
        axis.ticks.y=element_line(color="black"),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        plot.title=element_text(face="bold",hjust=0.5, size=big.font)
      )
  )

utar.vln[[1]] <- utar.vln[[1]]+
    scale_y_log10(labels=scales::scientific)

utar.vln[[2]] <- utar.vln[[2]]+
    scale_y_continuous(labels=scales::percent)

utar.vln <- utar.vln %>%
  wrap_plots(
    nrow=1
  )
utar.vln

# genomics origins pie ----

df.list <- 
  list(
    data.frame(
      sample = meta_nuc$sample[1],
      exonic =  c(50.14),
      intronic = c(43.05),
      intergenic = c(6.81),
      overlapping_exon = c(1.8)
    ),
    data.frame(
      sample = meta_nuc$sample[2],
      exonic =  c(35.11),
      intronic = c(24.91),
      intergenic = c(39.97),
      overlapping_exon = c(10.17)
    )
)

df.list <- lapply(df.list, reshape2::melt)
# colnames(df) <- c("Sample", "Locus", "pct.reads")

df.list <- lapply(
  df.list,
  FUN =function(df) df %>% mutate(
    csum = rev(cumsum(rev(value))), 
    pos = value/2 + lead(csum, 1),
    pos = if_else(is.na(pos), value/2, pos)
  )
)

origin.pie<-lapply(
  df.list,
  FUN = function(df) ggplot(
    df,
    aes(
      x=2,
      y=value,
      fill=variable
    )
  )+
    geom_col(
      width=1,
      color="white"
    )+
    ggrepel::geom_label_repel(
      aes(
        y = pos,
        label = paste0(round(value,digits = 2), "%"),
        color=Biotype
      ),
      color="white",
      max.overlaps = 100
    ) +
    coord_polar("y", start=0)+
    xlim(.2,2.5)+
    scTheme$pie+
    theme(
      legend.position = "right",
      plot.title=element_text(size=big.font,face="bold"),
      legend.title = element_text(size=small.font,face="bold")
    )+
    labs(fill="Genomic Locus")+
    scale_fill_manual(values = mckolors$rickandmorty_ggsci[c(3,5,2,4)])+
    facet_wrap(facets="sample")
) %>% wrap_plots(guides="collect")

# gene position ----
qualimap.list <- list()

for(i in 1:nrow(meta_nuc)){
  qualimap.list[[i]] <- fread(
    paste0(meta_nuc$data.dir.STARsolo[i],"/qualimap_out/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt")
  )
  colnames(qualimap.list[[i]]) <- c("position", "coverage")
  qualimap.list[[i]]$sample <- rep(meta_nuc$sample[i],nrow(qualimap.list[[i]]))
  qualimap.list[[i]]$type <- rep(meta_nuc$type[i],nrow(qualimap.list[[i]]))
  qualimap.list[[i]]$polyA <- rep(meta_nuc$polyA[i],nrow(qualimap.list[[i]]))
}
qualimap.df <- do.call(rbind,qualimap.list)

position.line <- ggplot(
  qualimap.df,
  aes(
    x=position/100,
    y=coverage,
    group=sample,
    # linetype=type,
    color=polyA
  )
) + 
  geom_line(size=1)+
  scTheme$scatter+
  theme(
    axis.line = element_line(color="black",size = line.width),
    axis.ticks = element_line(color="black",size = line.width),
    axis.text = element_text(color="black")
  )+
  labs(
    x="Normalized Gene Position",
    y="Coverage",
    color=""
  )+
  scale_color_manual(values=mckolors$txg[c(1,4)])+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::scientific)
position.line


# Wrap & save ----
wrap_plots(
  c2c.vln,
  plot_spacer(),
  wrap_plots(
    position.line,
    origin.pie,
    utar.vln,
    ncol=1,
    heights = c(1.5,1.5,1)
  ),
  nrow=1,
  widths=c(4,0.1,2)
)

ggsave(
  file = "spTotal/figures/single_nuc_v2.pdf",
  device="pdf",
  width = 18*2, 
  height = 8*2,
  units = "cm", 
  dpi=400
)
