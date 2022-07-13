# Rarefaction analysis of Visium and STRS libraries ----
# Load data ----
# meta_down,down.df

# Set themes for plotting -----
source("/home/dwm269/DWM_utils/plotting_utils/scThemes.R")
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

# Plots ----

wrap_plots(
  ggplot( # nCounts per nReads
    cbind(meta_down,down.df),
    aes(
      x = nReads,
      y = nUMIs,
      group = sample,
      fill=chemistry,
      color=chemistry
    )
  )+
    geom_line(size=1) +
    geom_point(
      aes(
        size=nSpots,
        shape=tissue
      ),
      color="black",
      alpha=0.4
    )+
    scale_shape_manual(values=c(21,24))+
    # scale_size_continuous(breaks=unique(meta_down$nReads)%>%sort())+
    scale_y_continuous(limits=c(0,NA))+
    scale_color_manual(values=mckolors$txg[c(1,4)])+
    scale_fill_manual(values=mckolors$txg[c(1,4)])+
    scTheme$scatter,
  
  # ggplot( # nFeatures per nReads
  #   cbind(meta_down,down.df),
  #   aes(
  #     x = nReads,
  #     y = nFeatures,
  #     group = sample,
  #     fill=chemistry,
  #     color=chemistry
  #   )
  # )+
  #   geom_line(size=1) +
  #   geom_point(
  #     aes(
  #       size=nReads,
  #       shape=tissue
  #     ),
  #     color="black",
  #     alpha=0.4
  #   )+
  #   scale_shape_manual(values=c(21,24))+
  #   scale_size_continuous(breaks=unique(meta_down$nReads)%>%sort())+
  #   scale_y_continuous(limits=c(0,NA))+
  #   scale_color_manual(values=mckolors$txg[c(1,4)])+
  #   scale_fill_manual(values=mckolors$txg[c(1,4)])+
  #   scTheme$scatter,
  
  ggplot(
    cbind(meta_down,down.df),
    aes(
      x = nReads,
      y = nUMIs.unspliced,
      group = sample,
      fill=chemistry,
      color=chemistry
    )
  )+
    geom_line(size=1) +
    geom_point(
      aes(
        size=nSpots,
        shape=tissue
      ),
      color="black",
      alpha=0.4
    )+
    scale_shape_manual(values=c(21,24))+
    # scale_size_continuous(breaks=unique(meta_down$nReads)%>%sort())+
    scale_y_continuous(limits=c(0,NA))+
    scale_color_manual(values=mckolors$txg[c(1,4)])+
    scale_fill_manual(values=mckolors$txg[c(1,4)])+
    scTheme$scatter,
  
  # ggplot( # 10x flavor of saturation
  #   cbind(meta_down,down.df),
  #   aes(
  #     x = nReads,
  #     y = 1-(nUMIs/nReads),
  #     group = sample,
  #     fill=chemistry,
  #     color=chemistry
  #   )
  # )+
  #   geom_line(size=1) +
  #   geom_point(
  #     aes(
  #       size=nReads,
  #       shape=tissue
  #     ),
  #     color="black",
  #     alpha=0.4
  #   )+
  #   scale_shape_manual(values=c(21,24))+
  #   scale_size_continuous(breaks=unique(meta_down$nReads)%>%sort())+
  #   scale_y_continuous(limits=c(0,NA))+
  #   scale_color_manual(values=mckolors$txg[c(1,4)])+
  #   scale_fill_manual(values=mckolors$txg[c(1,4)])+
  #   scTheme$scatter,
  nrow=1,
  guides="collect"
)&
  labs(
    shape="Tissue",
    color="Method",
    fill="Method"
  )

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_rarefaction_v2.pdf",
  device="pdf",
  units="cm",
  width = 16*2,
  height = 7*2
)
