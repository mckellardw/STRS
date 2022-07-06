# Spot deconvolution using scMuscle + BayesPrism

# Load data ----
if(!exists("skm.list")){
  load("/workdir/dwm269/totalRNA/spTotal/robjs/skm_vis_BP_scMuscle_v2.RData")
}

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

# BayesPrism maps ----
tmp.titles = stringr::str_remove_all(meta_skm$sample,pattern = "Vis_") %>%
  stringr::str_remove_all(pattern = "_TA") %>%
  stringr::str_remove_all(pattern = "-SkM")

tmp.titles <- c(
  "SkM\n(2dpi)\nVisium",
  "SkM\n(5dpi)\nVisium",
  "SkM\n(7dpi)\nVisium",
  "SkM\n(uninjured)\nSTRS",
  "SkM\n(2dpi)\nSTRS",
  "SkM\n(5dpi)\nSTRS",
  "SkM\n(7dpi)\nSTRS"
)

bp.map <- visListPlot(
  seu.list = skm.list,
  assay = "celltype.bp",
  sample.titles = tmp.titles,
  pt.size = 0.1,
  nrow = 1,
  min.value = 10^-2,
  font.size = small.font,
  axis.title.angle.y=0,
  colormap.same.scale = T,
  colormap = rev(RColorBrewer::brewer.pal(n=11,name="Spectral")),
  # features=GetAssayData(skm.list[[1]],assay="celltype.bp")%>% rownames()%>%unlist()
  features = c(
    "FAPs-(Pro-remodeling)",
    "Monocyte-(Patrolling)",
    "Quiescent-MuSCs",
    "Activated-MuSCs",
    "Committed-Myoblasts",
    "Fusing-Myocytes",
    "Myonuclei-(Type-IIb)",
    "Myonuclei-(Type-IIx)"
  ),
  alt.titles = c(
    "FAPs\n(Pro-remodeling)",
    "Monocyte\n(Patrolling)",
    "Quiescent\nMuSCs",
    "Activated\nMuSCs",
    "Committed\nMyoblasts",
    "Fusing\nMyocytes",
    "Myonuclei\n(Type-IIb)",
    "Myonuclei\n(Type-IIx)"
  )
)&coord_fixed(1/1.6)
  # scale_color_continuous(
  #   labels=scales::percent
  # )

# PCA plots ----
pca.list <- lapply(
  c("injury.zones","timepoint","polyA"),
  FUN=function(GROUP) DimPlot(
    strs.merged,
    cells=sample(Cells(strs.merged)),
    reduction = "pca",
    group.by=GROUP
  )+
    scTheme$umap+
    theme(
      plot.title=element_text(size=small.font, face="bold"),
      legend.position = "bottom",
      legend.direction = "vertical",
      plot.margin = unit(rep(0,4),"cm")
      # legend.text.align = "center"
    )+
    coord_fixed()+
    labs(
      x=paste0("BayesPrism_PC_1"),
      y=paste0("BayesPrism_PC_2")
    )
)

pca.list[[1]] <- pca.list[[1]] + 
  scale_color_manual(values=mckolors$primary[c(2,1,3)]) +
  ggtitle("Injury Zones")
pca.list[[2]] <- pca.list[[2]] + 
  scale_color_manual(values=viridis(n=4,option = "inferno"))+
  ggtitle("Injury Timepoint")
pca.list[[3]] <- pca.list[[3]] + 
  scale_color_manual(values=mckolors$txg[c(1,4)])+
  ggtitle("Method")

wrap_plots(
  pca.list,
  nrow=1
)


# Wrap final figure ----

wrap_plots(
  bp.map,
  plot_spacer(),
  wrap_plots(
    pca.list,
    nrow=1
  ),
  ncol=1,
  heights = c(7,0.1,3.25)
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/Fig_S_deconvolution_v1.pdf",
  device="pdf",
  units="cm",
  width = 15*2,
  height = 14*2
)
