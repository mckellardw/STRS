
# Prep heart.list ----
heart.list <- lapply(
  heart.list,
  FUN=function(SEU){ 
    SEU$AnatomicalRegion <- factor(
      SEU$AnatomicalRegion,
      levels=c(
        "Atria",
        "Inflamed atria",
        "Ventricle",
        "Inflamed ventricle",
        "Border zone",
        "Myocarditic region",
        "Cavity"
      )
    )
    return(SEU)
  }
)
# Cluster maps ----
tmp.i <- 1:4 # horizontal - all polyA 
cluster.map <- lapply(
  heart.list[tmp.i],
  FUN=function(VIS) DimPlot(
    VIS,
    reduction = "space",
    pt.size=0.3,
    cols = list(
      Ventricle=mckolors$colblind_8[8],
      Atria=mckolors$colblind_8[7],
      `Inflamed ventricle`=mckolors$colblind_8[3],
      `Inflamed atria`=mckolors$colblind_8[4],
      `Myocarditic region`=mckolors$colblind_8[1],
      `Border zone`=mckolors$colblind_8[2],
      Cavity=mckolors$colblind_8[6]
    ),
    group.by = c(
      "AnatomicalRegion"
      # "kallisto_collapsed_snn_res.1"
    )
  )+
    ggtitle(
      VIS$sample[1]%>%
        stringr::str_remove("Heart-")%>%
        stringr::str_remove("Pro_")
    )+
    scTheme$umap+
    theme(
      axis.title = element_blank(),
      axis.line =  element_blank(),
      plot.title=element_text(size=small.font,color="black",hjust=0.5),
      plot.margin = unit(rep(0,4), units="in")
    )
)

cluster.map[[1]] <- cluster.map[[1]] + ggtitle("Standard - Mock")
cluster.map[[2]] <- cluster.map[[2]] + ggtitle("Standard - Infected")
cluster.map[[3]] <- cluster.map[[3]] + ggtitle("Standard - Mock")
cluster.map[[4]] <- cluster.map[[4]] + ggtitle("Standard - Infected")

cluster.map%>% 
  wrap_plots(
    nrow=1,
    guides="collect"
  )&coord_fixed(ratio=1.6)

# Correlation of gene expression w/ reo.counts ----
library(ggpmisc)
heart.list[[4]]@active.assay <- "kallisto_collapsed"
tmp.feat <- c(
  "Ccl2", "Cxcl9","Gzma","AW112010",
  "Fabp3", "Slc25a4", "Cox8b", "Fhl2"
)
tmp.feat <- c("Cxcl11")
correlation.scatter <- lapply(
  tmp.feat,
  FUN=function(FEAT){
    FeatureScatter(
      heart.list[[4]],
      cells = Cells(heart.list[[4]])[!heart.list[[4]]$AnatomicalRegion%in%c("Cavity","Atria","Inflamed atria")],
      shuffle = T,jitter = F,
      feature1 = "nCount_xGen.kallisto.log2p1",
      feature2 = FEAT,
      plot.cor = T,
      pt.size = 0.6,
      cols=list(
        Ventricle=mckolors$colblind_8[8],
        Atria=mckolors$colblind_8[7],
        `Inflamed ventricle`=mckolors$colblind_8[3],
        `Inflamed atria`=mckolors$colblind_8[4],
        `Myocarditic region`=mckolors$colblind_8[1],
        `Border zone`=mckolors$colblind_8[2],
        Cavity=mckolors$colblind_8[6]
      ),
      group.by = "AnatomicalRegion"
    )+
      guides(
        color = guide_legend(override.aes = list(size=2))
      )+
      geom_smooth(
        # formula='y~x',
        method="gam",
        color="black"
      )+
      # stat_poly_eq(
      #   method = "lm",
      #   aes(
      #     label = paste(..eq.label.., ..rr.label.., sep = "~~~")
      #   ), 
      #   parse = TRUE
      # ) + 
      scTheme$scatter+
      theme(
        axis.title.y = element_text(
          face="italic",
          hjust=0.5,
          size=small.font
        ),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(rep(0,4),"cm")
      )%>%
      return()
  }
)

correlation.scatter <- wrap_plots(
  correlation.scatter,
  guides="collect",
  nrow=1
)&theme(
  legend.position = "bottom"
)#&coord_fixed(ratio=2)
correlation.scatter
#
# Spatial GE plots ----
tmp.feat <- c(
  "AA474408",
  "AW112010",
  "2410006H16Rik",
  # "Bc1",
  "Ly6a2",
  "Cxcl11", 
  "Mx2"
)
heart.feat.maps <- visListPlot(
  heart.list,
  sample.titles = c(
    "Visium\nControl","Visium\nInfected",
    "STRS\nControl", "STRS\nInfected"
    ),
  reduction = "space",
  assay="kallisto_collapsed",
  slot="data",
  pt.size = 0.3,
  features=tmp.feat,
  axis.title.angle.y = 0,
  legend.position = "bottom",
  combine=T,
  nrow = 1,
  colormap = "viridis",
  colormap.direction = 1,
  colormap.same.scale = F
)&coord_fixed(
  ratio = 1.6
)
heart.feat.maps

# wrap and save ----
wrap_plots(
  wrap_plots(
    plot_spacer(),
    cluster.map%>% 
      wrap_plots(
        ncol=1,
        guides="collect"
      )&coord_fixed(ratio=1.6)&NoLegend()&theme(plot.title=element_blank()),
    heart.feat.maps,
    nrow=1,
    widths=c(1,1.1,7)
  ),
  wrap_plots(
    correlation.scatter,
    plot_spacer(),
    nrow=1,
    widths = c(1,0.6)
  ),
  heights=c(4.5,2),
  ncol=1
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/Fig_Heart_v2.pdf",
  device="pdf",
  units="cm",
  width = 20*2,
  height = 20*2
)
