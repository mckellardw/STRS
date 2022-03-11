

# Injury zone clusters ----
tmp.i <- c(1:3,5:7) # horizontal - ctrl & polyA
 # tmp.i <- c(1,5,2,6,3,7) # vertical - ctrl & polyA
tmp.i <- 4:7 # horizontal - all polyA 
cluster.map <- lapply(
  skm.list[tmp.i],
  FUN=function(VIS) DimPlot(
    VIS,
    reduction = "space",
    pt.size=0.4,
    cols = list(
      injury_zone=mckolors$primary[1],
      injury_border=mckolors$primary[2],
      myofiber=mckolors$primary[3]
    ),
    group.by = c(
      "injury.zones"
      # "kallisto_collapsed_snn_res.1"
    )
  )+
    ggtitle(
      VIS$sample[1]%>%
              stringr::str_remove("SkM-")%>%
              stringr::str_remove("Pro_")
      )+
    scTheme$umap+
    theme(
      axis.title = element_blank(),
      axis.line =  element_blank(),
      plot.title=element_text(color="black",hjust=0.5),
      plot.margin = unit(rep(0,4), units="in")
    )
)

cluster.map%>% 
  wrap_plots(
    guides="collect"
  )&coord_fixed(ratio=1/1.6)


# cluster vlns ----
tmp.feat <- c(
  ## Canonical
  # "Myog",
  # "Ctss",
  # "Acta1",
  
  ## ncRNAs!
  "Gm13568",
  "Notch1"
  # "H19",
  # "Meg3",
  # "Itm2a"
  
)
lapply(
  c(
    tmp.feat
    # "kal.protein_coding",
    # "kal.lncRNA",
    # "kal.Mt_tRNA",
    # "kal.miRNA"
    # "kal.snoRNA"
    ),
  FUN=function(Y) ggplot(
    # skm.merged@meta.data,
    cbind(
      skm.merged@meta.data,
      t(
        GetAssayData(
          skm.merged,
          assay="kallisto_collapsed",
          slot="data"
        )[tmp.feat,]
      )
    ),
    aes_string(
      x=as.character("timepoint"),
      group=as.character("timepoint"),
      y=Y
    )
  )+
    geom_violin(
      # fill=mckolors$txg[4]
      aes(
        fill=injury.zones
      )
    )+
    geom_jitter(
      color=gray(0.42),
      alpha=0.1,
      size=0.1
    )+
    scTheme$vln+
    theme(
      plot.margin = unit(rep(0,4), units="in"),
      axis.title.x=element_blank(),
      axis.title.y=element_text(face="italic"),
      axis.ticks = element_line(color="black"),
      axis.line=element_line(color="black")
    )+
    scale_fill_manual(
      values = list(
        injury_zone=mckolors$primary[1],
        injury_border=mckolors$primary[2],
        myofiber=mckolors$primary[3]
      )
    )+
    scale_x_continuous(
      breaks = c(0,2,5,7)
    )+
    guides(
      color = guide_legend(override.aes = list(size=2))
    )+
    facet_wrap(
      facets = "injury.zones"
    )
)%>%
  wrap_plots(
    guides="collect",
    ncol=1
    )



# Notch1/Gm13568 ----
notch.map <- visListPlot(
  # skm.list[meta_skm$polyA=="yPAP"],
  # sample.titles = c("Uninjured", "2dpi", "5dpi","7dpi"),
  skm.list,
  sample.titles = c(
    "ctrl 2dpi", "ctrl 5dpi","ctrl 7dpi",
    "polyA Uninjured", "polyA 2dpi", "polyA 5dpi","polyA 7dpi"
    ),
  reduction="space",
  assay="kallisto_collapsed", 
  slot = 'data',
  pt.size=0.3,
  legend.position = "bottom",
  font.size = small.font,
  axis.title.angle.y=0,
  nrow = 1,
  combine = T,
  verbose=F,
  colormap = "viridis",
  # colormap = mckolors$Spectral%>%rev(),
  features = c(
    "Notch1",
    "Gm13568"
    
    # "Notch3",
    # "Gm17276"
  )
)&theme(
  legend.text = element_text(size=small.font-2)
)&coord_fixed(
  ratio=1/1.6
)
notch.map

# scatter plot showing Notch1 vs. Gm13568 expression
notch.scatter <- lapply(
  skm.list[meta_skm$polyA=="yPAP"],
  FUN=function(VIS){
    tmp.nonzeros <- GetAssayData(
      VIS,
      assay="kallisto_collapsed"
    )[c("Notch1","Gm13568"),]%>%
      apply(
        MARGIN=2,
        FUN=function(X) data.frame(
          NOTCH = X[1]>0 & X[2]==0,
          GM = X[1]==0 & X[2]>0,
          DP = X[1]>0 & X[2]>0
        )
      )%>%
      do.call(what=rbind)
    
    tmp.label <- paste0(
      "Notch1+: ",table(tmp.nonzeros$NOTCH)["TRUE"],"/",sum(tmp.nonzeros),"\n",
      "Gm13568+:",table(tmp.nonzeros$GM)["TRUE"],"/",sum(tmp.nonzeros),"\n",
      "Double Pos: ",table(tmp.nonzeros$DP)["TRUE"],"/",sum(tmp.nonzeros)
    )
    
    FeatureScatter(
      VIS,
      shuffle = T,jitter = F,
      feature1 = "Notch1",
      feature2 = "Gm13568",
      plot.cor = F,
      pt.size = 1,
      cols=list(
        injury_zone=mckolors$primary[1],
        injury_border=mckolors$primary[2],
        myofiber=mckolors$primary[3]
      ),
      group.by = "injury.zones"
    )+
      xlim(c(0,2))+
      ylim(c(0,6))+
      # ggtitle(VIS$sample[1])+
      geom_label(
        x=1,
        y=5.2,
        size = small.font/ggplot2::.pt,
        # size = small.font/4,
        label=tmp.label
      )+
      guides(
        color = guide_legend(override.aes = list(size=2))
      )+
      scTheme$scatter+
      theme(
        axis.title.y = element_text(
          face="bold.italic",
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

notch.scatter[[length(notch.scatter)]] <- notch.scatter[[length(notch.scatter)]]+
  theme(
    axis.title.y = element_text(
      face="bold.italic",
      hjust=0.5,
      size=small.font
    )
  )

notch.scatter <- wrap_plots(
  notch.scatter,
  guides="collect",
  nrow=1
)&theme(
  legend.position = "bottom"
)
notch.scatter

# wrap_plots
wrap_plots(
  cluster.map%>% 
    wrap_plots(
      guides="collect",
      nrow=1
    )&coord_fixed(ratio=1/1.6)&theme(
      legend.position="none"
    ),
  notch.scatter,
  # nrow=1,
  # widths=c(2,1)
  nrow=2, heights=c(1,1.75),
  guides="collect"
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/Fig_SkM_Notch1_v1.pdf",
  device="pdf",
  units="cm",
  width = 7*2,
  height = 5*2
)


# mir1 ----
tmp.feat=c(
  "mmu-miR-1a-3p",
  "mmu-miR-206-3p"
  # "mmu-miR-133a-3p/133b-3p",
)%>%sort()

wrap_plots(
  visListPlot(
    skm.list[meta_skm$polyA=="yPAP"],
    sample.titles = meta_skm$sample[meta_skm$polyA=="yPAP"],
    alt.titles = tmp.feat,
    reduction = "space",
    assay="mirge3",
    slot="data",
    pt.size = 0.1,
    # features="nCount_mirge3",
    features=tmp.feat,
    axis.title.angle.y = 0,
    legend.position = "bottom",
    combine=T,
    colormap = "plasma",
    colormap.direction = -1,
    colormap.same.scale = F
  )&coord_fixed(
    ratio = 1/1.6
  ),
  
  visListPlot(
    skm.list[meta_skm$polyA=="yPAP"],
    # sample.titles = meta_skm$sample[meta_skm$polyA=="yPAP"],
    sample.titles=rep("",4),
    reduction = "space",
    assay="mirge3",
    slot="data",
    pt.size = 0.1,
    features=c(
      "mir1_3p_targets1",
      "mir1_3p_targets_spliced1"
    ),
    alt.titles = c(
      "miR1 Target Score",
      "miR1 Target Score (spliced)"
      # "Mir133a Target Score"
    ),
    axis.title.angle.y = 0,
    legend.position = "bottom",
    combine=T,
    colormap = mckolors$Spectral%>%rev(),
    colormap.same.scale = F
  )&coord_fixed(
    ratio = 1/1.6
  )&theme(
    axis.title.x=element_blank()
  ),
  
  nrow=1,
  widths = c(2,2)
)





mir1.scatter <- lapply(
  skm.list[meta_skm$polyA=="yPAP"],
  FUN=function(VIS){
    
    FeatureScatter(
      VIS,
      shuffle = T,jitter = F,
      feature1 = "mmu-miR-1a-3p",
      # feature1 = "mir1_3p_targets_spliced1",
      feature2 = "mir1_3p_targets1",
      plot.cor = T,
      pt.size = pt.size,
      cols=list(
        injury_zone=mckolors$primary[1],
        injury_border=mckolors$primary[2],
        myofiber=mckolors$primary[3]
      ),
      group.by = "injury.zones"
    )+
      xlim(c(0,4.2))+
      ylim(c(-0.15,0.25))+
      labs(
        x="miR-1a-3p",
        y="miR1 Targets Score"
      )+
      # ggtitle(VIS$sample[1])+
      scTheme$scatter+
      theme(
        axis.title = element_text(
          face="italic",
          hjust=0.5,
          size=small.font
        ),
        plot.margin = unit(rep(0,4),"cm")
      )+
      guides(
        color = guide_legend(override.aes = list(size=2))
        )%>%
      return()
  }
)%>%
  wrap_plots(
    guides="collect",
    # ncol=1
    nrow=1
  )&theme(
    legend.position = "bottom"
  )
mir1.scatter

# Spliced vs. total counts miR target scores
lapply(
  
)





















