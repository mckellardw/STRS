
if(!exists("skm.list")){
  load("robjs/skm_list_v2.RData")
}

if(!exists("skm.merged")){
  skm.merged <- merge(
    skm.list[[1]],
    skm.list[2:length(skm.list)],
    add.cell.ids = meta_skm$sample
  )
}

# Injury zone cluster maps ----
# tmp.i <- c(1:3,5:7) # horizontal - ctrl & polyA
 # tmp.i <- c(1,5,2,6,3,7) # vertical - ctrl & polyA
tmp.i <- 4:7 # horizontal - all polyA 
cluster.map <- lapply(
  skm.list[tmp.i],
  FUN=function(VIS) DimPlot(
    VIS,
    reduction = "space",
    pt.size=0.6,
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
      plot.title=element_text(color="black",hjust=0.5, size=small.font),
      plot.margin = unit(rep(0,4), units="in")
    )
)

cluster.map%>% 
  wrap_plots(
    guides="collect"
  )&coord_fixed(ratio=1/1.6)


# injury line plots ----
tmp.feat <- c(
  ## Canonical
  "Myog",
  # "S100a4",
  # "Ctss",
  # "Acta1",
  
  ## injury_zone
  "Meg3",
  "Gm10076",
  # "Gm39459",
  # "Gm47283",
  "Snhg1",
  
  "n-R5s136",
  "Rpph1",
  "Gm25939"
  
  
)

injury.line <- lapply(
  tmp.feat,
  FUN=function(FEAT){
    tmp.df = cbind(
      skm.merged@meta.data[skm.merged$polyA=="yPAP",],
      t(
        GetAssayData(
          skm.merged,
          assay="kallisto_collapsed"
        )[c(FEAT,"Gapdh"),skm.merged$polyA=="yPAP"]
      )
    )
    
    ggplot(
      data=tmp.df,
      mapping=aes_string(
        x="timepoint",
        y=as.name(FEAT),
        # group=interaction(timepoint,injury.zones),
        color="injury.zones"
      )
    )+
      stat_summary(
        data=tmp.df[tmp.df$injury.zones=="myofiber",],
        fun=mean,
        color=mckolors$primary[3],
        geom="line",
        size=1.5,alpha=0.7
      )+
      stat_summary(
        data=tmp.df[tmp.df$injury.zones=="myofiber",],
        fun=mean,
        color=mckolors$primary[3],
        alpha=1,
        size=2,
        geom="point"
      )+
      stat_summary(
        data=tmp.df[tmp.df$injury.zones=="myofiber",],
        fun.data=mean_sdl,
        color=mckolors$primary[3],
        geom="errorbar",
        alpha=0.7,width=0.5
      )+
      
      
      stat_summary(
        data=tmp.df[tmp.df$injury.zones=="injury_border",],
        fun=mean,
        color=mckolors$primary[2],
        size=1.5,alpha=0.7,
        geom="line"
      )+
      stat_summary(
        data=tmp.df[tmp.df$injury.zones=="injury_border",],
        fun=mean,
        color=mckolors$primary[2],
        size=2,alpha=1,
        geom="point"
      )+
      stat_summary(
        data=tmp.df[tmp.df$injury.zones=="injury_border",],
        fun.data=mean_sdl,
        color=mckolors$primary[2],
        geom="errorbar",
        alpha=0.7,width=0.5
      )+
      
      
      stat_summary(
        data=tmp.df[tmp.df$injury.zones=="injury_zone",],
        fun=mean,
        color=mckolors$primary[1],
        geom="line",alpha=0.7,
        size=1.5
      )+
      stat_summary(
        data=tmp.df[tmp.df$injury.zones=="injury_zone",],
        fun=mean,
        color=mckolors$primary[1],
        geom="point",alpha=1,
        size=2
      )+
      stat_summary(
        data=tmp.df[tmp.df$injury.zones=="injury_zone",],
        fun.data=mean_sdl,
        color=mckolors$primary[1],
        geom="errorbar",
        alpha=0.7,width=0.5
      )+
      scTheme$scatter+
      theme(
        panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        plot.margin = unit(rep(0,4),"cm"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_text(face="bold.italic",size = small.font)
      )+
      scale_x_continuous(
        breaks = c(0,2,5,7)
      )+
      guides(
        color = guide_legend(override.aes = list(size=2))
      )+
      labs(
        y=FEAT
      ) %>%
      return()
  }
)

injury.line[[length(injury.line)]] <- injury.line[[length(injury.line)]] +
  theme(
    axis.text.x = element_text(color="black", size=small.font),
    axis.title.x = element_text(size=small.font, color="black", face="bold")
  )+
  labs(x="Injury Timepoint (dpi)")

injury.line <- wrap_plots(
  injury.line,
  ncol=1,
  guides="collect"
)&theme(
  legend.position = "right"
)
injury.line

#

# SkM GE maps----
tmp.feat <- c(
  "Meg3",
  "Gm10076",
  "Rpph1"
)  

skm.map <- visListPlot(
  skm.list[meta_skm$polyA=="yPAP"],
  sample.titles = c("Uninjured", "2dpi", "5dpi","7dpi"),
  reduction="space",
  assay="kallisto_collapsed", 
  slot = 'data',
  pt.size=0.6,
  legend.position = "right",
  font.size = small.font,
  axis.title.angle.y=90,
  nrow = length(tmp.feat),
  flip_row_col = T,
  combine = T,
  verbose=F,
  colormap = "viridis",
  features=tmp.feat,
  alt.titles = stringr::str_remove(tmp.feat,pattern="mmu-")
)&theme(
  legend.text = element_text(size=small.font)
)&coord_fixed(
  ratio=1/1.6
)
skm.map
#
### scatter plot showing Notch1 vs. Gm13568 expression ----
notch.scatter <- lapply(
  skm.list[meta_skm$polyA=="yPAP"],
  FUN=function(VIS){
    DefaultAssay(VIS) <- "kallisto_collapsed"
    
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
        injury_zone=alpha(mckolors$primary[1],0.6),
        injury_border=alpha(mckolors$primary[2],0.6),
        myofiber=alpha(mckolors$primary[3],0.6)
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
#
# Unspliced vs. spliced maps ----
tmp.feat = c(
  "Myo3a",
  "Itga9",
  "Fcgr3",
  "Iqcg",
  "Mylip"
)
i <- 4:7

skm.uns.map <- visListPlot(
  skm.list[i],
  sample.titles = c("Uninjured","2dpi","5dpi","7dpi"),
  reduction="space",
  assay="unspliced",
  slot = 'data',
  pt.size=0.5,
  legend.position = "bottom",
  font.size = small.font,
  axis.title.angle.y=0,
  nrow= 1,
  flip_row_col = F,
  combine = T,
  verbose=F,
  colormap = "viridis",
  features=tmp.feat,
  alt.titles = stringr::str_remove(tmp.feat,pattern="mmu-")
)&theme(
  legend.text = element_text(size=small.font)
)&coord_fixed(
  ratio=1/1.6
)
skm.uns.map



skm.spl.map <- visListPlot(
  skm.list[i],
  sample.titles = c("Uninjured","2dpi","5dpi","7dpi"),
  reduction="space",
  assay="kallisto_collapsed",
  slot = 'data',
  pt.size=0.5,
  legend.position = "bottom",
  font.size = small.font,
  axis.title.angle.y=0,
  nrow = 1,
  flip_row_col = F,
  combine = T,
  verbose=F,
  colormap = "viridis",
  features=tmp.feat,
  alt.titles = stringr::str_remove(tmp.feat,pattern="mmu-")
)&theme(
  legend.text = element_text(size=small.font)
)&coord_fixed(
  ratio=1/1.6
)
skm.uns.map

#
# MiR maps ----
tmp.feat = c(
  "mmu-miR-1a-3p",
  "mmu-miR-206-3p"
  # "mmu-miR-133a-3p/133b-3p"
)
i <- 4:7

skm.mir.map <- visListPlot(
  skm.list[i],
  sample.titles = c("Uninjured","2dpi","5dpi","7dpi"),
  reduction="space",
  assay="mirge3",
  slot = 'data',
  pt.size=0.6,
  legend.position = "right",
  font.size = small.font,
  axis.title.angle.y=90,
  nrow = 2,
  flip_row_col = T,
  combine = T,
  verbose=F,
  colormap = "inferno",
  features=tmp.feat,
  alt.titles = stringr::str_remove(tmp.feat,pattern="mmu-")
)&theme(
  legend.text = element_text(size=small.font)
)&coord_fixed(
  ratio=1/1.6
)
skm.mir.map

#
# miR correlation ----
tmp.feat <- mir.rpm[,36:43]%>%rowSums()%>%sort(decreasing =T)%>%head(n=100)%>%names()
tmp.df <- data.frame(
  "smRNAseq"=mir.rpm[tmp.feat,36:43]%>%rowMeans(),
  "STRS"=mir.rpm[tmp.feat,50:53]%>%rowMeans(),
  
  "miR"=stringr::str_remove(tmp.feat,pattern = "mmu-"),
  "three"=mir.fa[tmp.feat,"three_prime"],
  "three_2"=mir.fa[tmp.feat,"three_prime_2"],
  "length"=mir.fa[tmp.feat,"length"],
  "A_count"=mir.fa[tmp.feat,"A_count"],
  "U_count"=mir.fa[tmp.feat,"U_count"],
  "G_count"=mir.fa[tmp.feat,"G_count"],
  "C_count"=mir.fa[tmp.feat,"C_count"],
  "GC"=(mir.fa[tmp.feat,"G_count"]+mir.fa[tmp.feat,"C_count"])/mir.fa[tmp.feat,"length"],
  "five"=mir.fa[tmp.feat,"five_prime"],
  "five"=mir.fa[tmp.feat,"five_prime_2"]
)

mir.scatter <- ggplot(
  tmp.df,
  aes(
    x=smRNAseq,
    y=STRS
  )
)+
  geom_smooth(
    # formula='y~x',
    method="lm",
    color="black"
  )+
  geom_point(
    color=mckolors$txg[4],
    size=1,
    alpha=0.9
  )+
  ggpmisc::stat_poly_eq(
    method = "lm",
    aes(
      label = paste(..eq.label.., ..rr.label.., sep = "~~~")
    ),
    parse = TRUE
  ) +
  ggrepel::geom_text_repel(
    aes(label=miR),
    color=mckolors$txg[4],
    face="italic"
  )+
  xlim(c(0,20))+
  ylim(c(0,20))+
  scTheme$scatter+
  theme(
    legend.position="bottom"
  )
mir.scatter

#
# GM13568 vs. Notch1 over time ----
ggplot(
  cbind(
    skm.merged@meta.data[skm.merged$polyA=="yPAP",],
    t(
      GetAssayData(
        skm.merged,
        assay="kallisto_collapsed"
      )[c("Notch1","Gm13568"),skm.merged$polyA=="yPAP"]
    )
  ),
  aes(
    x=timepoint,
    group=timepoint,
    color=polyA
  )
)+
  geom_violin(
    aes(
      y=log2(Gm13568/Notch1+1)
    )
  )
  
  stat_summary(
    aes(
      group=1,
      y=Gm13568
      ),
    fun=mean,
    color="red", 
    geom="line",
    group=1
    )+
  stat_summary(
    aes(
      group=1,
      y=Notch1
    ),
    fun=mean,
    color="blue", 
    geom="line",
    group=1
  )+
  scTheme$scatter+
  labs(
    y="Mean Log-normalized\nExpression"
  )

## wrap_plots ----
# wrap_plots(
#   wrap_plots(
#     plot_spacer(),
#     cluster.map%>% 
#       wrap_plots(
#         guides="collect",
#         nrow=1
#       )&coord_fixed(ratio=1/1.6)&theme(
#         legend.position="none"
#       ),
#     skm.map,
#     skm.mir.map,
#     wrap_plots(
#       mir.scatter&coord_fixed(1),
#       plot_spacer()
#     ),
#     # notch.scatter&coord_fixed(ratio=2/6),
#     nrow=5, heights=c(1,1,3,2,2.4),
#     guides="collect"
#   ),
#   injury.line&NoLegend(),
#   plot_spacer(),
#   widths=c(1,0.25,0.2)
# )
  
wrap_plots(
  wrap_plots( #1st column
    plot_spacer(),
    cluster.map%>% 
      wrap_plots(
        guides="collect",
        nrow=1
      )&coord_fixed(ratio=1/1.6)&theme(
        legend.position="none"
      ),
    skm.map,
    skm.mir.map,
    nrow=4,
    heights=c(1,1,3,2)
    # guides="collect"
  ),
  
  wrap_plots(#2nd column
    wrap_plots(
      injury.line&NoLegend(),
      plot_spacer(),
      nrow=1,
      widths=c(1,0.3)
    ),
    mir.scatter&xlim(c(5,20)), #&coord_fixed(1)
    ncol=1,
    heights=c(5,2)
  ),
  
  widths=c(1,0.5)
)
ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/Fig_SkM_v6.pdf",
  device="pdf",
  units="cm",
  width = 18*2,
  height = 16*2
)

#
# ----



















