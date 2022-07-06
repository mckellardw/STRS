
# Gene biotype spatial maps ----
meta_vis <- read.csv("/workdir/dwm269/totalRNA/spTotal/resources/metadata_sheets/meta_sheet_visium.csv")
## Skeletal muscle ----
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
    assay="STARsolo_collapsed",
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
      # "kal.protein_coding",
      # "kal.rRNA",
      # "kal.Mt_rRNA",
      # "kal.miRNA",
      # "kal.lncRNA",
      # "kal.Mt_tRNA",
      # "kal.snoRNA",
      # "kal.snRNA",   
      # "kal.ribozyme",  
      # "kal.misc_RNA",              
      # "kal.scaRNA"
      # "Itgam",
      # "C1qa"
      
      "Yy1",
      "Myod1",
      "Myog",
      "Myh1",
      "Yaf2"
    )
  )&coord_fixed(ratio=1/1.6)&
    scTheme$space&
    theme(
      legend.text = element_text(size=6*2)
    )
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_biotypes_SkM.pdf",
  device="pdf",
  units="cm",
  width = 24*2,
  height = 10*2
)

### Biotype pie charts ----
#Load lists of genes for each biotype
# load("/workdir/dwm269/totalRNA/spTotal/resources/gene_lists/biotype_gene_lists.RData")
gene.list <- list()

for(BT in unique(gtf.info$Biotype)){
  gene.list[[BT]] <- gtf.info$GeneSymbol[gtf.info$Biotype == BT] %>% unique() 
}

gene.list <- gene.list[c(
  "protein_coding","rRNA","Mt_rRNA",
  "miRNA","lncRNA","Mt_tRNA",
  "snoRNA","snRNA","ribozyme","misc_RNA", 
  "scaRNA","TEC"
  # ,"sRNA","scRNA"          
)]

# "transcribed_processed_pseudogene"  "transcribed_unprocessed_pseudogene" "unprocessed_pseudogene"             
# "pseudogene"    "processed_pseudogene",   "unitary_pseudogene",                  
# [17] "polymorphic_pseudogene"             "transcribed_unitary_pseudogene"     "translated_unprocessed_pseudogene"  "TR_V_gene"                         
# [21] "TR_V_pseudogene"                    "TR_D_gene"                          "TR_J_gene"                          "TR_C_gene"                         
# [25] "TR_J_pseudogene"                    "IG_LV_gene"                         "IG_V_gene"                          "IG_V_pseudogene"                   
# [29] "IG_J_gene"                          "IG_C_gene"  
# [33] "IG_C_pseudogene"                    "IG_D_gene"                          "IG_D_pseudogene"                    "IG_pseudogene"   

i = c(
  1:3, # SkM ctrl
  # 4,6,7,5,8,9, # SkM SUPERase
  10:13 # # SkM Protector
)

bt.perc <- lapply(
  vis.list[i],
  FUN = function(SEU){
    expr.vec <- GetAssayData(
      SEU, 
      assay= "kallisto_collapsed",
      slot="counts"
    ) %>% 
      rowSums()
    
    total.counts <- sum(expr.vec)
    
    bt.list <- list()
    for(BT in names(gene.list)){
      tmp.vec <- expr.vec[gene.list[[BT]]]
      tmp.vec <- tmp.vec[!is.na(tmp.vec)]
      
      
      bt.list[[BT]] <- sum(tmp.vec) / total.counts * 100
    }
    
    out.df <- as.data.frame(bt.list)%>%
      reshape2::melt()
    
    colnames(out.df) <- c("Biotype", "percent")
    
    out.df$Sample <- rep(SEU$sample[1],nrow(out.df))
    out.df$Method <- rep(SEU$polyA[1],nrow(out.df))
    
    out.df <- out.df %>% mutate(
      csum = rev(cumsum(rev(percent))), 
      pos = percent/2 + lead(csum, 1),
      pos = if_else(is.na(pos), percent/2, pos)
    )
    
    return(out.df)
  }
) 
# %>% do.call(what = rbind)

tmp.colors <- mckolors$ldw29[1:length(gene.list)] #rainbow(n=length(names(gene.list)))
names(tmp.colors) <- names(gene.list)

# skm.donuts <-
  lapply(
    bt.perc,
    FUN = function(df) ggplot(
      df[df$percent>0,],
      aes(
        x=2,
        y=percent,
        fill=Biotype
      )
    )+
      geom_col(
        width=1,
        color="white"
      )+
      ggrepel::geom_label_repel(
        aes(
          y = pos, 
          label = paste0(round(percent,digits = 2), "%"),
          color=Biotype
        ),
        color="white",
        max.overlaps = 100
      ) +
      coord_polar("y", start=0)+
      xlim(.2,2.5)+
      scTheme$pie+
      theme(
        legend.position = "right"
      )+
      scale_fill_manual(values = tmp.colors)+
      scale_color_manual(values = tmp.colors)
  )%>%
wrap_plots(
  nrow=2,
  guides="collect"
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_biotypes_pie_skm_v1.pdf",
  device="pdf",
  units="cm",
  width = 18*2,
  height = 8*2
)

#
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
# &scale_color_viridis_c(
#     option = "plasma",
#     labels=scales::percent
#   )

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_biotypes_heart_v2.pdf",
  device="pdf",
  units="cm",
  width = 24*2,
  height = 11*2
)

### Biotype pie charts (heart) ----
i = c(
  14:17
)

bt.perc <- lapply(
  vis.list[i],
  FUN = function(SEU){
    expr.vec <- GetAssayData(
      SEU, 
      assay= "kallisto_collapsed",
      slot="counts"
    ) %>% 
      rowSums()
    
    total.counts <- sum(expr.vec)
    
    bt.list <- list()
    for(BT in names(gene.list)){
      tmp.vec <- expr.vec[gene.list[[BT]]]
      tmp.vec <- tmp.vec[!is.na(tmp.vec)]
      
      
      bt.list[[BT]] <- sum(tmp.vec) / total.counts * 100
    }
    
    out.df <- as.data.frame(bt.list)%>%
      reshape2::melt()
    
    colnames(out.df) <- c("Biotype", "percent")
    
    out.df$Sample <- rep(SEU$sample[1],nrow(out.df))
    out.df$Method <- rep(SEU$polyA[1],nrow(out.df))
    
    out.df <- out.df %>% mutate(
      csum = rev(cumsum(rev(percent))), 
      pos = percent/2 + lead(csum, 1),
      pos = if_else(is.na(pos), percent/2, pos)
    )
    
    return(out.df)
  }
) 
# %>% do.call(what = rbind)

tmp.colors <- mckolors$ldw29[1:length(gene.list)] #rainbow(n=length(names(gene.list)))
names(tmp.colors) <- names(gene.list)

# heart.donuts <-
lapply(
  bt.perc,
  FUN = function(df) ggplot(
    df[df$percent>0,],
    aes(
      x=2,
      y=percent,
      fill=Biotype
    )
  )+
    geom_col(
      width=1,
      color="white"
    )+
    ggrepel::geom_label_repel(
      aes(
        y = pos, 
        label = paste0(round(percent,digits = 2), "%")
        # color=Biotype
      ),
      color="white",
      max.overlaps = 100
    ) +
    coord_polar("y", start=0)+
    xlim(.2,2.5)+
    scTheme$pie+
    theme(
      legend.position = "right"
    )+
    scale_fill_manual(values = tmp.colors)+
    scale_color_manual(values = tmp.colors)
) %>%
  wrap_plots(
    nrow=2,
    guides="collect"
  )

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_biotypes_pie_heart_v1.pdf",
  device="pdf",
  units="cm",
  width = 18*2,
  height = 8*2
)
#
# xGen and STRS efficiency ----
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

# ----