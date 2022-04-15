
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

# kallisto box plots -----
if(!exists("vis.merged")){
  vis.merged <- merge(
    vis.list[[1]],
    vis.list[2:length(vis.list)],
    add.cell.ids = meta_vis$sample
  )
}

samples.include <- c(
  "CTRL-SkM-D2","yPAP-Pro_SkM-D2",
  "T1L_D7PI" ,"yPAP-Pro_Heart-D7T1L"
)
# CELLS.KEEP = sample(Cells(vis.merged)[vis.merged$rnase_inhib!="SUPER"])
CELLS.KEEP = sample(Cells(vis.merged)[vis.merged$sample%in%samples.include])

tmp.plot <- lapply(
  list(
    "kal.protein_coding",
    "pct_kallisto_unspliced",
    
    "kal.rRNA",
    "kal.miRNA"
  ),
  FUN = function(Y){ 
    vis.merged$sample <- factor(
      vis.merged$sample,
      levels=samples.include
    )
    ggplot(
    vis.merged@meta.data[CELLS.KEEP,],
    aes(
      # x=tissue,
      x=sample,
      fill=polyA,
      color=polyA
    )
  )+
    geom_boxplot(
      aes_string(
        y=Y
      ),
      outlier.color = NA,
      outlier.size = 0,
      alpha=0.8
    )+
    geom_jitter(
      aes_string(
        y=Y
      ),
      color=gray(0.6),
      size=0.01,
      alpha=0.1
    )+
    scale_y_continuous(
      limits = c(0, NA),
      labels = scales::percent
    )+
    scTheme$scatter+
    theme(
      axis.text.x=element_text(angle=45,hjust=1),
      panel.grid.minor = element_blank()
    )+
    aes(stroke=pt.stroke)+
    scale_color_manual(
      values=mckolors$txg[c(1,4)]
    )+
    scale_fill_manual(
      values=mckolors$txg[c(1,4)]
    )%>%
      return()
  }
) 
tmp.plot[c(1:2)] <- lapply(
  tmp.plot[c(1:2)],
  FUN = function(PLOT) PLOT + 
    theme(
      axis.text.x = element_blank()
    )
)
tmp.plot[[1]] <- tmp.plot[[1]] + labs(y="mRNA")
tmp.plot[[2]] <- tmp.plot[[2]] + labs(y="Unspliced")
tmp.plot[[3]] <- tmp.plot[[3]] + labs(y="rRNA")
tmp.plot[[4]] <- tmp.plot[[4]] + labs(y="miRNA")

wrap_plots(
  tmp.plot,
  nrow=3,
  guides="collect"
)&theme(
  axis.title.x = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "bottom"
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/Fig1_box_v1.pdf",
  device="pdf",
  units="cm",
  width = 5*2,
  height = 8*2
)

# Reo/xGen maps ----
wrap_plots(
  visListPlot(
    list(heart.list[[2]]),
    sample.titles = meta_heart$sample[c(2)],
    reduction = "space",
    pt.size = 0.1,
    font.size=small.font,
    features=c("reo.log2p1"),
    alt.titles = c("log2(Reovirus UMIs+1)","log2(xGen Reovirus UMIs+1)"),
    axis.title.angle.y = 0,
    combine=F,
    colormap = "magma",
    colormap.direction = -1,
    colormap.same.scale = T
  )%>%
    wrap_plots(
    )&theme(
      legend.position="right",
      axis.title.y = element_blank()
    )&coord_fixed(
      ratio = 1.6
    ),
  
  visListPlot(
    list(heart.list[[4]]),
    sample.titles = meta_heart$sample[c(4)],
    reduction = "space",
    pt.size = 0.1,
    font.size=small.font,
    features=c("reo.log2p1","nCount_xGen.kallisto.log2p1"),
    alt.titles = c("log2(Reovirus UMIs+1)","log2(xGen Reovirus UMIs+1)"),
    axis.title.angle.y = 0,
    combine=F,
    colormap = "magma",
    colormap.direction = -1,
    colormap.same.scale = T
  )%>%
    wrap_plots(
    )&theme(
      legend.position="right",
      axis.title.y = element_blank()
    )&coord_fixed(
      ratio = 1.6
    ),
  nrow=1,
  guides="collect"
)



visListPlot(
  heart.list[c(2,4)],
  sample.titles = meta_heart$sample[c(2,4)],
  reduction = "space",
  pt.size = 0.1,
  font.size=small.font,
  features=c("reo.log2p1"),
  alt.titles = c("log2(Reovirus UMIs+1)"),
  axis.title.angle.y = 0,
  combine=F,
  colormap = "magma",
  colormap.direction = -1,
  colormap.same.scale = T
)%>%
  wrap_plots(
  )&theme(
    legend.position="right",
    axis.title.y = element_blank()
  )&coord_fixed(
    ratio = 1.6
  )
ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/Fig1_ReoMaps_v2.pdf",
  device="pdf",
  units="cm",
  width = 3*2,
  height = 5*2
)

# vln counts by tissue ----
out.plot <- list()
for(TISSUE in c("muscle", "heart")){
  out.plot[[TISSUE]] <- lapply(
    list(
      "nFeature_kallisto",
      "nCount_kallisto",
      "pct_kallisto_unspliced"
      # "nCount_TAR"
    ),
    FUN = function(Y) ggplot(
      vis.merged@meta.data[sample(Cells(vis.merged)[vis.merged$rnase_inhib!="SUPER" & vis.merged$tissue ==TISSUE]),],
      aes(
        x=sample,
        fill=polyA,
        color=polyA
      )
    )+
      geom_violin(
        aes_string(
          y=Y
        ),
        alpha=0.8
      )+
      geom_jitter(
        aes_string(
          y=Y
        ),
        color=gray(0.6),
        size=0.1,
        alpha=0.1
      )+
      scale_y_continuous(
        limits = c(0, NA)
      )+
      scTheme$scatter+
      theme(
        axis.text.x=element_text(angle=45,hjust=1),
        panel.grid.minor = element_blank()
      )+
      # scale_y_continuous(labels=scales::percent)+
      scale_color_manual(
        values=mckolors$txg[c(1,4)]
      )+
      scale_fill_manual(
        values=mckolors$txg[c(1,4)]
      )
  ) %>%
    wrap_plots(
      guides="collect",
      nrow=1
      
    )
}

wrap_plots(
  out.plot,
  ncol=1,
  guides="collect"
)&theme(
  axis.title.x = element_blank()
)




# example ncRNA dotplot ----
samples.include <- c(
  "CTRL-SkM-D2",
  "yPAP-Pro_SkM-D2",
  "T1L_D7PI",
  "yPAP-Pro_Heart-D7T1L"
)

CELLS.KEEP = sample(Cells(vis.merged)[vis.merged$sample %in% samples.include])

Idents(vis.merged)<-"sample"
#TODO- fix x-axis order
# add more genes?
#  better snoRNA
# examples from other biotypes
# remove ctrl samples?
VlnPlot(
  vis.merged,
  features=c(
    "Rny3",
    "7SK",
    "ENSMUSG00002075551",
    "Gm42826",
    "mt-Ta",
    "Bc1",
    "Snord13",
    "Mir6236"
  ),
  pt.size=0,
  group.by="sample",
  assay = "kallisto_collapsed",
  idents = samples.include,
  combine = F
)%>%
  lapply(
    FUN = function(X) X +
      scTheme$vln +
      theme(
        axis.title.x = element_blank(),
        legend.position="none",
        plot.title = element_text(face="bold.italic",hjust=0.5),
        axis.text.x = element_text(hjust=1,angle=45),
        axis.line=element_line(color="black")
      )
  ) %>%
  wrap_plots(
    nrow=1,
    guides="collect"
  )

grepGenes(heart.list[[4]],assay="kallisto_collapsed", pattern="")[1:1000]%>%
  lapply(
    FUN=function(X) c(X,head(gtf.info$Biotype[gtf.info$GeneSymbol==X], n=1))
  ) %>%
  lapply(
    FUN=function(X) if(X[2]=="protein_coding"){
      return(NULL)
    }else(
      return(X)
    )
  )%>%
  do.call(what=rbind)

## Genes Dot plot! ----
DotPlot(
  vis.merged,
  features=c(
    "Gapdh","Ckm",
    "Malat1","Neat1",
    "Bc1",
    # "Rn18s-rs5",
    "Gm42826",
    "Gm37357",
    "Rny1",
    "Rny3",
    "7SK",
    # "ENSMUSG00002075551",
    "Rpph1",
    "mt-Ta",
    "mt-Th",
    "Snord118",
    "Snord35b",
    "Mir6236"
  ),
  scale = F,
  group.by="sample",
  assay = "kallisto_collapsed",
  idents = samples.include
)+
  scale_color_viridis_c(option="viridis")+
  scTheme$dot

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/Fig1_dot_v1.pdf",
  device="pdf",
  units="cm",
  width = 10*2,
  height = 3*2
)

## Biotypes Dot plot! ----
DotPlot(
  vis.merged,
  features=c(
    # "kal.protein_coding",
    "kal.rRNA",
    "kal.Mt_rRNA",
    "kal.miRNA",
    "kal.lncRNA",
    "kal.Mt_tRNA",
    "kal.snoRNA",
    "kal.snRNA",   
    "kal.ribozyme",  
    "kal.misc_RNA",              
    "kal.scaRNA",
    "kal.pseudogene",
    "kal.polymorphic_pseudogene",
    "kal.TEC"
  ),
  scale = F,
  group.by="sample",
  assay = "kallisto_collapsed",
  idents = samples.include
)+
  scale_color_viridis_c(option="plasma")+
  scTheme$dot

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/Fig1_biotype_dot_v1.pdf",
  device="pdf",
  units="cm",
  width = 10*2,
  height = 6*2
)

# example ncRNA maps ----
tmp.feat <- c(
  "Gapdh",
  # "Malat1",
  "ENSMUSG00002075551",
  "Gm42826",
  "7SK",
  "mt-Th",
  # "Gm22980",
  # "Snord38a",
  # "Gm26330",
  # "Snord55"
  "Rny1",
  "Snord118"
)

wrap_plots(
  visListPlot(
    vis.list[c(1,11)],
    sample.titles = meta_vis$sample[c(1,11)],
    reduction = "space",
    assay="kallisto_collapsed",
    pt.size = 0.3,
    font.size=small.font,
    features=tmp.feat,
    # alt.titles = c("log2(Reovirus UMIs+1)","log2(xGen Reovirus UMIs+1)"),
    axis.title.angle.y = 0,
    combine=T,nrow=1,
    colormap = "viridis",
    colormap.direction = 1,
    colormap.same.scale = F
  )&theme(
      legend.position="bottom",
      axis.title.y = element_blank()
    )&coord_fixed(
      ratio = 1/1.6
    ),
  
  visListPlot(
    vis.list[c(15,17)],
    sample.titles = meta_vis$sample[c(15,17)],
    reduction = "space",
    assay="kallisto_collapsed",
    pt.size = 0.3,
    font.size=small.font,
    features=tmp.feat,
    # alt.titles = c("log2(Reovirus UMIs+1)",tmp.feat"log2(xGen Reovirus UMIs+1)"),
    axis.title.angle.y = 0,
    combine=T,nrow=1,
    colormap = "viridis",
    colormap.direction = 1,
    colormap.same.scale = F
  )&theme(
      legend.position="bottom",
      axis.title.y = element_blank()
    )&coord_fixed(
      ratio = 1.6
    ),
  nrow=2,
  heights=c(1,2.2)
  # guides="collect"
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/Fig1_ncMaps_v1.pdf",
  device="pdf",
  units="cm",
  width = 14*2,
  height = 9*2
)

# ----



















