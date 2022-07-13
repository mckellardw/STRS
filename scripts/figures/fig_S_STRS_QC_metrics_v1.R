# Rarefaction analysis of Visium and STRS libraries ----
# Load libs ----
library(Matrix)
library(dplyr)
library(Seurat)

library(ggplot2)
library(patchwork)
library(pals)
library(viridis)

source("/home/dwm269/DWM_utils/sc_utils/seurat_helpers/seutils.R")
source("/home/dwm269/DWM_utils/sc_utils/seurat_helpers/seuplots.R")
# Load data ----


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

# Load preprocessed data ----
load("/local/workdir/dwm269/totalRNA/spTotal/robjs/vis_list_v5.RData")

# SkM Plots ----
i = c(
  1:3, # SkM ctrl
  # 4,6,7,5,8,9, # SkM SUPERase
  10:13 # # SkM Protector
)

skm.qc.map <- visListPlot(
  vis.list[i],
  sample.titles = c(
    "SkM\n2dpi\nVisium",
    "SkM\n5dpi\nVisium",
    "SkM\n7dpi\nVisium",
    "SkM\nUninjured\nSTRS",
    "SkM\n2dpi\nSTRS",
    "SkM\n5dpi\nSTRS",
    "SkM\n7dpi\nSTRS"
  ),
  assay="STARsolo_collapsed",
  reduction="space",
  # slot = 'counts',
  pt.size=0.2,
  legend.position = "bottom",
  font.size = small.font,
  axis.title.angle.y=0,
  nrow = 1,
  combine = T,
  verbose=F,
  colormap="inferno",
  # colormap = mckolors$RdYlBl %>% rev(),
  features = c(
    "nCount_kallisto_collapsed",
    "nFeature_kallisto_collapsed",
    "pct_kallisto_spliced",
    "pct_kallisto_unspliced"
  )
)&coord_fixed(ratio=1/1.6)&
  scTheme$space&
  theme(
    legend.text = element_text(size=small.font)
  )

# Heart Plots ----
i = c(
  14:17
)

heart.qc.map <- visListPlot(
  vis.list[i],
  sample.titles = c(
    "Heart\nUninfected\nVisium",
    "Heart\nInfected\nVisium",
    "Heart\nUninfected\nSTRS",
    "Heart\nInfected\nSTRS"
  ),
  assay="STARsolo_collapsed",
  reduction="space",
  # slot = 'counts',
  pt.size=0.2,
  legend.position = "bottom",
  font.size = small.font,
  axis.title.angle.y=0,
  nrow = 1,
  combine = T,abs.heights = T,
  verbose=F,
  colormap="inferno",
  # colormap = mckolors$RdYlBl %>% rev(),
  features = c(
    "nCount_kallisto_collapsed",
    "nFeature_kallisto_collapsed",
    "pct_kallisto_spliced",
    "pct_kallisto_unspliced"
  )
)&coord_fixed(ratio=1.6)&
  scTheme$space&
  theme(
    legend.text = element_text(size=small.font)
  )

# QC scatter plots ----
i=c(
  1:3, # SkM ctrl
  # 4,6,7,5,8,9, # SkM SUPERase
  10:13, # # SkM Protector
  14:17
)
vis.merged <- merge(
  vis.list[[1]],
  vis.list[i],
  add.cell.ids = meta_vis$sample[c(1,i)]
)
vis.merged$sample <- factor(
  vis.merged$sample %>% 
    stringr::str_remove_all(pattern = "yPAP-")%>% 
    stringr::str_remove_all(pattern = "CTRL-"),
  levels = c(
    "yPAP-Pro_SkM-D0",
    "CTRL-SkM-D2",
    "yPAP-Pro_SkM-D2",
    "CTRL-SkM-D5",
    "yPAP-Pro_SkM-D5",
    "CTRL-SkM-D7",
    "yPAP-Pro_SkM-D7",
    "mock_D7PI",
    "yPAP-Pro_Heart-mock",
    "T1L_D7PI",
    "yPAP-Pro_Heart-D7T1L"
  )
)

# Count features & UMIs for each biotype
for(BT in unique(gtf.info$Biotype)){
  cat(BT,": ")
  tmp.mat <- GetAssayData(vis.merged,assay="kallisto_collapsed")
  tmp.feat <- gtf.info$GeneSymbol[gtf.info$Biotype==BT] %>%unique()
  cat(length(tmp.feat),"\n")
  
  vis.merged[[paste0("nCount_kc_",BT)]] <- colSums(tmp.mat[tmp.feat%in%rownames(tmp.mat),])
  vis.merged[[paste0("nFeature_kc_",BT)]] <- apply(
    tmp.mat[tmp.feat%in%rownames(tmp.mat),],
    2,
    FUN = function(SPOT) SPOT[SPOT>0]%>%length()
  )
  
  rm(tmp.mat,tmp.feat)
}


# median gene number of spot, median gene number of each gene type of spot, and total gene detected
bt.box <- lapply(
  c(
    "nCount_kallisto_collapsed",
    "nFeature_kallisto_collapsed",
    "nCount_kc_protein_coding",
    "nFeature_kc_protein_coding",
    "nCount_kc_rRNA",
    "nFeature_kc_rRNA",
    "nCount_kc_lncRNA",
    "nFeature_kc_lncRNA",
    "nCount_kc_miRNA",
    "nFeature_kc_miRNA",
    "nCount_kc_snoRNA",
    "nFeature_kc_snoRNA"
  ),
  FUN = function(FEAT) ggplot(
    vis.merged@meta.data,
    aes(
      x=sample,
      fill=polyA
    )
  )+
    geom_boxplot(
      color="black",
      aes_string(
        y=FEAT
      ),
      outlier.alpha = 0.5
    )+
    scTheme$scatter+
    theme(
      axis.title.y = element_text(size=small.font),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )+
    labs(
      fill="Method",
      y=stringr::str_remove(FEAT,"kc_")%>%
        stringr::str_remove("_collapsed")
    )+
    scale_y_continuous(labels=scales::scientific)+
    scale_fill_manual(values=mckolors$txg[c(1,4)])
)

bt.box[[length(bt.box)-1]] <- bt.box[[length(bt.box)-1]] + theme(
  axis.text.x = element_text(color = "black",angle=45,hjust=1,vjust=1)
)
bt.box[[length(bt.box)]] <- bt.box[[length(bt.box)]] + theme(
  axis.text.x = element_text(color = "black",angle=45,hjust=1,vjust=1)
)

bt.box <- bt.box %>%
  wrap_plots(
    guides="collect",
    ncol=2
  )&theme(legend.position = "bottom")
bt.box

# wrap & save  ----
wrap_plots(
  wrap_plots( # 10cm
    skm.qc.map,
    plot_spacer(),
    heart.qc.map,
    ncol=1,
    heights=c(2,0.05,2.3)
  ),
  wrap_plots(
    bt.box,
  plot_spacer(),
  ncol=1,
  heights = c(20,0.1)
  ),
  nrow=1,
  widths=c(10,9)
)
ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_STRSqc_v2.pdf",
  device="pdf",
  units="cm",
  width = 19*2,
  height = 20*2
)
