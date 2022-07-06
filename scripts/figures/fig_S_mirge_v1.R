
require("pheatmap")
require("ggplot2")

# Load processed miRNA counts ----
#TODO

# biotype bar plots ----
tmp.vars = c(
  "Filtered.miRNA.Reads", "Hairpin.miRNAs",
  "mature.tRNA.Reads","primary.tRNA.Reads",
  "snoRNA.Reads"
  # "rRNA.Reads",
  # "ncRNA.others"
  # "mRNA.Reads",
  # "Remaining.Reads"     
)

tmp <- reshape2::melt(
  df.anno,
  id.vars=c("sample","chemistry"), 
  measure.vars=tmp.vars,
  variable.name="biotype"
)
tmp$sample <- factor(
  tmp$sample,
  levels=levels(meta_smRNA$sample)
)

skm.bar <- ggplot(
  tmp,
  aes(
    x=sample,
    y=value,
    fill=biotype
  )
)+
  geom_col(
    position="fill"
  )+
  scale_fill_manual(values=mckolors$polychrome[c(6:12,14:32)])+
  scale_y_continuous(labels = scales::percent)+
  labs(
    color="non-mRNA biotypes"
  )+
  scTheme$bar+
  theme(
    axis.text.x=element_text(angle=45, hjust=1, vjust=1),
    axis.text.y = element_text(angle=0),
    legend.position="right",
    legend.title = element_text(face="bold"),
    axis.title = element_blank()
  )

# top miR heatmap (heart) ----

tmp.feat <- mir.rpm[,14:17]%>%rowSums()%>%sort(decreasing =T)%>%head(n=50)%>%names()

tmp.col.nums <- c(
  1:9,#female smRNA
  # 10:13, #male smRNA
  14:15, #mock smRNA
  16:17, #T1L smRNA
  20:21, #ctrl visium
  18:19 # Protector visium
  # 1:21 # all samples
)
tmp.df = data.frame(meta_smRNA[tmp.col.nums,c(
  "source",
  "chemistry"
  # "polyA"
)])
rownames(tmp.df)<-meta_smRNA$sample[tmp.col.nums]

heart.heat <- pheatmap::pheatmap(
  t(mir.rpm[tmp.feat, tmp.col.nums]),
  annotation_row = tmp.df,
  border_color = "black",fontsize = small.font,
  labels_col= stringr::str_remove(tmp.feat,"mmu-"),
  # cellwidth = 10,
  color = viridis(42, option="inferno"),
  cluster_rows = F,
  cluster_cols = F
)%>% ggplotify::as.ggplot()
heart.heat
#
# miR maps with heart data ----
tmp.feat = c(
  "mmu-miR-1a-3p",
  "mmu-let-7c-5p"
  # "mmu-miR-206-3p",
  # "mmu-miR-133a-3p/133b-3p"
)

heart.mir.map <- visListPlot(
  heart.list,
  sample.titles = stringr::str_remove_all(meta_heart$sample,pattern = "Vis_") %>%
    stringr::str_remove_all(pattern ="_SkM"),
  reduction="space",
  assay="mirge3",
  slot="data",
  pt.size=0.4,
  legend.position = "bottom",
  font.size = small.font,
  axis.title.angle.y=0,
  nrow = 1,
  combine = T,
  verbose=F,
  colormap = "inferno",
  features=tmp.feat,
  alt.titles = stringr::str_remove(tmp.feat,"mmu-")
)&theme(
  legend.text = element_text(size=small.font)
)&coord_fixed(
  ratio=1.6
)
heart.mir.map


# top miR heatmap (SkM) ----

tmp.feat <- mir.rpm[,36:43]%>%rowSums()%>%sort(decreasing =T)%>%head(n=50)%>%names()

tmp.col.nums <- c(
  22:31,#female smRNA
  # 32:35, #male smRNA
  36:43, # Ntx smRNA 
  54:56, #ctrl visium
  # 44,46,47,45,48,49, #SUPERase visium
  50:53 # Protector visium
  # 18:44 # all samples
  )

tmp.df = data.frame(meta_smRNA[tmp.col.nums,c(
  "source",
  "chemistry"
  # "polyA"
  )])
rownames(tmp.df)<-meta_smRNA$sample[tmp.col.nums]

skm.heat <- pheatmap::pheatmap(
  t(mir.rpm[tmp.feat, tmp.col.nums]),
  annotation_row  = tmp.df,
  labels_col= stringr::str_remove(tmp.feat,"mmu-"),
  border_color = "black",
  color = viridis(42,option = "inferno"),
  fontsize = small.font,
  # cellwidth = 10,
  cluster_rows = F,
  cluster_cols = F
) %>% ggplotify::as.ggplot()
skm.heat
#
# miR maps with SkM data ----
tmp.feat = c(
  "mmu-miR-1a-3p",
  "mmu-miR-206-3p",
  "mmu-miR-133a-3p/133b-3p"
)

skm.mir.map <- visListPlot(
  skm.list,
  sample.titles = stringr::str_remove_all(meta_skm$sample,pattern = "Vis_") %>%
    stringr::str_remove_all(pattern ="_SkM"),
  reduction="space",
  assay="mirge3",
  slot = 'data',
  pt.size=0.4,
  legend.position = "bottom",
  font.size = small.font,
  axis.title.angle.y=0,
  nrow = 1,
  combine = T,
  verbose=F,
  colormap = "inferno",
  features=tmp.feat
)&theme(
  legend.text = element_text(size=small.font)
)&coord_fixed(
  ratio=1/1.6
)
skm.mir.map

# wrap plots and save ----

wrap_plots(
  heart.heat,
  skm.heat,
  ncol=1,
  heights = c(1,1.6),
  guides='collect'
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_miR_v5.pdf",
  device="pdf",
  units="cm",
  width = 20*2,
  # width=11*2,
  height = 19*2
)


# miR correlation - SkM ----

tmp.feat <- mir.rpm[,36:43]%>%rowSums()%>%sort(decreasing =T)%>%head(n=250)%>%names()

mir.rpm[tmp.feat,36:43]%>%rowSums()

mirbase <- read.csv("/workdir/dwm269/genomes/mirbase/mirna.txt",header=F,sep="\t")
mirbase <- mirbase[grepl(mirbase$V3,pattern = "mmu"),]

mir.fa <- read.csv("/workdir/dwm269/genomes/mm39_all/STAR_smRNA/smRNA_Mus_musculus.fa",header=F)
tmp <- list()
for(i in seq(1,nrow(mir.fa),2)){
  tmp[[i]] <- data.frame(
    "mir"=stringr::str_remove(mir.fa[i,1],pattern = ">"),
    "seq"=mir.fa[i+1,1]
  )
}
mir.fa <- do.call(rbind, tmp)
rownames(mir.fa)<- mir.fa$mir

#Convert reference sequences from DNA to RNA (T->U)
mir.fa$seq <- stringr::str_replace_all(mir.fa$seq, pattern="T",replacement = "U")

# Add metadata on miR sequences
mir.fa$three_prime <- lapply(
  mir.fa$seq,
  FUN=function(SEQ) substr(SEQ,start=nchar(SEQ),stop=nchar(SEQ))
) %>% unlist()
mir.fa$three_prime_2 <- lapply(
  mir.fa$seq,
  FUN=function(SEQ) substr(SEQ,start=nchar(SEQ)-1,stop=nchar(SEQ))
) %>% unlist()

mir.fa$five_prime <- lapply(
  mir.fa$seq,
  FUN=function(SEQ) substr(SEQ,start=1,stop=1)
) %>% unlist()
mir.fa$five_prime_2 <- lapply(
  mir.fa$seq,
  FUN=function(SEQ) substr(SEQ,start=1,stop=2)
) %>% unlist()

mir.fa$length <- lapply(
  mir.fa$seq,
  FUN=function(SEQ) nchar(SEQ)
) %>% unlist()

mir.fa$A_count <- lapply(
  mir.fa$seq,
  FUN=function(SEQ) stringr::str_count(SEQ,pattern = "A")
) %>% unlist()

mir.fa$U_count <- lapply(
  mir.fa$seq,
  FUN=function(SEQ) stringr::str_count(SEQ,pattern = "U")
) %>% unlist()

mir.fa$G_count <- lapply(
  mir.fa$seq,
  FUN=function(SEQ) stringr::str_count(SEQ,pattern = "G")
) %>% unlist()

mir.fa$C_count <- lapply(
  mir.fa$seq,
  FUN=function(SEQ) stringr::str_count(SEQ,pattern = "C")
) %>% unlist()

tmp.df <- data.frame(
  "smRNAseq"=mir.rpm[tmp.feat,36:43]%>%rowMeans(),
  "STRS"=mir.rpm[tmp.feat,50:53]%>%rowMeans(),
  "Isakova2020"=mir.rpm[tmp.feat,22:31]%>%rowMeans(),
  
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

## Plot with mir names
mir.scatter <- ggplot(
  tmp.df,
  aes(
    x=smRNAseq,
    y=Isakova2020
  )
)+
  geom_abline()+
  geom_point(
    aes_string(),
    color=mckolors$txg[4],
    alpha=0.7
  )+
  geom_smooth( # all of the top miRNAs
    # formula='y~x',
    method="lm",
    # color=mckolors$txg[1]
    color="black"
  )+
  ggpmisc::stat_poly_eq(
    method = "lm",
    aes(
      label = paste(..eq.label.., ..rr.label.., sep = "~~~")
    ),
    parse = TRUE
  ) +
  
  # geom_smooth( #only co-detected miRNAs
  #   data = tmp.df[tmp.df$STRS>0,],
  #   # formula='y~x',
  #   method="lm",
  #   # color=mckolors$txg[1]
  #   color="blue"
  # )+
  # ggpmisc::stat_poly_eq(
  #   data = tmp.df[tmp.df$STRS>0,],
  #   method = "lm",
  #   aes(
  #     label = paste(..eq.label.., ..rr.label.., sep = "~~~")
  #   ),
  #   color="blue",
  #   label.y = 30,
  #   parse = TRUE
  # ) +
  
  ggrepel::geom_text_repel(
    aes(label=miR),
    color=mckolors$txg[4]
  )+
  xlim(c(0,20))+
  ylim(c(0,20))+
  scale_color_manual(values=mckolors$colblind_8)+
  scTheme$scatter+
  theme(
    legend.position="bottom"
  )
mir.scatter
#
## plots with colors by sequence characteristics ----
mir.scatter.list <- list()
mir.scatter.list[1:2]<- lapply( #discrete variables...
    c(
      "three","five"
    ),
    FUN=function(COLOR)
      ggplot(
        tmp.df,
        aes(
          x=smRNAseq,
          y=STRS
        )
      )+
      geom_abline()+
      geom_point(
        aes_string(
          color=COLOR
        ),
        alpha=0.7
      )+
      # ggrepel::geom_text_repel(aes(label=miR))+
      xlim(c(0,20))+
      ylim(c(0,20))+
      scale_color_manual(values=mckolors$colblind_8)+
      scTheme$scatter+
      theme(
        legend.position="bottom"
      )
  )
  
mir.scatter.list[3:4]<-lapply( #continuous variables
  c("GC","length"),
  FUN=function(COLOR)
    ggplot(
      tmp.df,
      aes(
        x=smRNAseq,
        y=STRS
      )
    )+
    geom_abline()+
    geom_point(
      aes_string(
        color=COLOR
      ),
      alpha=0.7
    )+
    # ggrepel::geom_text_repel(aes(label=miR))+
    xlim(c(0,20))+
    ylim(c(0,20))+
    scale_color_viridis(option="plasma")+
    scTheme$scatter+
    theme(
      legend.position="bottom"
    )
)


# wrap & save ----
wrap_plots(
  mir.scatter,
  wrap_plots(
    mir.scatter.list,
    nrow=1
  )&theme(
    axis.title=element_blank()
  ),
  ncol = 1,
  heights=c(4,1)
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_miR_correlation_v2.pdf",
  device="pdf",
  units="cm",
  width = 16*2,
  height = 13*2
)
