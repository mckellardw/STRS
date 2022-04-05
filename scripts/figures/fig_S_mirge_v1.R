
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

tmp.feat = c(
  "mmu-let-7a-5p",
  "mmu-let-7b-5p",
  "mmu-let-7c-5p",
  "mmu-let-7d-5p",
  "mmu-let-7f-5p",
  "mmu-let-7g-5p",
  "mmu-let-7i-5p",
  "mmu-miR-1a-3p",
  "mmu-miR-21a-5p",
  "mmu-miR-22-3p",
  "mmu-miR-26a-5p",
  "mmu-miR-30a-5p",
  "mmu-miR-30e-5p",
  "mmu-miR-133a-3p/133b-3p",
  "mmu-miR-142a-3p" ,
  "mmu-miR-145a-5p",
  "mmu-miR-322-5p",
  "mmu-miR-3473b/3473e"
)

tmp.feat <- mir.rpm[,14:17]%>%rowSums()%>%sort(decreasing =T)%>%head(n=50)%>%names()

tmp.col.nums <- c(
  # 1:5,10:13,#male smRNA
  6:9, #female smRNA
  14:15, #mock smRNA
  16:17, #T1L smRNA
  20:21, #ctrl visium
  18:19 # Protector visium
  # 1:21 # all samples
)
tmp.df = data.frame(meta_smRNA[tmp.col.nums,c("source","chemistry","polyA")])
rownames(tmp.df)<-meta_smRNA$sample[tmp.col.nums]

heart.heat <- pheatmap::pheatmap(
  mir.rpm[tmp.feat, tmp.col.nums],
  annotation_col = tmp.df,
  border_color = "black",fontsize = small.font,
  labels_row = stringr::str_remove(tmp.feat,"mmu-"),
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
  colormap = "viridis",
  features=tmp.feat,
  alt.titles = stringr::str_remove(tmp.feat,"mmu-")
)&theme(
  legend.text = element_text(size=small.font)
)&coord_fixed(
  ratio=1.6
)
heart.mir.map


# top miR heatmap (SkM) ----

tmp.feat = c(
  "mmu-let-7a-5p",
  "mmu-let-7b-5p",
  "mmu-let-7c-5p",
  "mmu-let-7d-5p",
  "mmu-let-7f-5p",
  "mmu-let-7g-5p",
  "mmu-let-7i-5p",
  "mmu-miR-1a-3p",
  "mmu-miR-16-5p",
  "mmu-miR-21a-5p",
  "mmu-miR-22-3p",
  "mmu-miR-26a-5p",
  "mmu-miR-30a-5p",
  "mmu-miR-126a-5p",
  "mmu-miR-133a-3p/133b-3p",
  "mmu-miR-142a-3p" ,
  "mmu-miR-143-3p",
  "mmu-miR-206-3p",
  "mmu-miR-451a",
  "mmu-miR-486a-5p/486b-5p"
)

tmp.feat <- mir.rpm[,36:43]%>%rowSums()%>%sort(decreasing =T)%>%head(n=50)%>%names()

tmp.col.nums <- c(
  # 22:26,32:35,#male smRNA
  27:31, #female smRNA
  36:43, # Ntx smRNA 
  54:56, #ctrl visium
  # 44,46,47,45,48,49, #SUPERase visium
  50:53 # Protector visium
  # 18:44 # all samples
  )

tmp.df = data.frame(meta_smRNA[tmp.col.nums,c("source","chemistry","polyA")])
rownames(tmp.df)<-meta_smRNA$sample[tmp.col.nums]

skm.heat <- pheatmap::pheatmap(
  mir.rpm[tmp.feat, tmp.col.nums],
  annotation_col = tmp.df,
  labels_row = stringr::str_remove(tmp.feat,"mmu-"),
  border_color = "black",
  fontsize = small.font,
  cluster_rows = F,
  cluster_cols = F
) %>% ggplotify::as.ggplot()
  skm.heat

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
  colormap = "viridis",
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
  nrow=1,
  widths = c(1,1.2),
  guides='collect'
)

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_miR_v1.pdf",
  device="pdf",
  units="cm",
  width = 20*2,
  # width=11*2,
  height = 18*2
)


# miR correlation - SkM ----

tmp.feat <- mir.rpm[,36:43]%>%rowSums()%>%sort(decreasing =T)%>%head(n=200)%>%names()

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

#COnvert reference sequences from DNA to RNA (T->U)
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

mir.scatter <- list()
mir.scatter[1:2]<- lapply( #discrete variables...
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
  
mir.scatter[3:4]<-lapply( #continuous variables
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

wrap_plots(
  mir.scatter,
  ncol = 1
)


ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_miR_correlation_v1.pdf",
  device="pdf",
  units="cm",
  width = 4*2,
  height = 15*2
)
