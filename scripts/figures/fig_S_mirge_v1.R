
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
  # "mmu-miR-706",
  "mmu-miR-3473b/3473e"
)

tmp.df = data.frame(meta_smRNA[1:17,c("chemistry","polyA")])
rownames(tmp.df)<-meta_smRNA$sample[1:17]

heart.heat <- pheatmap::pheatmap(
  mir.rpm[tmp.feat, 1:17],
  annotation_col = tmp.df,
  border_color = "black",
  cluster_rows = F,
  cluster_cols = F
)%>% ggplotify::as.ggplot()

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
  "mmu-miR-1a-1-5p",
  "mmu-miR-1a-2-5p",
  "mmu-miR-10a-5p",
  "mmu-miR-10b-5p",
  "mmu-miR-16-5p",
  "mmu-miR-21a-5p",
  "mmu-miR-22-3p",
  "mmu-miR-26a-5p",
  "mmu-miR-30a-5p",
  "mmu-miR-30e-5p",
  "mmu-miR-125b-1-3p", "mmu-miR-125b-2-3p", "mmu-miR-125b-5p",
  "mmu-miR-126a-5p",
  "mmu-miR-133a-3p/133b-3p",
  "mmu-miR-142a-3p" ,
  "mmu-miR-143-3p",
  "mmu-miR-190a-3p",  "mmu-miR-190a-5p", 
  "mmu-miR-206-3p",
  "mmu-miR-451a",
  "mmu-miR-486a-5p/486b-5p"
  # "mmu-miR-706",
)

tmp.col.nums <- c(18:44)

tmp.df = data.frame(meta_smRNA[tmp.col.nums,c("chemistry","polyA")])
rownames(tmp.df)<-meta_smRNA$sample[tmp.col.nums]

skm.heat <- pheatmap::pheatmap(
  mir.rpm[tmp.feat, tmp.col.nums],
  annotation_col = tmp.df,
  border_color = "black",
  cluster_rows = F,
  cluster_cols = F
) %>% ggplotify::as.ggplot()

# wrap plots and save 