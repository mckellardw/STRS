

# ctrl vs. polyA - gene-by-gene comparison ----
samples.include <- c(
  "CTRL-SkM-D2",
  "yPAP-Pro_SkM-D2",
  "T1L_D7PI",
  "yPAP-Pro_Heart-D7T1L"
)

all.features <- lapply(
  vis.list,
  FUN=function(SEU) Features(SEU, assay = "kallisto_collapsed")
) %>%
  unlist() %>%
  unique()
nc.features <- gtf.info$GeneSymbol[gtf.info$GeneSymbol %in% all.features & gtf.info$Biotype!="protein_coding"]
pc.features <- gtf.info$GeneSymbol[gtf.info$GeneSymbol %in% all.features & gtf.info$Biotype=="protein_coding"]

tmp <-AverageExpression(
  vis.merged,
  group.by = "sample",
  slot="counts",
  assays="kallisto_collapsed"
)$kallisto_collapsed[,samples.include]

get_biotype <- function(GENE){ 
  GENE=stringr::str_split(GENE,pattern="\\.")[[1]][1]
  out=gtf.info$Biotype[gtf.info$GeneSymbol==GENE]%>%head(n=1)
  if(is.character(out)){
    return(out)
  }else{
    print(GENE)
    return("Unknown")
  }
}


df <- data.frame(
  gene=rep(rownames(tmp),2),
  ctrl_expression = c(tmp[,"CTRL-SkM-D2"], tmp[,"T1L_D7PI"]),
  ypap_expression = c(tmp[,"yPAP-Pro_SkM-D2"], tmp[,"yPAP-Pro_Heart-D7T1L"]),
  # biotype = rep(lapply(rownames(tmp), get_biotype)%>%unlist(),2),
  tissue = c(rep("SkM",nrow(tmp)),rep("Heart",nrow(tmp)))
)
df$pc <- df$gene %in% pc.features

#plot!
ggplot(
  df[sample(rownames(df)),],
  aes(
    x=log1p(ctrl_expression),
    y=log1p(ypap_expression),
    # color=pc,
    shape=tissue
  )
)+
  geom_abline()+
  geom_point(
    alpha=0.6,
    size=1.5,
    color=mckolors$txg[4]
  )+
  # geom_smooth(
  #   method="lm",
  #   formula='y ~ x'
  # )+
  ggrepel::geom_text_repel(
    data=df[abs(log2(df$ctrl_expression/df$ypap_expression))>2 & (log1p(df$ctrl_expression)>1|log1p(df$ypap_expression)>2) | df$ctrl_expression>2 | df$ypap_expression>1.2,],
    max.overlaps = 30,
    size=small.font/ggplot2::.pt,
    aes(label=gene)
  )+
  scTheme$scatter+
  theme(
    legend.position="none",
    panel.grid.minor = element_blank()
  )+
  labs(
    x="log10(Visium+1)",
    y="log10(STRS+1)"
  )+
  facet_grid(
    rows=vars(tissue),
    cols=vars(pc)
  )

ggsave(
  filename="/workdir/dwm269/totalRNA/spTotal/figures/FigS_genexgene_v1.pdf",
  device="pdf",
  units="cm",
  width = 18*2,
  height = 18*2
)

#