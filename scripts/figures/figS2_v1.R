


suppressMessages(
  visListPlot(
    heart.list[c(2,4)],
    sample.titles = meta_heart$sample[c(2,4)],
    reduction = "space",
    assay="xGen.kallisto",
    slot="counts",
    pt.size = 0.1,
    features=reo.genes,
    # features=Features(heart.list[[4]],assay="kallisto")%>%tail(),
    # alt.titles = c("log2(Reovirus UMIs+1)","log2(xGen Reovirus UMIs+1)"),
    axis.title.angle.y = 0,
    legend.position = "right",
    combine=T,ncol = 4,
    colormap = "magma",
    colormap.direction = -1,
    colormap.same.scale = F
  )%>%
    wrap_plots(
      guides="collect"
    )&coord_fixed(
      ratio = 1.6
    )
)