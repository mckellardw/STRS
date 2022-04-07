# This script contains helper functinos that I commonly use in my Seurat analyses.
# Development versions of these and other Seurat helpers can be found here: https://github.com/mckellardw/DWM_utils/tree/main/sc_utils/seurat_helpers

########################################
## Helpers for grabbing/searching feature names
########################################

# Quickly return genes/feature names from a Seurat object
Features <- function(
  SEU,
  assay=NULL
){
  if(is.null(assay)){
    assay <- SEU@active.assay
  }

  return(rownames(GetAssayData(SEU,assay=assay)))
}

# Check gene names for a pattern, using grep
grepGenes <- function(
  SEU,
  pattern=NULL, # pattern to look for
  assay=NULL,
  filter.pattern=NULL, # pattern to remove
  sort.by = c("expression","abc"),
  verbose=T
){
  if(is.null(pattern)){
    if(verbose){message("Need a pattern to look for!")}
    return(NULL)
  }
  if(is.list(SEU)){
    if(verbose){message("Don't pass a list!")}
    return(NULL)
  }
  if(is.null(assay)){
    assay <- SEU@active.assay
  }

  genes = SEU@assays[[assay]]@counts@Dimnames[[1]]



  if(length(pattern)>1){
    if(verbose){
      message(paste0("Found ", length(genes), " features in the assay '",assay,"'..."))
      message(paste0("Looking for multiple patterns in these features..."))
    }
    out.genes <- lapply(
      pattern,
      FUN = function(PAT) genes[grep(pattern=PAT,x=genes)] #get genes
    )
    out.genes <- unlist(out.genes)
  }else{
    if(verbose){
      message(paste0("Found ", length(genes), " features in the assay '",assay,"'..."))
      message(paste0("Looking for '", pattern, "' in these features..."))
    }
    out.genes <- genes[grep(pattern=pattern,x=genes)] #get genes
  }

  if(length(out.genes)==0){
    message("Nothing found!\n")
    return(NULL)
  }

  if(!is.null(filter.pattern)){ # filter out filter.pattern
    if(verbose){message(paste0("Removing features containing '", filter.pattern,"' from output..."))}
    for(fp in filter.pattern){
      out.genes <- out.genes[!grepl(pattern=fp, x=out.genes)]
    }
  }

  if(is.null(sort.by[1])){
    if(verbose){message("Output genes not sorted...")}
  }else if(length(out.genes)==1){
    # Do nothing
  }else if(sort.by[1]=="expression"){
    if(verbose){message("Output genes sorted by expression...")}

    out.genes <- GetAssayData(
      SEU,
      assay=assay,
      slot="counts"
    )[out.genes,] %>%
        rowSums() %>% # sum counts across all cells
        sort(decreasing=T) %>% # sort genes according to total expression
        names() # get gene names back
  }else if(sort.by[1]=="abc"){
    if(verbose){message("Output genes sorted alphabetically...")}

    out.genes <- sort(out.genes, decreasing=T)
  }

  # Return matching gene names!
  return(
    out.genes
  )
}

# Convert ensembl IDs to gene IDs using a biomaRt reference
ens2gene <- function(
  ens=NULL, # vector of ensembl IDs to convert
  biomart.info=NULL, # biomaRt database
  ens.colname="ensembl_gene_id",
  gene.colname="mgi_symbol",
  ncores=1,
  force.unique=F, # switch to make gene names unique (adds ".1", ".2", etc.)
  verbose=F
){
  require(dplyr)

  if(is.null(ens)){
    message("Need ensembl IDs to convert!")
    return(NULL)
  }
  if(is.null(biomart.info)){
    message("Need biomaRt reference for conversion!")
    return(NULL)
  }
  if(!ens.colname %in% colnames(biomart.info)){
    message("ensembl ID column not found in biomaRt reference. Check input for 'ens.colname'")
    return(NULL)
  }
  if(!gene.colname %in% colnames(biomart.info)){
    message("Gene ID column not found in biomaRt reference. Check input for 'gene.colname'")
    return(NULL)
  }

  if(ncores==1){
    out <- lapply(
      ens,
      FUN = function(ENS){
        tmp = biomart.info[biomart.info[[ens.colname]]==ENS,]
        if(nrow(tmp==1)){ # only found one entry that matches
          feat = tmp[1,gene.colname]
        }else if(length(unique(tmp[,gene.colname])==1)){ # only found one unique gene name that matches
          feat = tmp[1,gene.colname]
        }else if(length(unique(tmp[,gene.colname])>1)){ # only found one unique gene name that matches
          if(verbose){message(paste0("Found multiple matches for ", ENS,", returning the first one I see..."))}
          feat = tmp[1,gene.colname]
        }else{
          if(verbose){cat(paste0("Nothing found for ", ENS,", returning ensembl ID"))}
          feat = ENS
        }

        if(feat==""){
          if(verbose){message(paste0("No gene name found in '", gene.colname,"', for",ENS, ", returning ensembl ID."))}
          feat = ENS
        }

        return(feat)
      }
    ) %>%
      unlist()
  }else if(ncores>1){
    require(parallel)
    if(verbose){message(paste0("Running on ", ncores," threads..."))}

    out <- mclapply(
      ens,
      FUN = function(ENS){
        tmp = biomart.info[biomart.info[[ens.colname]]==ENS,]
        if(nrow(tmp==1)){ # only found one entry that matches
          feat = tmp[1,gene.colname]
        }else if(length(unique(tmp[,gene.colname])==1)){ # only found one unique gene name that matches
          feat = tmp[1,gene.colname]
        }else if(length(unique(tmp[,gene.colname])>1)){ # only found one unique gene name that matches
          if(verbose){message(paste0("Found multiple matches for ", ENS,", returning the first one I see..."))}
          feat = tmp[1,gene.colname]
        }else{
          if(verbose){cat(paste0("Nothing found for ", ENS,", returning ensembl ID"))}
          feat = ENS
        }

        if(feat==""){
          if(verbose){message(paste0("No gene name found in '", gene.colname,"', for",ENS, ", returning ensembl ID."))}
          feat = ENS
        }

        return(feat)
      },
      mc.cores = ncores
    ) %>%
      unlist()
  }

  if(force.unique){
    out <- make.unique(out)
  }

  return(
    out
  )
}


# Collapse cell/nuclei/spot counts for multimapped genes - these are features with a period ("\\.") in their name
#TODO: parallelize!
collapseMultimappers <- function(
    SEU,
    assay=NULL,
    new.assay.name=NULL,
    verbose=F
){

  if(is.null(new.assay.name)){
    new.assay.name = paste0(assay,"_collpased")
    message("Using ",new.assay.name, " as new.assay.name...")
  }
  if(is.null(new.assay.name)){
    message("Need new.assay.name!")
    return(SEU)
  }

  SEU@active.assay <- assay

  multi.feats <- grepGenes( #Find genes with a period in them
    SEU,
    assay = assay,
    pattern="\\.",
    sort.by="abc",
    verbose=verbose
  )
  if(length(multi.feats)==0){
    message("No multimappers found!")
    return(SEU)
  }

  multi.patterns <- stringr::str_split( #extract actual gene names
    multi.feats,
    pattern = "\\.",
    n = 2
  ) %>%
    lapply(FUN=function(X) X[1]) %>%
    unlist() %>%
    unique()

  if(verbose){
    message(paste0("Found ", length(multi.patterns), " multimappers and ", length(multi.feats)," loci..."))
  }

  # Collapse counts for each gene
  mat.multi <- GetAssayData(
    SEU,
    assay=assay,
    slot="counts"
  )

  collapsed.list <- lapply(
    multi.patterns,
    FUN=function(X){
      tmp.genes = rownames(mat.multi)[grep(rownames(mat.multi),pattern=X)]
      tmp.mat = mat.multi[tmp.genes,]

      if(length(tmp.genes)==1){
        return(tmp.mat)
      }else{
        return(colSums(tmp.mat))
      }
    }
  )
  collapsed.mat <- do.call(rbind, collapsed.list) %>% as.sparse()
  rownames(collapsed.mat) <- multi.patterns

  # Add new assay with collapsed counts + the rest of the genes
  if(verbose){cat(paste0("Adding back ", nrow(collapsed.mat), " features...\n"))}

  solo.feats <- rownames(SEU)[!rownames(SEU)%in%c(multi.feats,multi.patterns)]

  out.mat <- rbind(
    GetAssayData(SEU,assay=assay, slot="counts")[solo.feats,],
    collapsed.mat
  )
  SEU[[new.assay.name]] <- CreateAssayObject(counts=out.mat)

  SEU@active.assay <- new.assay.name

  # Return Seurat object!
  return(SEU)
}


getSpatialLocation <- function(
  SEU,
  whitelist="/home/dwm269/DWM_utils/align_pipes/10x_kallisto/resources/barcodes_10x/visium-v1_coordinates.txt",
  assay="RNA",
  reduction.name = "space",
  verbose=F
){
  require(Seurat)
  require(dplyr)

  # Check assays
  if(!assay %in% Assays(SEU)){
    assay=SEU@active.assay
    if(verbose){message("Using the default assay...")}
  }

  # Read in whitelist coordinates
  bc.coords <- read.table(
    whitelist,
    row.names = 1,
    col.names = c(
      "X",
      "Y"
    )
  )

  # Build reduction based on spot barcodes & whitelist
  bcs <- Cells(SEU)
  if(verbose){
    message(paste0(
      table(bcs %in% rownames(bc.coords))["TRUE"], " out of ", length(bcs), "barcodes found in whitelist..."
    ))
  }

  #TODO- check/remove '-1'

  tmp.mat <- lapply(
    bcs,
    FUN=function(BC) bc.coords[BC,]
  ) %>%
    do.call(what=rbind)

  colnames(tmp.mat) <- paste0(reduction.name, 1:2)
  rownames(tmp.mat) <- bcs

  # Add reduc to seurat obj
  SEU[[reduction.name]] <- CreateDimReducObject(
    embeddings=as.matrix(tmp.mat),
    assay = assay,
    key = paste0(reduction.name,"_")
  )

  return(SEU)
}

########################################
## General Seurat workflow helpers
########################################
# Calculate the number of PCs that contain some proportion (default is 95%) of the variance
npcs <- function(
  SEU,
  var.total=0.95,
  reduction="pca"
){
  if(is.null(SEU@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }

  tmp.var <- (SEU@reductions[[reduction]]@stdev)^2
  var.cut <- var.total*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }

  return(n.pcs)
}



# Preprocessing wrapper function
#   (1) NormalizeData(SEU) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
#   (2) FindNeighbors %>% RunUMAP, FindClusters
seuPreProcess <- function(
  SEU=NULL,
  assay='RNA',
  n.pcs=50,
  res=0.8,
  verbose=F
){
  if(is.null(SEU)){
    cat("Need a Seurat object to preprocess!\n")
    return(NULL)
  }
  if(!assay %in% Assays(SEU)){
    cat(paste0(assay, " not found in the seurat object! Not preprocessed.\n"))
    return(SEU)
  }

  # NormalizeData(SEU) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  pca.name = paste0('pca_', assay)
  pca.key = paste0(pca.name,'_')
  umap.name = paste0('umap_', assay)

  SEU = NormalizeData(
    SEU
  ) %>% FindVariableFeatures(
    assay = assay,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = verbose
  ) %>% ScaleData(
    assay = assay
  ) %>% RunPCA(
    assay = assay,
    reduction.name = pca.name,
    reduction.key = pca.key,
    verbose = verbose,
    npcs = n.pcs
  )

  #find pcs to use
  n.pcs.use = npcs(
    SEU=SEU,
    var.total = 0.95,
    reduction = pca.name
  )

  # FindNeighbors %>% RunUMAP, FindClusters
  SEU <- FindNeighbors(
    SEU,
    reduction = pca.name,
    dims = 1:n.pcs.use,
    force.recalc = TRUE,
    verbose = verbose
  ) %>% RunUMAP(
    reduction = pca.name,
    dims = 1:n.pcs.use,
    verbose = verbose,
    reduction.name=umap.name
  )

  SEU@reductions[[umap.name]]@misc$n.pcs.used <- n.pcs.use

  SEU <- FindClusters(
    object = SEU,
    resolution = res,
    verbose = verbose
  )
  # SEU[[paste0('RNA_res.',res)]] <- as.numeric(SEU@active.ident)
  gc()
  return(
    tryCatch(
      SEU,
      error=function(e) NULL
    )
  )
}

# Add a new ident seurat metadata filed) based on a list of cell types
#
#     object:     seurat object
#     old.idents: name of the idents metadata you will be assigning cell types to
#     new.idents: vector of cell types, in order of cluster number
#     newName:    string of the new idents name
#
AddCellTypeIdents <- function(
  SEU=NULL, old.name, new.name=NULL, new.idents, verbose=FALSE
){
  old.idents = as.vector(names(table(SEU[[old.name]])))

  if(is.null(new.name)){
    cat("**Need a new.name for the new idents**\n")
  }else{
    SEU[[new.name]] <- as.vector(SEU[[old.name]])
    for(i in 1:length(old.idents)){
      if(verbose){cat("Adding ", new.idents[i],"...", sep = "")}
      SEU[[new.name]][ SEU[[new.name]]==old.idents[i] ] <- new.idents[i]
      if(verbose){cat("Done!\n", sep = "")}
    }
  }

   return(SEU)

}


# Add gene biotype percentages to a seurat object, given a biomaRt object.
seu_biotypes <- function(
  SEU,
  biomart=NULL, # biomaRt or data.frame containing gene biotypes
  gene.colname,
  biotype.colname,
  add.as=c("metadata","assay"), # how percent features should be added
  assay="RNA",
  prefix="pct.",
  scale=c(1,100),
  verbose=TRUE
){

  if(is.null(biomart)){
    message("Need a list of gene biotypes! Nothing done.")
    return(SEU)
  }

  if(add.as[1]=="assay"){
    message("add.as='assay' is not yet implemented")
    new.assay.name = paste0(assay,".biotypes")
    return(SEU)
  }

  if(verbose){message(paste0("Adding gene biotype percentage values as ", add.as, " ..."))}

  biotypes = unique(biomart[[biotype.colname]])
  for(biotype in biotypes){
    tmp.mart =  biomart[biomart[[biotype.colname]]==biotype,] #subset out current gene biotype
    tmp.feat = unique( #get unique gene names which are present in the SEU object
      tmp.mart[[gene.colname]][tmp.mart[[gene.colname]]%in%Features(SEU, assay=assay)]
    )

    if(length(tmp.feat)==0){
      if(verbose){message(paste0("  No ", biotype, " genes found..."))}
    }else{
      if(add.as[1]=="metadata"){
        SEU <- PercentageFeatureSet(
          SEU,
          col.name = paste0(prefix, biotype),
          assay = assay,
          features=tmp.feat
        )
        if(scale[1]==1){
          SEU[[paste0(prefix, biotype)]] <- SEU[[paste0(prefix, biotype)]]/100
        }else if(scale[1]==100){
          #Do nothing...
        }else{
          message("Given scale was not found. Scaling to 100...")
        }
      }
      if(add.as[1]=="assay"){
        #TODO
        message("Not implemented yet...")
      }
    }
  }
  return(SEU)
}

# Plotting Utils ----

# Generate feature plots given a list of Visium/Seurat objects
visListPlot <- function(
  seu.list,
  features=NULL,
  alt.titles=NULL, # alternative titles for genes/features being passed
  sample.titles=NULL, # Sample IDs (y-axis labels)
  assay='Spatial',
  reduction="space",
  slot="data",
  legend.position="bottom",
  pt.size=1,
  font.size=8,
  axis.title.angle.y=90, #y-axis title angle (parameter for ggplot::element_text)
  combine=TRUE,
  abs.heights=TRUE, # Use absolute heights to size each subplot
  nrow=NULL,
  ncol=NULL,
  flip_row_col=F, #Default is samples in rows, features in cols. Set to `T` for samples in cols, features in rows
  colormap="viridis", # either a viridis option, a vector of colors, or a list of options corresponding to `features`
  colormap.direction=1,
  colormap.same.scale=F, #whether (T) or not (F) to set all features to the same colormap scale
  na.value=gray(0.69), # color for na.value (spot where gene is not detected)
  verbose=FALSE
){
  require(dplyr)
  require(Seurat)
  require(ggplot2)
  require(viridis)

  if(is.null(seu.list)){message("List of Seurat objects to plot is NULL!")}

  if(verbose){cat(paste0("Plotting Visium data, using the assay ",assay,"!\n"))}

  # Check colormap
  if(length(colormap)==1){ # single colormap option
    if(!colormap %in% c("magma","inferno","plasma","viridis","cividis","A","B","C","D")){
      message(paste0("Color palette '", colormap, "' not found!"))
      return(NULL)
    }else{
      colormap = rep(
        x=colormap,
        length=length(features)
      )%>%
        as.list()
    }
  }else if(length(colormap)>1 & !is.list(colormap)){
    colormap = rep(
      x=list(colormap),
      length=length(features)
    )
  }else if(is.list(colormap)){ # list for feature-specific colormaps
    if(length(colormap)!= length(features)){
      message(paste0("Wrong number of color palettes!"))
      return(NULL)
    }
  }

  # Alternate feature names
  if(is.null(alt.titles)){
    alt.titles=features
  }

  # Check slot, assay length(s)
  if(length(assay)==1){
    assay = rep(assay,length(features))
  }else if(length(assay) != length(features)){
    message("`assay` and `features` lengths don't match!")
  }

  if(length(slot)==1){
    slot = rep(slot,length(features))
  }else if(length(slot) != length(features)){
    message("`slot` and `features` lengths don't match!")
  }

  # Check for genes
  tmp.features = paste0("tmp.",features) # place-holder name for features; allows assay-specific feature plotting in FeaturePlot

  seu.list <- lapply(
    seu.list,
    FUN = function(SEU){
      for(i in 1:length(features)){
        if(features[i] %in% colnames(SEU@meta.data)){
          SEU@meta.data[,tmp.features[i]] <- SEU@meta.data[,features[i]]

        }else if(!features[i] %in% Features(SEU, assay = assay[i])){
          # Stick a bunch of zeroes into the metadata to act like an undetected gene
          SEU@meta.data[,tmp.features[i]] <- rep(0,nrow(SEU@meta.data))

        }else{
          SEU@meta.data[,tmp.features[i]] <- GetAssayData(
            SEU,
            assay=assay[i],
            slot=slot[i]
          )[features[i],]
        }
      }
      return(SEU)
    }
  )

  # Get expression limits for each gene, across all datasets
  gene.lims <- mapply(
    FUN = function(FEAT, ASS, SLOT){
      out.max <- lapply(
        seu.list,
        FUN = function(SEU){
          if(FEAT %in% Features(SEU, assay=ASS)){
            return(max(GetAssayData(SEU,assay=ASS,slot=SLOT)[FEAT,]))
          }else if(FEAT %in% colnames(SEU@meta.data)){
            return(max(SEU@meta.data[,FEAT]))
          }else{
            if(verbose){message(FEAT, " not found!")}
            return(0)
          }
        }
      ) %>%
        unlist() %>%
        max()

      return(c(10^-100, out.max))
    },
    SIMPLIFY = F,
    FEAT = features,
    ASS = assay,
    SLOT = slot
  )

  # Set colormap scale to the global min/max of all features
  if(colormap.same.scale){
    gene.lims <- lapply(
      gene.lims,
      FUN=function(X) c( 10^-100, max(unlist(gene.lims)) )
    )
  }

  #DEPRECATED- just use coord_fixed()!
  # Get plot heights
  #TODO- build in coord_fixed()
  if(abs.heights){
    heights <- lapply(
      seu.list,
      FUN=function(SEU) abs(diff(range(SEU@reductions[[reduction]]@cell.embeddings[,2])))
    ) %>% unlist()
    if(verbose){
      message(paste0("Using these plot heights:"))
      print(heights)
    }
  }else{
    heights <- rep(1,length(features))
  }

  # Plot
  plot.list <- list()
  for(i in 1:length(features)){
    tmp <- lapply(
      seu.list,
      FUN = function(SEU){
        tmp.plot = FeaturePlot(
          SEU,
          slot = slot[i],
          features = tmp.features[i],
          pt.size = pt.size,
          reduction = reduction
        ) +
          theme(
            plot.margin = unit(rep(0,4), "inches"),
            axis.ticks = element_blank(),
            axis.text=element_blank(),
            axis.title = element_blank(),
            axis.line=element_blank(),
            plot.title = element_blank(),
            legend.position=legend.position,
            legend.title = element_text(size=font.size,face="bold", hjust=0.5),
            legend.text = element_text(size=font.size,face="bold")
          )

        # set colormap
        if(length(colormap[[i]])==1){
          if(verbose){message(paste0("Using viridis palette ", colormap))}
          tmp.plot = tmp.plot +
            scale_color_viridis(
              limits = unlist(gene.lims[i]),
              option = colormap[[i]],
              direction = colormap.direction,
              na.value = na.value
            )
        }else{
          if(verbose){message(paste0("Using custom color gradient"))}

          # Flip colormap if direction is `-1`
          if(colormap.direction==-1){
            colormap=rev(colormap[[i]])
          }

          tmp.plot = tmp.plot +
            scale_color_gradientn(
              limits=unlist(gene.lims[i]),
              colors=colormap[[i]],
              na.value = na.value
            )
        }
        return(tmp.plot)
      }
    )

    if(!flip_row_col){
      tmp[[1]] <- tmp[[1]] +
        theme(
          plot.title = element_text(
            size=font.size,
            face="bold.italic",
            vjust=1
          )
        ) +
        labs(title=alt.titles[i])
    }

    plot.list[[i]] <- tmp
  }

  if(flip_row_col){ #samples=columns; features=rows

    # Add sample titles
    for(i in 1: length(plot.list[[1]])){
      plot.list[[1]][[i]] <- plot.list[[1]][[i]] +
        theme(
          plot.title = element_text(
            size=font.size,
            face="bold",
            vjust=1
          )
        ) +
        labs(title=sample.titles[i])
    }

    # Feature titles on y-axis
    for(j in 1:length(plot.list)){
      plot.list[[j]][[1]] <- plot.list[[j]][[1]] +
        theme(
          axis.title.y = element_text(
            size=font.size,
            face="bold.italic",
            color="black",
            hjust=0.5,
            vjust=0.5,
            angle=axis.title.angle.y
          )
        )

      if(!is.null(alt.titles)){ # add feature titles
        plot.list[[j]][[1]] <- plot.list[[j]][[1]] +
          labs(y=alt.titles[j])
      }
    }

    #Wrap
    plot.list <- lapply(
      plot.list,
      FUN = function(X){
        wrap_plots(
          X,
          nrow=1,
          heights=heights,
          guides="collect"
        )&theme(
          legend.position=legend.position,
          legend.margin = margin(0,0,0,0,"inches")
        )
      }
    )

  }else{ #samples=rows; features=columns
    # Feature axis title
    for(i in 1:length(plot.list[[1]]) ){
      plot.list[[1]][[i]] <- plot.list[[1]][[i]] +
        theme(
          axis.title.y = element_text(
            size=font.size,
            face="bold",
            color="black",
            hjust=0.5,
            vjust=0.5,
            angle=axis.title.angle.y
          )
        )

      if(!is.null(sample.titles)){ # add sample titles
        plot.list[[1]][[i]] <- plot.list[[1]][[i]] +
          labs(y=sample.titles[i])
      }
    }

    # Sample axis title
    plot.list <- lapply(
      plot.list,
      FUN = function(X){
        wrap_plots(
          X,
          ncol=1,
          heights=heights,
          guides="collect"
        )&theme(
          legend.position=legend.position,
          legend.margin = margin(0,0,0,0,"inches")
        )
      }
    )
  }
  if(verbose){cat("Done plotting Visium data!\n")}

  if(combine){
    return(
      wrap_plots(
        plot.list,
        nrow=nrow,
        ncol=ncol,
        guides=NULL
      )
    )
  }else{
    return(plot.list)
  }
}


# Plot co-expression of any two genes/features
## Built on top of visListPlot()
visCoMap <- function(
  seu.list,
  features=NULL, # pair of genes
  alt.titles=NULL, # alternative titles for genes/features being passed
  sample.titles=NULL, # Sample IDs (y-axis labels)
  assay='Spatial', # Either a single assay, or a pair of assays (useful for spliced vs. unspliced plotting)
  reduction="space",
  slot="data", # Either a single slot, or a pair of slots
  legend.position="bottom",
  pt.size=1,
  font.size=8,
  axis.title.angle.y=90, #y-axis title angle (parameter for ggplot::element_text)
  combine=TRUE,
  abs.heights=TRUE, # Use absolute heights to size each subplot
  nrow=NULL,
  ncol=NULL,
  colormap=NULL, # either a viridis option, a vector of colors, or a list of options corresponding to `features`
  colormap.direction=1,
  colormap.same.scale=F, #whether (T) or not (F) to set all features to the same colormap scale
  na.value=gray(0.69), # color for na.value (spot where gene is not detected)
  comap.fxn = prod,
  coex.name = NULL, # plot title for computed co-expression values
  include.scatter = F,
  scatter.group.by= 'orig.ident',
  scatter.theme= NULL, # default is theme_minimal()
  verbose=FALSE
){
  require(dplyr)
  require(Seurat)
  require(ggplot2)
  require(viridis)

  if(is.null(seu.list)){message("List of Seurat objects to plot is NULL!")}

  if(verbose){cat(paste0("Plotting Visium data, using the assay ",assay,"!\n"))}

  # Check colormap
  if(is.null(colormap)){
    colormap=list(
      mckolors$blue_ramp,
      mckolors$red_ramp,
      mckolors$purple_ramp
    )
  }

  # Alternate feature names
  if(is.null(alt.titles)){
    alt.titles=features
  }

  # Check `features` length
  if(length(features)!=2){
    message(paste0("Only supports plotting co-=expression of 2 features... for now..."))
  }

  # Check slot, assay length(s)
  if(length(assay)==1){
    assay = rep(assay,2)
  }
  if(length(slot)==1){
    slot = rep(slot,2)
  }

  # Compute co-expression values and add back to Seurat objects (as entry in meta.data)
  if(is.null(coex.name)){
    coex.name = "coexpression.tmp.values"
  }

  maps.out <- lapply(
    seu.list,
    FUN = function(SEU){
      tmp.coex <- list()

      DefaultAssay(SEU) <- assay[1]
      if(features[1] %in% Features(SEU,assay=assay[1])){
        tmp.coex[[1]] <- FetchData(
          object = SEU,
          vars = features[1],
          slot = slot[1]
        )
      }else{
        if(verbose){message(paste0(features[1], " is missing from `", assay[1], "`"))}
        tmp.coex[[1]] <- rep(0, ncol(SEU))

      }

      DefaultAssay(SEU) <- assay[2]
      if(features[2] %in% Features(SEU,assay=assay[2])){
        tmp.coex[[2]] <- FetchData(
          object = SEU,
          vars = features[2],
          slot = slot[2]
        )
      }else{
        if(verbose){message(paste0(features[2], " is missing from `", assay[2], "`"))}
        tmp.coex[[2]] <- rep(0, ncol(SEU))
      }

      tmp.coex <- do.call(
        cbind,
        tmp.coex
      ) %>%
        apply(
          # X = tmp.coex,
          MARGIN = 1,
          FUN = prod # Function for co-expression
        )

      # SEU <- AddMetaData(
      #   object = SEU,
      #   metadata = tmp.coex,
      #   col.name = coex.name
      # )
      SEU@meta.data[,coex.name] <- tmp.coex
      return(SEU)
    }
  ) %>%
    visListPlot( # Plot all three with visListPlot
      # seu.list=seu.list,
      features=c(features,coex.name), # pass two features plus the co-expression values
      alt.titles=c(alt.titles, coex.name),
      sample.titles=sample.titles,
      reduction=reduction,
      assay=c(assay,assay[1]), #extra assay entry for the coex data
      slot=c(slot,slot[1]),#extra assay entry for the coex data
      legend.position=legend.position,
      pt.size=pt.size,
      font.size=font.size,
      axis.title.angle.y=axis.title.angle.y,
      combine=combine,
      abs.heights=abs.heights,
      nrow=nrow,
      ncol=ncol,
      colormap=colormap,
      colormap.direction=colormap.direction,
      colormap.same.scale=colormap.same.scale,
      na.value=na.value,
      verbose=verbose
    )

  #Scatter plot
  if(include.scatter){
    if(is.null(scatter.theme)){
      scatter.theme <- theme_minimal()
    }

    # Get expression limits for each gene, across all datasets
    gene.lims <- mapply(
      FUN = function(FEAT, ASS, SLOT){
        out.max <- lapply(
          seu.list,
          FUN = function(SEU){
            if(FEAT %in% Features(SEU, assay=ASS)){
              return(max(GetAssayData(SEU,assay=ASS,slot=SLOT)[FEAT,]))
            }else if(FEAT %in% colnames(SEU@meta.data)){
              return(max(SEU@meta.data[,FEAT]))
            }else{
              if(verbose){message(FEAT, " not found!")}
              return(0)
            }
          }
        ) %>%
          unlist() %>%
          max()

        return(c(0, out.max))
      },
      SIMPLIFY = F,
      FEAT = features,
      ASS = assay,
      SLOT = slot
    )

    scatter.out <- lapply(
      seu.list,
      FUN=function(VIS){

        # tmp.nonzeros <- GetAssayData(
        #   VIS,
        #   assay=assay,
        # )[c(features[1],features[2]),]%>%
        #   apply(
        #     MARGIN=2,
        #     FUN=function(X) data.frame(
        #       FEAT1 = X[1]>0 & X[2]==0,
        #       FEAT2 = X[1]==0 & X[2]>0,
        #       DP = X[1]>0 & X[2]>0
        #     )
        #   )%>%
        #   do.call(what=rbind)

        # tmp.label <- paste0(
        #   "Notch1+: ",table(tmp.nonzeros$NOTCH)["TRUE"],"/",sum(tmp.nonzeros),"\n",
        #   "Gm13568+:",table(tmp.nonzeros$GM)["TRUE"],"/",sum(tmp.nonzeros),"\n",
        #   "Double Pos: ",table(tmp.nonzeros$DP)["TRUE"],"/",sum(tmp.nonzeros)
        # )

        # FeatureScatter(
        #   VIS,
        #   shuffle = T,
        #   jitter = F,
        #   feature1 = features[1],
        #   feature2 = features[2],
        #   plot.cor = F,
        #   pt.size = pt.size,
        #   cols=c("#c42a2e", "#d8a837", "#2e6198"),#TODO: parameterize
        #   group.by = scatter.group.by
        # )+
        # geom_label(
        #   x=1,
        #   y=5.2,
        #   size = small.font/ggplot2::.pt,
        #   # size = small.font/4,
        #   label=tmp.label
        # )+
        #
        # tmp.df <- cbind(
        #   VIS@meta.data,
        #   t(
        #     GetAssayData(
        #       VIS,
        #       assay=assay[1],
        #       slot=slot[1]
        #     )[features[1],]
        #   ),
        #   t(
        #     GetAssayData(
        #       VIS,
        #       assay=assay[2],
        #       slot=slot[2]
        #     )[features[2],]
        #   )
        # )
        #
        # #Plot
        # ggplot(
        #   tmp.df[sample(rownames(tmp.df)),],
        #   aes_string(
        #     x=features[1],
        #     y=features[2]
        #   )
        # )+
        #   geom_point(
        #     alpha=0.7
        #   )+

        #For undetected features...
        for(i in 1:length(features)){
          if(!features[i] %in% Features(VIS,assay=assay[i])){
            features[i] <- paste0("tmp_",features[i])# in case gene starts with a number...
            VIS@meta.data[[features[i]]]<- rep(0, length(Cells(VIS)))
          }
        }

        tmp.plot <- FeatureScatter(
          VIS,
          shuffle = T,
          jitter = F,
          feature1 = features[1],
          feature2 = features[2],
          plot.cor = T,
          pt.size = pt.size,
          cols = "black",
          # cols=c("#c42a2e", "#d8a837", "#2e6198"),#TODO: parameterize
          group.by = scatter.group.by
        )+
          geom_smooth(
            method="lm",
            formula = "y ~ x",
            color="black"
          )+
          xlim(gene.lims[[features[1]]])+
          ylim(gene.lims[[features[2]]])+
          labs(
            x=stringr::str_remove(features[1],pattern="tmp_"),
            y=stringr::str_remove(features[2],pattern="tmp_")
          )+
          scale_color_manual(
            values = c("#c42a2e", "#d8a837", "#2e6198"),#TODO: parameterize
          )+
          guides(
            color = guide_legend(override.aes = list(size=2))
          )+
          scatter.theme+
          theme(
            axis.title.y = element_text(
              face="bold.italic",
              hjust=0.5,
              size=font.size
            ),
            axis.title.x = element_blank(),
            legend.position = "right",
            plot.margin = unit(rep(0,4),"cm")
          )+
          coord_fixed(
            ratio=max(gene.lims[[features[1]]])/max(gene.lims[[features[2]]])
          )

        tmp.plot$layers[[1]]$aes_params$alpha <- 0.5

        return(tmp.plot)
      }
    )

    scatter.out[[length(scatter.out)]] <- scatter.out[[length(scatter.out)]]+
      theme(
        axis.title.x = element_text(
          face="bold.italic",
          hjust=0.5,
          size=font.size
        )
      )

    scatter.out <- wrap_plots(
      scatter.out,
      ncol=1,
      heights=rep(1,length(scatter.out)),
      guides="collect"
    )&theme(
      legend.position="none"#TODO
    )

    return(
      # scatter.out
      wrap_plots(
        maps.out,
        scatter.out,
        nrow=1,
        widths=c(3,1.5)
      )
    )
  }else{
    return(maps.out)
  }
}
