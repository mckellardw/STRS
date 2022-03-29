
# Skeletal muscle samples ----
mirge3.dir = "/workdir/dwm269/totalRNA/data/sp_mirge3/"
mirge.list <- list()
for(i in 1:length(skm.list)){ #
  cat(paste0("Concatenating mirge outs for ", meta_skm$sample[i],"...\n"))
  barcodes = Cells(skm.list[[i]])
  tmp.dir = paste0(mirge3.dir,meta_skm$capture.area[i])
  
  tmp.list <- lapply(
    barcodes,
    FUN=function(CB){
      tmp.file = paste0(tmp.dir,"/",CB,"/min10/miR.Counts.csv")
      
      if(file.exists(tmp.file)){
        X = read.csv(tmp.file)
        
        out = matrix(X$subset_R2, ncol=1)
        rownames(out) <- X$miRNA
        colnames(out) <- CB
        
        #reorder for easy merging
        out <- out[order(rownames(out)),]
        
        return(out)
      }else{
        # message(paste0("File for ", CB, " not found..."))
        
        cat("|")
        out <- rep(0,1927)
        
        return(out)
      }
    }
  )
  
  # Merge into a sparse matrix
  if(!is.null(tmp.list)){
    mirge.list[[i]] <- do.call(cbind, tmp.list)%>%
      as.sparse()
    colnames(mirge.list[[i]]) <- barcodes
  }
  
  # Save the concatenated matrix to disk
  # write_sparse(
  #   x = mirge.list[[i]],
  #   path = paste0(meta_skm$data.dir.kallisto[i],"/mirge3_matrix/",meta_skm$sample[i])
  # )
}

# Heart samples ----
mirge3.dir = "/workdir/dwm269/totalRNA/data/sp_mirge3/"
mirge.list <- list()
for(i in 1:length(heart.list)){ #
  cat(paste0("Concatenating mirge outs for ", meta_heart$sample[i],"...\n"))
  barcodes = Cells(heart.list[[i]])
  # barcodes = read.csv("/home/dwm269/DWM_utils/align_pipes/10x_kallisto/resources/barcodes_10x/visium-v1.txt",header = F)
  tmp.dir = paste0(mirge3.dir,meta_heart$capture.area[i])
  
  tmp.list <- lapply(
    barcodes,
    FUN=function(CB){
      tmp.file = paste0(tmp.dir,"/",CB,"/miR.Counts.csv")
      
      if(file.exists(tmp.file)){
        X = read.csv(tmp.file)
        
        out = matrix(X$subset_R2, ncol=1)
        rownames(out) <- X$miRNA
        colnames(out) <- CB
        
        #reorder for easy merging
        out <- out[order(rownames(out)),]
        
        return(out)
      }else{
        # message(paste0("File for ", CB, " not found..."))
        
        cat("|")
        out <- rep(0,1927)
        
        return(out)
      }
    }
  )
  
  # Merge into a sparse matrix
  if(!is.null(tmp.list)){
    mirge.list[[i]] <- do.call(cbind, tmp.list)%>%
      as.sparse()
    colnames(mirge.list[[i]]) <- barcodes
  }
  
  # Save the concatenated matrix to disk
  # write_sparse(
  #   x = mirge.list[[i]],
  #   overwrite = T,
  #   path = paste0(meta_heart$data.dir.kallisto.REO[i],"/mirge3_matrix/",meta_heart$sample[i])
  # )
}

  
# Whole capture areas ----
mirge3.dir = "/workdir/dwm269/totalRNA/data/sp_mirge3/"
mirge.list <- list()
for(i in 1:length(whole.list)){ #
  cat(paste0("Concatenating mirge outs for ", meta_whole$sample[i],"...\n"))
  barcodes = read.csv("~/DWM_utils/align_pipes/10x_kallisto/resources/barcodes_10x/visium-v1.txt")%>%unlist()
  tmp.dir = paste0(mirge3.dir,meta_whole$capture.area[i])
  
  tmp.list <- lapply(
    barcodes,
    FUN=function(CB){
      tmp.file = paste0(tmp.dir,"/",CB,"/min10/miR.Counts.csv")
      
      if(file.exists(tmp.file)){
        X = read.csv(tmp.file)
        
        out = matrix(X$subset_R2, ncol=1)
        rownames(out) <- X$miRNA
        colnames(out) <- CB
        
        #reorder for easy merging
        out <- out[order(rownames(out)),]
        
        return(out)
      }else{
        # message(paste0("File for ", CB, " not found..."))
        
        cat("|")
        out <- rep(0,1927)
        
        return(out)
      }
    }
  )
  
  # Merge into a sparse matrix
  if(!is.null(tmp.list)){
    mirge.list[[i]] <- do.call(cbind, tmp.list)%>%
      as.sparse()
    colnames(mirge.list[[i]]) <- barcodes
  }
  
  # Save the concatenated matrix to disk
  # write_sparse(
  #   x = mirge.list[[i]],
  #   path = paste0(meta_skm$data.dir.kallisto[i],"/mirge3_matrix/",meta_whole$sample[i])
  # )
}


