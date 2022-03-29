
library(dplyr)
library(parallel)

read.counts = c(
  100000,
  250000,
  500000,
  1000000,
  2500000,
  5000000,
  10000000,
  25000000,
  50000000
)

seeds <- c(
  3,
  42,
  100,
  142
)

# Get all permutations of seeds and read counts
run.permuts <- expand.grid(seeds, read.counts)%>%
  apply(
    MARGIN = 1, 
    function(X) data.frame(sd=X[1], rd=as.integer(X[2]))
    # simplify = T
  )

# Prefix for R1 and R2 fqs
fqs.in <- c(
  "/workdir/dwm269/muscle_visium/data/geo_vis579/Vis5A/Vis5A_S2_L001_",
  "/workdir/dwm269/muscle_visium/data/geo_vis579/Vis7B/Vis7B_S2_L001_",
  "/workdir/dwm269/muscle_visium/data/geo_vis579/Vis9A/Vis9A_S2_L001_",
  "/workdir/dwm269/totalRNA/data/fastq/Vis_CTRL_Heart_MM/spatialRNAseq_T1L_Heart_D7PI_S1_L001_",
  "/workdir/dwm269/totalRNA/data/fastq/Vis_yPAP3_2K/Vis_yPAP_3A_S1_L001_",
  "/workdir/dwm269/totalRNA/data/fastq/Vis_yPAP3_2K/Vis_yPAP_3B_S2_L001_",
  "/workdir/dwm269/totalRNA/data/fastq/Vis_yPAP3_2K/Vis_yPAP_3C_S3_L001_",
  "/workdir/dwm269/totalRNA/data/fastq/Vis_yPAP3_2K/Vis_yPAP_3D_S4_L001_",
  "/workdir/dwm269/totalRNA/data/fastq/Vis_CTRL_Heart_MM/spatialRNAseq_Mock_Heart_D7PI_S1_L001_"
)

fqs.out <- c(
  "/workdir/dwm269/totalRNA/data/fastq/strs_subsample/Vis5A_",
  "/workdir/dwm269/totalRNA/data/fastq/strs_subsample/Vis7B_",
  "/workdir/dwm269/totalRNA/data/fastq/strs_subsample/Vis9A_",
  "/workdir/dwm269/totalRNA/data/fastq/strs_subsample/ctrl_T1L_",
  "/workdir/dwm269/totalRNA/data/fastq/strs_subsample/3A_",
  "/workdir/dwm269/totalRNA/data/fastq/strs_subsample/3B_",
  "/workdir/dwm269/totalRNA/data/fastq/strs_subsample/3C_",
  "/workdir/dwm269/totalRNA/data/fastq/strs_subsample/3D_",
  "/workdir/dwm269/totalRNA/data/fastq/strs_subsample/ctrl_mock_"
)

out.X <- mclapply(
  FUN=function(X){
    for(i in 1:length(fqs.in)){
      if(!file.exists(paste0(fqs.out[i],X$sd,"_",X$rd,"_R1.fq")) & !file.exists(paste0(fqs.out[i],X$sd,"_",X$rd,"_R1.fq.gz"))){
        system(
          paste0(
            "seqtk sample -s", X$sd, " ",fqs.in[i],"R1_001.fastq.gz ",X$rd," > ", fqs.out[i],X$sd,"_",X$rd,"_R1.fq"
          )
        )
      }
      if(!file.exists(paste0(fqs.out[i],X$sd,"_",X$rd,"_R2.fq")) & !file.exists(paste0(fqs.out[i],X$sd,"_",X$rd,"_R2.fq.gz"))){
        system(
          paste0(
            "seqtk sample -s", X$sd, " ",fqs.in[i],"R2_001.fastq.gz ",X$rd," > ", fqs.out[i],X$sd,"_",X$rd,"_R2.fq"
          )
        )
      }
    }
    
    return(X)
  },
  X=run.permuts,
  mc.cores = 16
)
