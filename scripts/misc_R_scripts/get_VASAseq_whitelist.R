
# Build VASAdrop whitelist
vasa.ids <- read.csv("spTotal/resources/metadata_sheets/VASAdrop_samples.txt")

DATADIR = "/workdir/dwm269/totalRNA/VASAseq_2022/GEO/"

id.list <- list()
for(ID in unlist(vasa.ids)){
  id.list[[ID]] <- read.csv(
    paste0(DATADIR,ID,"_total.UFICounts.tsv"), 
    sep="\t",#header = F,
    nrows = 1
  ) %>% colnames()
}

whitelist <- unlist(id.list) %>% sort() %>% unique()

# Remove barcodes with 4 or more G's or A's in a row
for(i in 4:16){
  whitelist = whitelist[!grepl(
    pattern=paste0(rep("A",i),collapse = ""),
    whitelist
  )]
  
  whitelist = whitelist[!grepl(
    pattern=paste0(rep("G",i),collapse = ""),
    whitelist
  )]
}

write.table(
  whitelist,
  file="/home/dwm269/DWM_utils/align_pipes/VASAseq_kallisto/resources/VASAdrop_custom_whitelist.txt",
  quote = F,
  col.names=F,
  row.names = F
)
