# Code snippet to pull gene info from biomaRt

require(biomaRt)

if(!exists("mouse.info")){
  mouse = biomaRt::useMart(
    # biomart = "ENSEMBL_MART_MOUSE",
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "mmusculus_gene_ensembl"
    # host = "https://apr2020.archive.ensembl.org" ##for mm10 info
  )
  
  mouse.info <- biomaRt::getBM(
    attributes = c(
      "mgi_symbol", "gene_biotype","transcript_biotype", 
      # "description", "name_1006","definition_1006","external_gene_name",
      "ensembl_gene_id","ensembl_transcript_id","ensembl_peptide_id",
      "chromosome_name","transcript_length", "strand"
    ),
    uniqueRows = T,
    # host = "useast.ensembl.org",
    verbose = T,
    mart = mouse
  )
}

#Replace missing MGI symbols with ensembl IDs
for(i in 1:nrow(mouse.info)){
  if(mouse.info[i,"mgi_symbol"]==""){
    mouse.info[i,"mgi_symbol"] <- mouse.info[i,"ensembl_gene_id"]
  }
}