## biomart
# http://may2012.archive.ensembl.org

library("biomaRt")

mart = useMart('ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl', host="may2012.archive.ensembl.org")


attrs <- listAttributes(mart)

  attrs

  ens.gs.mm9 <- getBM(attributes = c("ensembl_gene_id", "external_gene_id"), mart = mart)
  
  WriteDataFrameAsTsv(data.frame.to.save = ens.gs.mm9, file.name.path = "references/ensembl_genesymbol_biomart_mm9")
  
  
  
  library("biomaRt")
  
  mart = useMart('ENSEMBL_MART_ENSEMBL')
  mart <- useDataset(dataset="mmusculus_gene_ensembl", mart = mart)
  
  attrs <- listAttributes(mart)
  
  attrs
  
  ens.gs.mm10 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = mart)
  
  WriteDataFrameAsTsv(data.frame.to.save = ens.gs.mm10, file.name.path = "references/ensembl_genesymbol_biomart_mm10")
  

   
  

