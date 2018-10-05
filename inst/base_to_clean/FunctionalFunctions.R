enrichKEGGFunction <- function(gene.list.to.enrich, gene.type = "ENSEMBL", organism.db = c("org.Rn.eg.db", "org.Mm.eg.db"), organism.code = c("rno", "mmu"), threshold=0.05, functional.folder, filename) {
  # require(assign("organism.db", organism.db))
  require("clusterProfiler")
  
  kegg.folder <- UpdateFolderPath(functional.folder, "clusterProfiler", "KEGG")
  filename <- UpdateFilename(filename, "KEGG_PATHWAYS")
  sign.genes.entrez <- bitr(gene.list.to.enrich, fromType = gene.type, toType = c("ENTREZID", "SYMBOL"), OrgDb = organism.db)
  # total.genes.entrez <- bitr(total.gene.list, fromType = gene.type, toType = c("ENTREZID", "SYMBOL"), OrgDb = organism.db)
  kegg.results <- enrichKEGG(gene = sign.genes.entrez$ENTREZID,
                             organism = organism.code,
                             # universe = total.genes.entrez$ENTREZID,
                             pvalueCutoff = threshold)
  ksdf <- as.data.frame(kegg.results)
  ksdf$geneID <- gsub("/", "; ", ksdf$geneID)
  WriteDataFrameAsTsv(data.frame.to.save = ksdf, file.name.path = file.path(kegg.folder, filename))
  GenerateAndSaveNetwork(kegg.results, kegg.folder, filename)
  return(ksdf)
}


enrichGOFunction <- function(gene.list.to.enrich, gene.type = "ENSEMBL", organism.db = c("org.Rn.eg.db", "org.Mm.eg.db"), organism.code = c("rno", "mmu"), threshold=0.05, functional.folder, filename, ontology =c("CC", "MF", "BP")) {
  # require(organism.db)
  require("clusterProfiler")
  
  go.folder <- UpdateFolderPath(functional.folder, "clusterProfiler", "GO")
  filename <- UpdateFilename(filename, paste0("clusterProfiler_GO_",ontology))
  sign.genes.entrez <- bitr(gene.list.to.enrich, fromType = gene.type, toType = c("ENTREZID", "SYMBOL"), OrgDb = organism.db)
  
  ego <- enrichGO(gene = sign.genes.entrez$ENTREZID,
                  keytype = "ENTREZID",
                  # universe      = rownames(thoracic.de.uqua.notnorm2.03d$de.not.na),
                  OrgDb         = organism.db,
                  ont           = ontology,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = threshold,
                  qvalueCutoff  = threshold,
                  readable      = TRUE)
  
  ksdf <- as.data.frame(ego)
  ksdf$geneID <- gsub("/", "; ", ksdf$geneID)
  WriteDataFrameAsTsv(data.frame.to.save = ksdf, file.name.path = file.path(go.folder, filename))
  GenerateAndSaveHierarchicalGO(go.results=ego, go.folder=go.folder, filename=filename, ontology=ontology)
  return(ksdf)
}



GenerateAndSaveNetwork <- function(kegg.results, kegg.folder, filename) {
  require("clusterProfiler")
  
  plot.functional.folder <- UpdateFolderPath(kegg.folder, "KEGG_MAP")
  filename <- UpdateFilename(filename, "kegg_map.pdf")
  filename <- file.path(plot.functional.folder, filename)
  
  pdf(filename)
  enrichMap(kegg.results)
  dev.off()
  message(filename, " Saved on disk!")
}



GenerateAndSaveHierarchicalGO <- function(go.results, go.folder, filename, ontology) {
  require("clusterProfiler")
  plot.functional.folder <- UpdateFolderPath(go.folder, paste0("GO_", ontology, "_Hierarchy"))
  filename <- UpdateFilename(filename, "_Hierarchy.pdf")
  filename <- file.path(plot.functional.folder, filename)
  
  pdf(filename)
  # tryCatch({
    plotGOgraph(go.results)
    # message("problem with plotgograph")
  # })
  dev.off()
  message(filename, " Saved on disk!")
}

########################################################################## GProfiler

# exclude.electr.flag <- FALSE, go.label <- "GO:BP", include.graph.flag <- TRUE, path.label <- "KEGG", organism.name = "rnorvegicus", significant.flag = TRUE, min.isect.size=0, correction.method = "fdr"

enrichPathwayGProfiler <- function(gene.list.to.enrich, organism.name = c("rnorvegicus", "mmusculus"), exclude.electr.flag = FALSE, path.label = c("KEGG", "REAC"), significant.flag = FALSE, min.isect.size = 0, correction.method = "fdr", functional.folder, filename) {
  # require(assign("organism.db", organism.db))
  require("gProfileR")
  
  path.folder <- UpdateFolderPath(functional.folder, "gProfiler", path.label)
  filename <- UpdateFilename(filename, path.label, "_PATHWAYS")
  
  ksdf <- gprofiler(query = gene.list.to.enrich, organism = organism.name, significant = significant.flag, min_isect_size = min.isect.size, correction_method = correction.method, src_filter = path.label, exclude_iea = exclude.electr.flag)
  ksdf$intersection <- gsub(",", "; ", ksdf$intersection)
  
  ksdf <- ksdf[order(ksdf$p.value), ]
  
  print(head(ksdf))
  
  if(dim(ksdf)[1]!=0) {
    if(significant.flag) {
      WriteDataFrameAsTsv(data.frame.to.save = ksdf, file.name.path = file.path(path.folder, paste0(filename, "_significant")))
    } else {
      WriteDataFrameAsTsv(data.frame.to.save = ksdf[which(ksdf$significant==TRUE),], file.name.path = file.path(path.folder, paste0(filename, "_significant")))
      WriteDataFrameAsTsv(data.frame.to.save = ksdf, file.name.path = file.path(path.folder, paste0(filename, "_all")))
    }
  } else {
    message("no results produced for ", filename)
  }
  
  
  return(ksdf)
}


enrichGOGProfiler <- function(gene.list.to.enrich, organism.name = c("rnorvegicus", "mmusculus"), exclude.electr.flag = FALSE, ontology = c("BP", "MF", "CC"),  significant.flag = FALSE, min.isect.size = 0, correction.method = "fdr", functional.folder, filename) {
  # require(organism.db)
  require("gProfileR")
  
  go.folder <- UpdateFolderPath(functional.folder, "gProfiler", "GO")
  filename <- UpdateFilename(filename, paste0("gProfiler_GO_", ontology))
  
  gprof.ontology <- paste0("GO:", ontology)
  
  ego <- gprofiler(query = gene.list.to.enrich, organism=organism.name, significant = significant.flag, min_isect_size = min.isect.size, correction_method = correction.method, src_filter = gprof.ontology, exclude_iea = exclude.electr.flag)
  
  ego$intersection <- gsub(",", "; ", ego$intersection)
  ego<-ego[order(ego$p.value),]
  print(head(ego))
  
  if(dim(ego)[1]!=0) {
    if(significant.flag) {
      WriteDataFrameAsTsv(data.frame.to.save = ego, file.name.path = file.path(go.folder, paste0(filename, "_significant")))
    } else {
      WriteDataFrameAsTsv(data.frame.to.save = ego[which(ego$significant==TRUE),], file.name.path = file.path(go.folder, paste0(filename, "_significant")))
      WriteDataFrameAsTsv(data.frame.to.save = ego, file.name.path = file.path(go.folder, paste0(filename, "_all")))
    }
  } else {
    message("no results produced for ", filename)
  }
  
  return(ego)
}

PlotKeggMapTimeCoursePathview <- function(interested.genes.kegg.lfc, kegg.id, specie, gene.idtype, output.dir, prefix, plot.col.key.flag=TRUE, is.signal.flag=TRUE, low.color.list = c("gene"="red", cpd="red"), mid.color.list= c("gene"="grey", cpd="grey"), high.color.list = c("gene"="green", cpd="green")) {
  
  require("pathview")
  pathview(gene.data = interested.genes.kegg.lfcfc, pathway.id = kegg.id, species = specie, 
           gene.idtype = gene.idtype, kegg.dir = output.dir, out.suffix=prefix, 
           plot.col.key=plot.col.key.flag,
           low = low.color.list,
           mid = mid.color.list,
           high = high.color.list,
           is.signal=is.signal.flag,
           # keys.align="y", 
           kegg.native=TRUE, muti.state=TRUE, same.layer=TRUE)  
}

PlotListKeggMapTimeCoursePathview <- function(interested.genes.kegg.lfc, kegg.id.list, specie, gene.idtype, output.dir, prefix, plot.col.key.flag=TRUE, is.signal.flag=TRUE, low.color.list = c("gene"="red", cpd="red"), mid.color.list= c("gene"="grey", cpd="grey"), high.color.list = c("gene"="green", cpd="green")) {
  for(kegg.id in kegg.id.list) {
    PlotKeggMapTimeCoursePathview(interested.genes.kegg.lfc, kegg.id, specie, gene.idtype, output.dir, prefix, plot.col.key.flag, is.signal.flag, low.color.list, mid.color.list, high.color.list)
  }
}
