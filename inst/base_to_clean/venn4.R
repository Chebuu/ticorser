Venn4de <- function(x,y,z,w,label1,label2,label3,label4, title="Venn Diagram", intersection.flag=TRUE, plot.dir=NULL, enrich.lists.flag=FALSE){
  
  require(VennDiagram)
  require(tiff)
  # print("capabilities(tiff)")
  # print(capabilities("tiff"))
  
  a = x
  b = y
  c = z
  d = w
  Lists <- list(a,b,c,d)
  
  if(is.null(plot.dir)) {
    stop("Please provide a directory where to save VENN results!")
  }
  names(Lists) <- c(label1, label2, label3, label4)
  ab <- intersect(a,b)     #common gene names
  cd <- intersect(c,d)
  abcd <- intersect(ab,cd)
  
  plot.combined.string <- paste(label1,"_",label2,"_",label3,"_", label4, sep="")
  
  if(intersection.flag) {
    if(is.null(plot.dir)) {
      stop("Please provide a directory where to save VENN results!")
    }
    
    res.path <- file.path(plot.dir, "Venn4", paste0(plot.combined.string,"_gene_lists"))
    dir.create(res.path, recursive = TRUE, showWarnings = FALSE)
    
    outputName=paste(plot.combined.string,"_genes_in_intersection.txt",sep="")
    filepathname <- file.path(res.path, filename)
    write.table(abcd, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
  }
  
  outputName=paste(label1,"_",label2,"_",label3,"_",label4,"_VennDiagramDE.tiff",sep="")
  file.path.name <- file.path(plot.dir, "Venn4")
  dir.create(file.path.name, recursive = TRUE, showWarnings = FALSE)
  file.path.name <- paste0(file.path.name, "/", outputName)
  
  # filepathname=paste(file.path.name,outputName,sep="")
  
  print(venn.diagram(x=Lists, filename = file.path.name,  col = "transparent",fill=c("blue","red","green","pink"), alpha = 0.50, cex = 1.5, fontfamily ="serif", fontface = "bold", cat.cex = 1.2, cat.pos =0,cat.dist = 0.07,cat.fontfamily = "serif", rotation.degree = 270, margin = 0.2))
  
  message("The ", outputName, " file has been saved in the ", file.path.name, " folder!", sep="")
   
}
