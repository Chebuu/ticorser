## gene ontology heatmap
setwd("/Volumes/Elements/time_course/")

cerv.thor.go.mf<-ReadDataFrameFromTsv("cervical_thoracic/cervical_thoracic_results/results_t/pathway/cerv_thor_DAVID_GO_MF_ALL.txt", row.names.col = NULL)
interested.go <- cerv.thor.go.mf[1:2,]


protein.binding <- interested.go[1,]
glycosaminoglycan.binding <- interested.go[2,]


cerv.thor.counts <- ReadDataFrameFromTsv("cervical_thoracic/cervical_thoracic_results/results_t/counts/cervical_thoracic__filtered_counts__Proportion__normalized__ruvg.tsv")


protein.binding.genes <- unlist(strsplit(as.vector(protein.binding$Genes), split = ", "))

protein.binding.counts <- cerv.thor.counts[which(rownames(cerv.thor.counts) %in% protein.binding.genes),]

cerv.thor.design <- ReadDataFrameFromTsv("cervical_thoracic/design_file/cervical_thoracic_design_file.txt")

glycosaminoglycan.binding.genes <- unlist(strsplit(as.vector(glycosaminoglycan.binding$Genes), split = ", "))
glycosaminoglycan.binding.counts <- cerv.thor.counts[which(rownames(cerv.thor.counts) %in% glycosaminoglycan.binding.genes),]

c.fc <- CalculateLogFoldChangeForEachTime(glycosaminoglycan.binding.counts, cerv.thor.design)

c.fc <- CalculateLogFoldChangeForEachTime(protein.binding.counts, cerv.thor.design)

CalculateLogFoldChangeForEachTime <- function(subset.counts, design.matrix) {
  times <- unique(design.matrix$Times)
  
  c.fc <- c()
  i <- 1
  for(time in times) {
    #print(time)
    sampl.names <- rownames(design.matrix)[which(design.matrix$Times == time )]
    sub.cols <- subset.counts[,which(colnames(subset.counts) %in% sampl.names)]
    t1 <- apply(sub.cols[,1:5], 1, mean ) 
    t2 <- apply(sub.cols[,6:10], 1, mean )
    fc <- log(t1/t2)
    c.fc <- cbind(c.fc, fc)
    colnames(c.fc)[i] <- time
    i<- i+1
  }
  return(c.fc)
} 

# times <- unique(cerv.thor.design$Times)
# 
# c.fc <- c()
# i <- 1
# for(time in times) {
#   print(time)
#   sampl.names <- rownames(cerv.thor.design)[which(cerv.thor.design$Times == time )]
#   sub.cols <- protein.binding.counts[,which(colnames(protein.binding.counts) %in% sampl.names)]
#   t1 <- apply(sub.cols[,1:5], 1, mean ) 
#   t2 <- apply(sub.cols[,6:10], 1, mean )
#   fc <- log(t1/t2)
#   c.fc <- cbind(c.fc, fc)
#   colnames(c.fc)[i] <- time
#   i<- i+1
# }
# 
# png("./GO-MF_Protein_Binding_Heatmap_dendrogram_gene_symbols_noscale.png",    # create PNG for the heat map        
#     width = 5*300,        # 5 x 300 pixels
#     height = 5*300,
#     res = 300,            # 300 pixels per inch
#     pointsize = 8) 
# 
# 
# 
# lmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
# # lmat = rbind(c(3,0),c(1,2),c(4,0))
# lwid = c(0.5, 4, 0.5)
# lhei = c(1,4,1)
# # lmat=rbind(c(2),c(3),c(1),c(4))
# # lhei=c(1,1,9,0)
# # lwid=c(1)
# 
# hm<-heatmap.2(c.fc, col = pal, scale="none",
#           dendrogram = "row", ## manage dendrogram
#           density.info = "none", ## manage density in legend
#           trace="none", ## manage trace along columns
#           margins=c(6, 6), ## manage margins
#           Colv=NA, ## manage clustering on columns
#           main="GO-MF Protein binding"
#           ,lmat=lmat
#           , lhei=lhei
#           , lwid = lwid
#           , cexRow = 0.45
#           # , key=FALSE
#           )
# dev.off()

PlotHeatmap <- function(fold.changes, scale="row", dendrogram="none", title="heatmap", row.labels.percent=0.5, save.png=FALSE, file.name=NULL) {
  if(save.png ) {
    if(!is.null(file.name)){
      # pdf(file.name,    # create PNG for the heat map        
      #     width = 5*300,        # 5 x 300 pixels
      #     height = 5*300)       
      png(file.name,    # create PNG for the heat map
          width = 5*300,        # 5 x 300 pixels
          height = 5*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8)
    } else {
      error("Please set a file.name for the png file!")
    }
  }
  require("gplots")
  pal <- colorRampPalette(c("red", "black", "green"))(n = 1000)
  
  lmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
  lwid = c(0.5, 4, 0.5)
  lhei = c(1,4,1)
  
  heatmap.2(
            fold.changes
            , col = pal
            , scale=scale, 
            dendrogram = dendrogram, ## manage dendrogram
            density.info = "none", ## manage density in legend
            trace="none", ## manage trace along columns
            margins=c(6, 6), ## manage margins
            Colv=NA, ## manage clustering on columns
            main=title
            ,lmat=lmat
            , lhei=lhei
            , lwid = lwid
            , cexRow = row.labels.percent
            # , key=FALSE
  )
  
  if(save.png) {
    dev.off()
  }
}

# 
# gage::geneData(genes = rownames(protein.binding.counts), exprs = log(protein.binding.counts + 1), ref = control, samp = treated, outname = outname, txt = TRUE, heatmap = TRUE, 
#                cexRow=2.2, cexCol=2.2, margins = c(10, 10), limit = 3, scatterplot = TRUE, dendrogram="none", 
#                pdf.size = c(21,30), key=FALSE)


# 
# require("biomaRt")
# listMarts()
# ensembl=useMart("ensembl")

ens.gene.map<-ReadDataFrameFromTsv("downloaded_references/conversion_maps/rnor5.0_biomart_ensembls_associatedgenenames.txt", row.names.col = NULL)
head(ens.gene.map)

glycosaminoglycan.binding.genes.o <- sort(glycosaminoglycan.binding.genes)
ens.gene.map.o <- ens.gene.map[order(ens.gene.map$Ensembl.Gene.ID),]
mini.map<-ens.gene.map.o[which(ens.gene.map.o$Ensembl.Gene.ID %in% glycosaminoglycan.binding.genes.o),]

glycosaminoglycan.binding.counts.o <- glycosaminoglycan.binding.counts[order(rownames(glycosaminoglycan.binding.counts)),]
rownames(glycosaminoglycan.binding.counts.o) <- mini.map$Associated.Gene.Name
c.fc <- CalculateLogFoldChangeForEachTime(glycosaminoglycan.binding.counts.o, cerv.thor.design)
PlotHeatmap(c.fc, row.labels.percent = 1, save.png = TRUE, file.name = "Glycosaminoglycan_heatmap_gene_symbols.png", title = "GO-MF Glycosaminoglycan binding")

####


protein.binding.genes.o <- sort(protein.binding.genes)
ens.gene.map.o <- ens.gene.map[order(ens.gene.map$Ensembl.Gene.ID),]
mini.map<-ens.gene.map.o[which(ens.gene.map.o$Ensembl.Gene.ID %in% protein.binding.genes.o),]
dim(mini.map)
length(protein.binding.genes.o)
protein.binding.counts.o <- protein.binding.counts[order(rownames(protein.binding.counts)),]
rownames(protein.binding.counts.o) <- mini.map$Associated.Gene.Name
c.fc <- CalculateLogFoldChangeForEachTime(protein.binding.counts.o, cerv.thor.design)

PlotHeatmap(c.fc, row.labels.percent = 0.45, save.png = TRUE, file.name = "Proteinbinding_heatmap_gene_symbols.png", title = "GO-MF Protein binding")
WriteDataFrameAsTsv(rownames(c.fc), "Proteinbinding_heatmap_gene_symbols")
##############################################################################################################################
### thoracic
CalculateLogFoldChangeForEachTimeConditions <- function(subset.counts, design.matrix) {
  times <- unique(design.matrix$Times)
  conditions <- unique(design.matrix$Conditions)
  c.fc <- c()
  i <- 1
  for(time in times) {
    #print(time)
    t <- c()
    sampl.names <- rownames(design.matrix)[which(design.matrix$Times == time )]
    sub.cols <- subset.counts[,which(colnames(subset.counts) %in% sampl.names)]
    for(condition in conditions) {
      sampl.cond.names <- rownames(design.matrix)[which(design.matrix$Conditions == condition)]
      cond.cols <- sub.cols[,which(colnames(sub.cols) %in% sampl.cond.names), drop=FALSE]
      t <- cbind(t, apply(cond.cols, 1, mean ) )
    }
    fc <- log(t[,1]/t[,2])
    c.fc <- cbind(c.fc, fc)
    colnames(c.fc)[i] <- time
    i<- i+1
  }
  return(c.fc)
} 

thor.counts <- ReadDataFrameFromTsv("Thoracic/Thoracic_results/counts/FeatureCounts/proportion/ruvg/thoracic__filtered_counts__Proportion__normalized__ruvg.tsv")

thor.go <- ReadDataFrameFromTsv("Thoracic/Thoracic_results/pathway/proportion/ruvg/thoracic_proportion_ruvg_david_GO_ALL.txt", row.names.col = NULL)
head(thor.go)
thor.prot.bind <- thor.go[18,]


thor.prot.bind.genes <- unlist(strsplit(as.vector(thor.prot.bind$Genes), split = ", "))

thor.count.prot.bind <- subset(thor.counts, rownames(thor.counts) %in% thor.prot.bind.genes)

thor.design <- ReadDataFrameFromTsv("Thoracic/design_file/thoracic_design_file.txt")


thor.prot.bind.fc <- CalculateLogFoldChangeForEachTimeConditions(thor.count.prot.bind, thor.design)

thor.prot.bind.fc.gs <- ConvertRowNamesToGeneSymbols(thor.prot.bind.fc, ens.gene.map)
  

PlotHeatmap(thor.prot.bind.fc.gs, row.labels.percent = 0.38, save.png = TRUE, file.name = "Thoracic_Proteinbinding_heatmap_gene_symbols.png", title = "Thoracic GO-MF Protein binding")
WriteDataFrameAsTsv(rownames(thor.prot.bind.fc.gs), "Thoracic_Proteinbinding_heatmap_gene_symbols")

ConvertRowNamesToGeneSymbols <- function(data.frame, map) {
  
  data.frame.o <- data.frame[order(rownames(data.frame)),]
  map.o <- map[order(map$Ensembl.Gene.ID),]
  
  mini.map<-map.o[which(map.o$Ensembl.Gene.ID %in% rownames(data.frame.o)),]
  
  rownames(data.frame.o) <- mini.map$Associated.Gene.Name
  return(data.frame.o)
}


##############################################################################################################################
### cervical


cerv.counts <- ReadDataFrameFromTsv("Cervical/Cervical_results/counts/FeatureCounts/proportion/ruvg/cervical__filtered_counts__Proportion__normalized__ruvg.tsv")

cervical.go <- ReadDataFrameFromTsv("Cervical/Cervical_results/pathway/proportion/ruvg/cervical_proportion_ruvg_david_GO_ALL_ex.txt", row.names.col = NULL)

head(cervical.go)
cerv.prot.bind <- cervical.go[8,]

cerv.prot.bind.genes <- unlist(strsplit(as.vector(cerv.prot.bind$Genes), split = ", "))

cerv.count.prot.bind <- subset(cerv.counts, rownames(cerv.counts) %in% cerv.prot.bind.genes)

cerv.design <- ReadDataFrameFromTsv("Cervical/design_file/cervical_design_file.txt")


cerv.prot.bind.fc <- CalculateLogFoldChangeForEachTimeConditions(cerv.count.prot.bind, cerv.design)

cerv.prot.bind.fc.gs <- ConvertRowNamesToGeneSymbols(cerv.prot.bind.fc, ens.gene.map)


PlotHeatmap(cerv.prot.bind.fc.gs, row.labels.percent = 0.2, save.png = TRUE, file.name = "Cervical_Proteinbinding_heatmap_gene_symbols.pdf", title = "Cervical GO-MF Protein binding")

WriteDataFrameAsTsv(rownames(cerv.prot.bind.fc.gs), "Cervical_Proteinbinding_heatmap_gene_symbols")


######
cerv.glycosaminoglycan <- cervical.go[351,]

cerv.glyco.genes <- unlist(strsplit(as.vector(cerv.glycosaminoglycan$Genes), split = ", "))

cerv.count.glyco <- subset(cerv.counts, rownames(cerv.counts) %in% cerv.glyco.genes)


cerv.glyco.fc <- CalculateLogFoldChangeForEachTimeConditions(cerv.count.glyco, cerv.design)

cerv.glyco.fc.gs <- ConvertRowNamesToGeneSymbols(cerv.glyco.fc, ens.gene.map)


PlotHeatmap(cerv.glyco.fc.gs, row.labels.percent = 1, save.png = TRUE, file.name = "Cervical_glycoaminoglycan_heatmap_gene_symbols.png", title = "Cervical GO-MF Glycoaminoglycan binding")

WriteDataFrameAsTsv(rownames(cerv.glyco.fc.gs), "Cervical_Proteinbinding_heatmap_gene_symbols")

######
thor.gam <- thor.go[80,]
thor.glyco.genes <- unlist(strsplit(as.vector(thor.gam$Genes), split = ", "))

thor.count.glyco <- subset(thor.counts, rownames(thor.counts) %in% thor.glyco.genes)

thor.glyco.fc <- CalculateLogFoldChangeForEachTimeConditions(thor.count.glyco, thor.design)

thor.glyco.fc.gs <- ConvertRowNamesToGeneSymbols(thor.glyco.fc, ens.gene.map)


PlotHeatmap(thor.glyco.fc.gs, row.labels.percent = 1, save.png = TRUE, file.name = "Thoracic_glycoaminoglycan_heatmap_gene_symbols.png", title = "Thoracic GO-MF Glycoaminoglycan binding")

WriteDataFrameAsTsv(rownames(thor.glyco.fc.gs), "Thoracic_Proteinbinding_heatmap_gene_symbols")
####
# thor immune response
thor.go.imm.resp <- thor.go[5,]


thor.imm.res.genes <- unlist(strsplit(as.vector(thor.go.imm.resp$Genes), split = ", "))

thor.count.imm.resp <- subset(thor.counts, rownames(thor.counts) %in% thor.imm.res.genes)

thor.i.r.fc <- CalculateLogFoldChangeForEachTimeConditions(thor.count.imm.resp, thor.design)

thor.i.r.fc.gs <- ConvertRowNamesToGeneSymbols(thor.i.r.fc, ens.gene.map)


PlotHeatmap(thor.i.r.fc.gs, row.labels.percent = 1, save.png = TRUE, file.name = "Thoracic_immune_response_heatmap_gene_symbols.png", title = "Thoracic GO-BP Immune Response")
WriteDataFrameAsTsv(rownames(thor.i.r.fc.gs), "Thoracic_Immune_response_heatmap_gene_symbols")




