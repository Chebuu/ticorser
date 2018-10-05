require("clusterProfiler")

  genes.sign <- rownames(cervical.de.uqua.notnorm2.07d.noi95$SIGN)
length(genes)
genes.entrez.sign <- bitr(genes.sign, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = "org.Rn.eg.db")
head(genes.entrez.sign)
dim(genes.entrez)

ego <- enrichGO(gene          = genes,
                keytype = "ENSEMBL",
                universe      = rownames(thoracic.de.uqua.notnorm2.03d$de.not.na),
                OrgDb         = org.Rn.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

plotGOgraph(ego)
enrichMap(kegg.results)

enrichKEGGFunction <- function(gene.list.to.enrich, gene.type = "ENSEMBL", organism.db = "org.Rn.eg.db", organism.code = "rno", threshold=0.05, functional.folder, filename) {
  require(organism.db)
  require("clusterProfiler")
  
  kegg.folder <- UpdateFolderPath(functional.folder, "clusterProfiler", "KEGG")
  filename <- UpdateFilename(filename, "KEGG_PATHWAYS")
  sign.genes.entrez <- bitr(genes.list.to.enrich, fromType = gene.type, toType = c("ENTREZID", "SYMBOL"), OrgDb = organism.db)
  # total.genes.entrez <- bitr(total.gene.list, fromType = gene.type, toType = c("ENTREZID", "SYMBOL"), OrgDb = organism.db)
  kegg.results <- enrichKEGG(gene = sign.genes.entrez$ENTREZID,
                   organism = organism.code,
                   # universe = total.genes.entrez$ENTREZID,
                   pvalueCutoff = threshold)
  ksdf <- as.data.frame(kegg.results)
  ksdf$geneID <- gsub("/", ", ", ksdf$geneID)
  WriteDataFrameAsTsv(data.frame.to.save = ksdf, file.name.path = filename)
  GenerateAndSaveNetwork(kegg.results, functional.folder, filename)
  return(kegg.results)
}

GenerateAndSaveNetwork<- function(kegg.results, functional.folder, filename) {
  require("clusterProfiler")
  plot.functional.folder <- UpdateFolderPath(functional.folder, "KEGG_MAP")
  filename <- UpdateFilename(filename, "kegg_map.pdf")
  pdf(filename)
  enrichMap(kegg.results)
  dev.off()
}


GenerateAndSaveHierarchicalGO <- function(go.results, functional.folder, filename, ontology) {
  require("clusterProfiler")
  plot.functional.folder <- UpdateFolderPath(functional.folder, paste0("GO_", ontology, "_Hierarchy"))
  filename <- UpdateFilename(filename, "GO_Hierarchy.pdf")
  pdf(filename)
  plotGOgraph(go.results)
  dev.off()
}


enrichGOFunction <- function(gene.list.to.enrich, gene.type = "ENSEMBL", organism.db = "org.Rn.eg.db", organism.code = "rno", threshold=0.05, functional.folder, filename, ontology =c("CC", "MF", "BP")) {
  require(organism.db)
  require("clusterProfiler")
  
  go.folder <- UpdateFolderPath(functional.folder, "clusterProfiler", "GO")
  filename <- UpdateFilename(filename, paste0("clusterProfiler_GO_",ontology))
  sign.genes.entrez <- bitr(genes.list.to.enrich, fromType = gene.type, toType = c("ENTREZID", "SYMBOL"), OrgDb = organism.db)
  
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
  ksdf$geneID <- gsub("/", ", ", ego$geneID)
  WriteDataFrameAsTsv(data.frame.to.save = ego, file.name.path = filename)
  GenerateAndSaveHierarchicalGO(ego, functional.folder, filename, ontology)
  return(ego)
}



# kk <- enrichKEGG(gene = genes.entrez.sign$ENTREZID,
#                  organism     = 'rno',
#                  pvalueCutoff = 0.05)
# head(kk)
# 


require("GGally")
require("plotly")
df <- normalized.counts.prop.uqua[, which(colnames(normalized.counts.prop.uqua)%in%rownames(cervical.14d.design.matrix)) ]
df <- normalized.counts.prop.uqua[, which(colnames(normalized.counts.prop.uqua)%in%rownames(cervical.design.matrix)[1:15]) ]
ddf <- data.frame()
for(j in 1:dim(df)[2]) {
  samples <- rep(colnames(df)[j], times=dim(df)[1])
  genes <- rownames(df)
  values <- df[,j]
  ddf1 <- cbind(genes, samples, values)
  ddf <- rbind(ddf, ddf1)
}

# bob <- data.frame(lapply(ddf, as.character), stringsAsFactors=FALSE)
# head(bob)
# head(bob[,1])vab
# bob[,3] <- as.numeric(bob[,3])

PlotScatterPlotMatrix(normalized.counts.prop.uqua, cervical.14d.design.matrix, plot.folder = "./", prefix.plot = "aaa", show.plot.flag = FALSE, plotly.flag = FALSE, save.plot = TRUE)
PlotTimesBoxplot(data.frame.to.plot = normalized.counts.prop.uqua, design.matrix = cervical.14d.design.matrix, output.path = "./", prefix.plot = "aaa", show.plot.flag = FALSE, save.plot = TRUE, plotly.flag = FALSE)
p<-ggpairs(ddf)#, mapping = aes(gene=rownames(df)))#, mapping = aes_string(colour="Species"))

pdf("plot.pdf")
plot(p)
dev.off()



# ddf <- data.frame()
# ddf <- stack(df)

plot(p)
ggplotly(p)

head(iris)

rno.kegg<-download_KEGG("hsa")
head(rno.kegg)

length(unique(rno.kegg[[1]][,1]))



data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = "org.Hs.eg.db")
head(gene.df)
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

require("refGenome")
beg <- ensemblGenome()
gtf.file <- "downloaded_references/ensembl_site/Rattus_norvegicus.Rnor_5.0.79.gtf"
read.gtf(object = beg, gtf.file)

beg@ev$gtf[which(beg@ev$gtf$gene_biotype != "protein_coding"),]

head(beg@ev$gtf)
tail(beg@ev$gtf)

gtf.df <- beg@ev$gtf[, c("seqid", "source", "feature", "start", "end", "score", "strand", "frame")]
head(gtf.df)

# attribute <- 
  
  gtf.table<-read.table(file = gtf.file, header = FALSE, sep = "\t", quote = "", skip = 5, stringsAsFactors = FALSE)

  # chrs <- unique(gtf.table$V1)
  # sort(chrs)
 head(gtf.table)

gtf.attributes <- strsplit(gtf.table$V9, split = "; ")
head(gtf.attributes)


# aa <- grep(pattern = "gene_biotype \"protein_coding\"", x = gtf.attributes) 
gtf.gsubbed1 <- lapply(X = gtf.attributes, gsub, pattern = ";", replacement = "")
head(gtf.gsubbed1)
gtf.gsubbed2 <- lapply(X = gtf.gsubbed1, gsub, pattern = "\"", replacement = "")
head(gtf.gsubbed2)

gtf.gsubbed3 <- lapply(X = gtf.gsubbed2, gsub, pattern = " ", replacement = ":")
head(gtf.gsubbed3)

aa <- grep(pattern = "gene_biotype:protein_coding", x = gtf.gsubbed3) 
length(aa)
# gtf.gsubbed3[aa]
# dim(gtf.table)
# gtf.attributes[aa]

gtf.subsetted<-gtf.table[aa,]
tail(gtf.subsetted)

canonical.chromosomes <- c(1:20, "X")

dim(gtf.subsetted)

ind.row.canonical <- which(gtf.subsetted$V1 %in% canonical.chromosomes)
gtf.subsetted.canonical <- gtf.subsetted[ind.row.canonical,]
dim(gtf.subsetted.canonical)
sort(unique(gtf.subsetted.canonical$V1))


write.table(x = gtf.subsetted.canonical, file = "downloaded_references/ensembl_site/Rattus_norvegicus.Rnor_5.0.79_gene_biotype_protein_coding_canonical_chrs_NO_MT.gtf", sep = "\t", col.names = FALSE, row.names = FALSE, qmethod = "double", quote = FALSE)
gtf.file <- "downloaded_references/ensembl_site/Rattus_norvegicus.Rnor_5.0.79_gene_biotype_protein_coding_canonical_chrs_NO_MT.gtf"
gtf.table<-read.table(file = gtf.file, header = FALSE, sep = "\t", quote = "", skip = 5, stringsAsFactors = FALSE)
gtf.table$V1 <- paste0("chr",gtf.table$V1)
head(gtf.table)
gtf.table.u <- unique(gtf.table[,c(1,4,5)])
head(gtf.table[,c(1,4,5)])
head(gtf.table.u)
write.table(x = gtf.table.u, file = "downloaded_references/ensembl_site/rnor5_subsetted_for_coverage_unique.bed", sep = "\t", col.names = FALSE, row.names = FALSE, qmethod = "double", quote = FALSE)


bed.graph<-read.table("cervical_thoracic/all_bams/cov/cervical_03_cov.bedgraph", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)

bed.graph$V1[which(bed.graph$V1 %in% "MT")] <- "M"
head(bed.graph)
bed.graph$V1 <- paste0("chr", bed.graph$V1)

write.table(bed.graph, file="cervical_thoracic/all_bams/cov/cervical_03_cov1.bedgraph", sep = "\t", col.names = FALSE, row.names = FALSE, qmethod = "double", quote = FALSE)

bed.file<-read.table("cervical_thoracic/all_bams/cov/aa.bed", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
bed.file$V1[which(bed.file$V1 %in% "MT")] <- "M"
head(bed.file)
bed.file$V1 <- paste0("chr", bed.file$V1)
write.table(bed.file, file="cervical_thoracic/all_bams/cov/aa_no_N_chr.bed", sep = "\t", col.names = FALSE, row.names = FALSE, qmethod = "double", quote = FALSE)

arr.ind.40k <- which(normalized.counts.prop.uqua >40000, arr.ind = TRUE)
head(arr.ind.40k)
dim(arr.ind.40k)
sub.norm.prop.uqua.40k <- normalized.counts.prop.uqua[unique(arr.ind.40k[,1]), ]
sub.norm.prop.uqua.40k <- sub.norm.prop.uqua.40k[order(rownames(sub.norm.prop.uqua.40k)),]

over40k.ens.genes<-unique(rownames(normalized.counts.prop.uqua)[arr.ind.40k[,1]])

over40k.ens.symb.genes <- ens.symb.biomart.map[which(ens.symb.biomart.map$Ensembl.Gene.ID %in% over40k.ens.genes),]
over40k.ens.symb.genes <- over40k.ens.symb.genes[order(over40k.ens.symb.genes[,1]),]

sub.norm.prop.uqua.40k$symb <- over40k.ens.symb.genes$Associated.Gene.Name




install.packages("gProfileR")
require(gProfileR)
?gprofiler


query.names <- rownames(cerv.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN)

query.names.gs <- ens.symb.biomart.map[which(ens.symb.biomart.map$Ensembl.Gene.ID %in% query.names),]


mmm<-gconvert(query = query.names, organism = "rnorvegicus", target = "ENSG")

prova<-
head(prova)





exclude.electr.flag <- FALSE ##include IEA cause they are too much
go.label <- "GO:BP"
include.graph.flag <- TRUE
path.label <- "KEGG"
organism.name = "rnorvegicus"
significant.flag = TRUE
min.isect.size=0
correction.method = "fdr"

for (path.label in c("KEGG", "REACT")) {
  tryCatch({
    # gprofiler(query = query.names, organism="rnorvegicus", src_filter = path.label, )
  }, error=function(e) {
    print(e)
  }
  )
}

for (go.label in c("BP", "CC", "MF")) {
  # tryCatch({
    pp <- gprofiler(query = query.names, organism="rnorvegicus", significant = significant.flag, min_isect_size = min.isect.size, correction_method = correction.method, src_filter = go.label, exclude_iea = exclude.electr.flag)
    if(length(attr(pp, "networks")$BIOGRID)!=0) {
      print(go.label)
    }else{
      print("0")
    }
  # }, error=function(e) {
  #   print(e)
  # }
  # )
}




