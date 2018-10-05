source("https://bioconductor.org/biocLite.R")
biocLite("RUVSeq")
biocLite("zebrafishRNASeq")
biocLite("mixOmics")
library("RUVSeq")
library("zebrafishRNASeq")
library("mixOmics")

data("zfGenes")
head(zfGenes)
tail(zfGenes)
grep("^ERCC", rownames(zfGenes))

filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]

tail(filtered)
x <- as.factor(rep(c("Ctl", "Trt"), each=3))

mtx.set <- as.matrix(filtered)
pheno.set <- data.frame(x, row.names=colnames(filtered))

set <- newSeqExpressionSet(mtx.set, phenoData=pheno.set)
set 

thoracic.dataset <- thoracic.filtered.counts.wilcoxon

emp.genes <- thoracic.de.genes[is.na(thoracic.de.genes$padj),]
dim(emp.genes)

emp.counts.genes <- thoracic.filtered.counts.wilcoxon[ which(rownames(thoracic.filtered.counts.wilcoxon) %in% rownames(emp.genes)),]

# h<-hist(as.matrix(emp.counts.genes), breaks = 100)
# 
# dds<-apply(as.matrix(emp.counts.genes), 1, density)
# 
# density(as.numeric(emp.counts.genes[3,]))
# 
# str(dds)
# 
# bws <- lapply(dds, function(x) { x$bw} )
# bwsunl <- unlist(bws)
# hist(bwsunl, breaks = 10)

par(mfrow = c(53, 50))
for(i in 1:dim(emp.counts.genes)[1] ) {
  hist(emp.counts.genes[i,], breaks = 5, main=rownames(emp.counts.genes)[i])
}


emp.ctl.genes <- rownames(tail(emp.genes[order(emp.genes$pvalue>0),], n=75))

# emp.ctl.genes <- hk.est

# thoracic.design.matrix.reo <- thoracic.design.matrix[order(thoracic.design.matrix$Times),]
# thoracic.dataset.reo <- thoracic.dataset[, order(colnames(thoracic.dataset))]
# colnames(thoracic.dataset.reo)

thoracic.dataset.set <- newSeqExpressionSet( as.matrix(thoracic.dataset),
                                             phenoData= data.frame(
                                               thoracic.design.matrix$Conditions, row.names = rownames(thoracic.design.matrix)
                                             ))


ruved.thoracic.set <- RUVg(thoracic.dataset.set, emp.ctl.genes, k=1)


# ruved.thoracic.set.03 <- RUVg(thoracic.dataset.set[,1:8], emp.ctl.genes, k=1)
# ruved.thoracic.set.07 <- RUVg(thoracic.dataset.set[,9:16], emp.ctl.genes, k=1)
# ruved.thoracic.set.14 <- RUVg(thoracic.dataset.set[,17:24], emp.ctl.genes, k=1)
# ruved.thoracic.set.56 <- RUVg(thoracic.dataset.set[,25:30], emp.ctl.genes, k=1)
plotRLE(thoracic.dataset.set, outline=FALSE, ylim=c(-4,4))
plotRLE(ruved.thoracic.set, outline=FALSE, ylim=c(-4,4))
plotPCA(ruved.thoracic.set, cex=0.8)
PlotPCAFunction(ruved.thoracic.set@assayData$counts, design.dataframe = thoracic.design.matrix, scale = TRUE)
PlotPCAFunction(ruved.thoracic.set@assayData$normalizedCounts, design.dataframe = thoracic.design.matrix, scale = TRUE)

#####################
differences <- makeGroups(thoracic.design.matrix$Conditions)
set3 <- RUVs(thoracic.dataset.set, rownames(thoracic.dataset.set), k=1, differences)

# plotRLE(thoracic.dataset.set, outline=FALSE, ylim=c(-4,4))
plotRLE(set3, outline=FALSE, ylim=c(-4,4))
plotPCA(set3, cex=0.8)
PlotPCAFunction(set3@assayData$normalizedCounts, design.dataframe = thoracic.design.matrix, scale = FALSE)

# PlotPCAFunction(thoracic.filtered.counts.wilcoxon, design.dataframe = thoracic.design.matrix)
# PlotPCAFunction(as.data.frame(ruved.thoracic.set@assayData$counts), design.dataframe = thoracic.design.matrix)
PlotPCAFunction(ruved.thoracic.set.03@assayData$normalizedCounts, design.dataframe = thoracic.design.matrix)
PlotPCAFunction(thoracic.filtered.counts.wilcoxon.normalized, design.dataframe = thoracic.design.matrix)


#########################################################################
# 
# thoracic.filtered.counts.wilcoxon.normalized <- NormalizeData(thoracic.filtered.counts.wilcoxon, norm.type = "fqua")
# PlotPCAFunction(thoracic.filtered.counts.wilcoxon.normalized, design.dataframe = thoracic.design.matrix)
# boxplot(log(thoracic.filtered.counts.wilcoxon.normalized+1))
# 
# thoracic.filtered.counts.wilcoxon.normalized.03 <- NormalizeData(thoracic.filtered.counts.wilcoxon[,1:5], norm.type = "fqua")
# boxplot(thoracic.filtered.counts.wilcoxon.normalized.03+1)
# PlotPCAFunction(thoracic.filtered.counts.wilcoxon.normalized.03, design.dataframe = thoracic.design.matrix)
# 
# thoracic.filtered.counts.wilcoxon.just.de <- thoracic.filtered.counts.wilcoxon[-which(rownames(thoracic.filtered.counts.wilcoxon) %in% emp.ctl.genes), ]
# thoracic.filtered.counts.wilcoxon.normalized.justde.07 <- NormalizeData(thoracic.filtered.counts.wilcoxon.just.de[,9:13], norm.type = "fqua")
# PlotPCAFunction(log(thoracic.filtered.counts.wilcoxon.normalized.justde.07+1), design.dataframe = thoracic.design.matrix[9:13,], scale=FALSE)
# boxplot(log(thoracic.filtered.counts.wilcoxon.normalized.justde.07+1))
# 
# df<-thoracic.filtered.counts.wilcoxon.normalized.justde.07
# names(df[, apply(df, 1, var, na.rm=FALSE) ==  0])








