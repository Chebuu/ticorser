# 
# require("gdata")
# 
# qpcr.data <- read.xls(xls = "heatmap_files/qPCR_modified.XLSX", stringsAsFactors=FALSE, header=TRUE)
# 
# # qpcr.data.cerv.sham <- qpcr.data[, c(1:6)]
# 
# # colnames(qpcr.data.cerv.sham) <- c("GeneName", "EnsemblID", "C3d", "C1w", "C2w", "C8w")
# 
# # qpcr.data.cerv.sham <- qpcr.data.cerv.sham[-1,]
# 
# 
# 
# head(qpcr.data)
# 
# colnames(qpcr.data) <- c("GeneName", "EnsemblID", "C03d", "C07d", "C14d", "C56d", "T03d", "T07d", "T14d", "T56d", "NCT03d", "NCT07d", "NCT14d", "NCT56d", "logNCT03d", "logNCT07d", "logNCT14d", "logNCT56d")
# 
# qpcr.data <- qpcr.data[order(qpcr.data$EnsemblID),]
# qpcr.ens.gs.map <- qpcr.data[,c("EnsemblID", "GeneName")]
# 
# rownames(qpcr.data) <- qpcr.data$GeneName
# qpcr.data.rn5 <- qpcr.data[-which(rownames(qpcr.data) %in% c("Acta2", "LIF")),]
# qpcr.data.log2.rn5 <- log(qpcr.data.rn5[, c(3:14)] , base = 2)
# 
# qpcr.data.cerv.sham <- qpcr.data.log2.rn5[,c(1:4)]
# qpcr.data.thor.sham <- qpcr.data.log2.rn5[,c(5:8)]
# qpcr.data.cerv.thor <- qpcr.data.log2.rn5[,c(9:12)]
# 
# #### rna counts
# nn.counts <- normalized.counts.prop.uqua
# 
# nn.counts <- nn.counts[order(rownames(nn.counts)),]
# 
# ind.genes <- which(rownames(nn.counts) %in% qpcr.ens.gs.map$EnsemblID)
# nn.counts.genes <- nn.counts[ind.genes,]
# 
# nn.counts.genes$genes <- rownames(qpcr.data.rn5)
# rownames(nn.counts.genes) <- nn.counts.genes$genes
# nn.counts.genes$genes
# 
# 
# cerv.zscores <- ComputeLogFCZscores(nn.counts.genes, cervical.design.matrix)
# thor.zscores <- ComputeLogFCZscores(nn.counts.genes, thoracic.design.matrix)
# cerv.thor.zscores <- ComputeLogFCZscores(nn.counts.genes, cerv.thor.design.matrix)
# 
# qpcr.data.cerv.sham.z <- t(ComputeZscore(qpcr.data.cerv.sham, byRowCol = 1))
# qpcr.data.thor.sham.z <- t(ComputeZscore(qpcr.data.thor.sham, byRowCol = 1))
# qpcr.data.cerv.thor.z <- t(ComputeZscore(qpcr.data.cerv.thor, byRowCol = 1))
# 
# cervical.log.fc.z <- cbind(qpcr.data.cerv.sham.z, cerv.zscores$lfcz)
# thoracic.log.fc.z <- cbind(qpcr.data.thor.sham.z, thor.zscores$lfcz)
# cerv.thor.log.fc.z <- cbind(qpcr.data.cerv.thor.z, cerv.thor.zscores$lfcz)
# 
# 
# PlotHeatmap(fold.changes = cervical.log.fc.z, row.labels.percent = 0.9, scale = "none", title = "cervical pcr-rna")
# 
# PlotHeatmap(fold.changes = thoracic.log.fc.z, row.labels.percent = 0.9, scale = "none", title = "thoracic pcr-rna")
# PlotHeatmap(fold.changes = cerv.thor.log.fc.z, row.labels.percent = 0.9, scale = "none", title = "cerv-thor pcr-rna")
# 
# # 
# # 
# # PlotQpcrRnaHeatmap <- function(qpcr, rna, threshold=Inf, heatmap.title="heatmap", scale.heat="none") {
# #   qpcr[qpcr>threshold] <- threshold
# #   qpcr[qpcr<(-threshold)] <- (-threshold)
# #   
# #   rna[rna>threshold]<-threshold
# #   rna[rna<(-threshold)]<- (-threshold)
# #   
# #   heat.bind<-as.matrix(cbind(qpcr,rna))
# #   PlotHeatmap(fold.changes = heat.bind, row.labels.percent = 0.9, scale = scale.heat, title = heatmap.title)
# #   
# # }
# # 
# # ComputeCorrelation <- function(qpcr, rna, threshold=Inf){
# #   qpcr[qpcr>threshold] <- threshold
# #   qpcr[qpcr<(-threshold)] <- (-threshold)
# #   
# #   rna[rna>threshold]<-threshold
# #   rna[rna<(-threshold)]<- (-threshold)
# #   
# #   corr <- list()
# #   for(i in 1:dim(qpcr)[1]) {
# #     # print(as.numeric(binded.values[i,]))
# #     # print(as.numeric(binded.values[i,]))
# #     corr[[i]]<-cor(as.numeric(qpcr[i,]), as.numeric(rna[i,]))
# #     # print(corr[[i]])
# #     # print("---")
# #   }
# #   names(corr) <- rownames(qpcr)
# #   
# #   corr.unl <- unlist(corr)
# #   cat("---\nNAs:\n ")
# #   print( corr.unl[is.na(corr.unl)] )
# #   cat("---\nanti-correlated:\n")
# #   print( corr.unl[which(corr.unl < 0)] )
# #   cat("---\nSummary:\n")
# #   print(summary(corr.unl))
# #   boxplot(corr.unl)
# #   
# #   return(corr.unl)
# # }
# 
# 
# qpcr.log <- qpcr.data.cerv.sham
# rna.log <- cerv.zscores$lfc
# thr <- Inf
# 
# PlotQpcrRnaHeatmap(qpcr.log, rna.log, heatmap.title = "cerv qpcr-rna logfc", threshold = thr)
# logfc.corr <- ComputeCorrelation(qpcr.log, rna.log, threshold = thr)
# 
# qpcr.logz <- qpcr.data.cerv.sham.z
# rna.logz <- cerv.zscores$lfcz
# PlotQpcrRnaHeatmap(qpcr.logz, rna.logz, heatmap.title = "cerv qpcr-rna Z")
# logfcz.corr <- ComputeCorrelation(qpcr.logz, rna.logz)
# 
# ################################
# 
# cerv.03.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/03d/Cervical_deseq/DE_results/DeSeq/Cervical_03d_prop_uqua_DeSeq_all.tsv")
# cerv.07.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/07d/Cervical_deseq/DE_results/DeSeq/Cervical_07d_prop_uqua_DeSeq_all.tsv")
# cerv.14.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/14d/Cervical_deseq/DE_results/DeSeq/Cervical_14d_prop_uqua_DeSeq_all.tsv")
# cerv.56.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/56d/Cervical_deseq/DE_results/DeSeq/Cervical_56d_prop_uqua_DeSeq_all.tsv")
# # 
# # OrderByRownames <- function(named.df) {
# #   return(named.df[order(rownames(named.df)),])
# # }
# # 
# 
# cerv.03.deseq <- OrderByRownames(cerv.03.deseq)
# cerv.07.deseq <- OrderByRownames(cerv.07.deseq)
# cerv.14.deseq <- OrderByRownames(cerv.14.deseq)
# cerv.56.deseq <- OrderByRownames(cerv.56.deseq)
# 
# cerv.deseq.binded.logs <- cbind(cerv.03.deseq$log2FoldChange, cerv.07.deseq$log2FoldChange, cerv.14.deseq$log2FoldChange, cerv.56.deseq$log2FoldChange)
# cerv.deseq.binded.padj <- cbind(cerv.03.deseq$padj, cerv.07.deseq$padj, cerv.14.deseq$padj, cerv.56.deseq$padj)
# 
# 
# 
# 
# rownames(cerv.deseq.binded.logs) <- rownames(cerv.03.deseq)
# rownames(cerv.deseq.binded.padj) <- rownames(cerv.03.deseq)
# 
# colnames(cerv.deseq.binded.logs) <- c("03d", "07d", "14d", "56d")
# colnames(cerv.deseq.binded.padj) <- c("03d", "07d", "14d", "56d")
# 
# 
# cerv.deseq.binded.logs.genes <- cerv.deseq.binded.logs[which(rownames(cerv.deseq.binded.logs) %in% qpcr.ens.gs.map$EnsemblID),]
# 
# 
# qpcr.ens.gs.map <- qpcr.ens.gs.map[order(qpcr.ens.gs.map$EnsemblID),]
# 
# cerv.deseq.binded.logs.genes.c <- cbind(cerv.deseq.binded.logs.genes, qpcr.ens.gs.map[which(qpcr.ens.gs.map$EnsemblID %in% rownames(cerv.deseq.binded.logs.genes)),])
# 
# rownames(cerv.deseq.binded.logs.genes.c) <- cerv.deseq.binded.logs.genes.c$GeneName
# cerv.deseq.binded.logs.genes.c <- cerv.deseq.binded.logs.genes.c[,-c(5,6)]
# 
# qpcr.data.cerv.sham.o <- OrderByRownames(qpcr.data.cerv.sham)
# cerv.deseq.binded.logs.genes.c.o <- OrderByRownames(cerv.deseq.binded.logs.genes.c)
# 
# PlotQpcrRnaHeatmap(qpcr = qpcr.data.cerv.sham.o, rna = cerv.deseq.binded.logs.genes.c.o, heatmap.title = "PCR-DESeq2 logfc")
# 
# qpcr.data.cerv.sham.o.z <- t(ComputeZscore(qpcr.data.cerv.sham.o))
# cerv.deseq.binded.logs.genes.c.o.z <- t(ComputeZscore(cerv.deseq.binded.logs.genes.c.o))
# 
# PlotQpcrRnaHeatmap(qpcr = qpcr.data.cerv.sham.o.z, rna = cerv.deseq.binded.logs.genes.c.o.z, heatmap.title = "PCR-DESeq2 row z")
# 
# ###########
# 
# qpcr.max <- apply(X = qpcr.data.cerv.sham.o, MARGIN = 1, function(x){max(abs(x))})
# 
# max.cost.interval <- 100
# qpcr.100 <- (qpcr.data.cerv.sham.o / qpcr.max) * max.cost.interval
# 
# rna.max <- apply(X = cerv.deseq.binded.logs.genes.c.o, MARGIN = 1, function(x){max(abs(x))})
# deseq.100 <- (cerv.deseq.binded.logs.genes.c.o / rna.max ) *max.cost.interval
# 
# PlotQpcrRnaHeatmap(qpcr = qpcr.100, rna = deseq.100, heatmap.title = "PCR-DESeq2 lfc max 100")
# 
# ComputeCorrelation(qpcr = qpcr.data.cerv.sham.o, rna = cerv.deseq.binded.logs.genes.c.o)
# 
# ####  just fc
#  # cerv.qpcr.fca <- qpcr.data.rn5[, c(3:6)]
# cerv.qpcr.fc<- 2^qpcr.data.cerv.sham.o
# cerv.qpcr.max <- apply(X = cerv.qpcr.fc, MARGIN = 1, function(x){max(x)})
# cerv.qpcr.100 <- (cerv.qpcr.fc / cerv.qpcr.max) * max.cost.interval
# 
# 
# 
# cerv.deseq.binded.fc <- 2^cerv.deseq.binded.logs.genes.c.o
# cerv.deseq.max <- apply(X = cerv.deseq.binded.fc, MARGIN = 1, function(x){max(x)})
# cerv.deseq.100 <- (cerv.deseq.binded.fc / cerv.deseq.max ) *max.cost.interval
# 
# PlotQpcrRnaHeatmap(qpcr = cerv.qpcr.100, rna = cerv.deseq.100, heatmap.title = "PCR-DESeq2 fc max 100")
# 
# 
# ComputeCorrelation(qpcr = cerv.qpcr.fc, rna = cerv.deseq.binded.fc)
# dim(cerv.qpcr.100)
# 
# 
# 



####################### pcr

require("gdata")

qpcr.data <- read.xls(xls = "/media/dario/dati/time_course/heatmap_files/qPCR_modified.XLSX", stringsAsFactors=FALSE, header=TRUE)

head(qpcr.data)
colnames(qpcr.data) <- c("GeneName", "EnsemblID", "C03d", "C07d", "C14d", "C56d", "T03d", "T07d", "T14d", "T56d", "NCT03d", "NCT07d", "NCT14d", "NCT56d", "logNCT03d", "logNCT07d", "logNCT14d", "logNCT56d")
qpcr.data <- qpcr.data[order(qpcr.data$EnsemblID),]
qpcr.ens.gs.map <- qpcr.data[,c("EnsemblID", "GeneName")]
qpcr.ens.gs.map <- qpcr.ens.gs.map[order(qpcr.ens.gs.map$EnsemblID),]
rownames(qpcr.data) <- qpcr.data$GeneName
qpcr.data.rn5 <- qpcr.data[-which(rownames(qpcr.data) %in% c("Acta2", "LIF")),]

cerv.qpcr.fc <- qpcr.data.rn5[,c(3:6)]
save("cerv.qpcr.fc", file="heatmap_files/report_files/cervical_qpcr_fc.RData")
thor.qpcr.fc <- qpcr.data.rn5[,c(7:10)]
cerv.thor.qpcr.fc <- qpcr.data.rn5[,c(11:14)]

############## deseq

cerv.03.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/03d/Cervical_deseq/DE_results/DeSeq/Cervical_03d_prop_uqua_DeSeq_all.tsv")
cerv.07.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/07d/Cervical_deseq/DE_results/DeSeq/Cervical_07d_prop_uqua_DeSeq_all.tsv")
cerv.14.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/14d/Cervical_deseq/DE_results/DeSeq/Cervical_14d_prop_uqua_DeSeq_all.tsv")
cerv.56.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/56d/Cervical_deseq/DE_results/DeSeq/Cervical_56d_prop_uqua_DeSeq_all.tsv")

cerv.03.deseq <- OrderByRownames(cerv.03.deseq)
cerv.07.deseq <- OrderByRownames(cerv.07.deseq)
cerv.14.deseq <- OrderByRownames(cerv.14.deseq)
cerv.56.deseq <- OrderByRownames(cerv.56.deseq)

cerv.deseq.binded.logs <- cbind(cerv.03.deseq$log2FoldChange, cerv.07.deseq$log2FoldChange, cerv.14.deseq$log2FoldChange, cerv.56.deseq$log2FoldChange)
rownames(cerv.deseq.binded.logs) <- rownames(cerv.03.deseq)
colnames(cerv.deseq.binded.logs) <- c("03d", "07d", "14d", "56d")

cerv.deseq.binded.logs.genes <- cerv.deseq.binded.logs[which(rownames(cerv.deseq.binded.logs) %in% qpcr.ens.gs.map$EnsemblID),]
cerv.deseq.binded.logs.genes.c <- cbind(cerv.deseq.binded.logs.genes, qpcr.ens.gs.map[which(qpcr.ens.gs.map$EnsemblID %in% rownames(cerv.deseq.binded.logs.genes)),])
rownames(cerv.deseq.binded.logs.genes.c) <- cerv.deseq.binded.logs.genes.c$GeneName
cerv.deseq.binded.logs.genes.c <- cerv.deseq.binded.logs.genes.c[,-c(5,6)]

cerv.deseq.fc <- 2^cerv.deseq.binded.logs.genes.c
save("cerv.deseq.fc", file = "heatmap_files/report_files/cerv_deseq_fc.Rdata")

#### cervical fc heatmap
max.cost.interval=100
cerv.qpcr.max <- apply(X = cerv.qpcr.fc, MARGIN = 1, function(x){max(x)})
cerv.qpcr.100 <- (cerv.qpcr.fc / cerv.qpcr.max) * max.cost.interval

cerv.deseq.max <- apply(X = cerv.deseq.fc, MARGIN = 1, function(x){max(x)})
cerv.deseq.100 <- (cerv.deseq.fc / cerv.deseq.max ) *max.cost.interval

PlotQpcrRnaHeatmap(qpcr = cerv.qpcr.100, rna = cerv.deseq.100, heatmap.title = "PCR-DESeq2 fc max 100", color.palette = c("yellow", "red"))

ComputeCorrelation(qpcr = cerv.qpcr.100, rna = cerv.deseq.100)

#### cervical logfc heatmap
cerv.qpcr.lfc <- log2(cerv.qpcr.fc)
cerv.deseq.lfc <- log2(cerv.deseq.fc)

PlotQpcrRnaHeatmap(qpcr = cerv.qpcr.lfc, rna = cerv.deseq.lfc, heatmap.title = "PCR-DESeq2 lfc")
ComputeCorrelation(qpcr = cerv.qpcr.lfc, rna = cerv.deseq.lfc)

###########

no.genes <- c("Slc1a", "aifm1", "Hipk3", "S1PR3", "PDGFRB", "Stat3", "Plaur")
cerv.nogenes.qpcr.lfc <- cerv.qpcr.lfc[-which(rownames(cerv.qpcr.lfc) %in% no.genes),]
cerv.nogenes.deseq.lfc <- cerv.deseq.lfc[-which(rownames(cerv.deseq.lfc) %in% no.genes),]
PlotQpcrRnaHeatmap(qpcr = cerv.nogenes.qpcr.lfc, rna = cerv.nogenes.deseq.lfc, heatmap.title = "PCR-DESeq2 lfc")


#### cervical logfc scale row heatmap
PlotQpcrRnaHeatmap(qpcr = cerv.qpcr.lfc, rna = cerv.deseq.lfc, heatmap.title = "PCR-DESeq2 scale-row", scale.heat = "row")

#### cervical logfc row zscore heatmap
cerv.qpcr.lfc.z <- t(ComputeZscore(cerv.qpcr.lfc))
cerv.deseq.lfc.z <- t(ComputeZscore(cerv.deseq.lfc))
PlotQpcrRnaHeatmap(qpcr = cerv.qpcr.lfc.z, rna = cerv.deseq.lfc.z, heatmap.title = "PCR-DESeq2 lfc row z-score")
ComputeCorrelation(qpcr = cerv.qpcr.100, rna = cerv.deseq.100)


summary.genes <- as.data.frame(ddf[which(ddf$genes %in% rownames(cerv.qpcr.lfc)),], stringsAsFactors =FALSE)
summary.genes$values <- as.character(summary.genes$values)

summary.genes[which(summary.genes$values=="both"), "values"] <- "BOTH"

pp <- ggplot(summary.genes, aes(x=times, y=genes, fill=values)) + geom_dotplot(binaxis='y', stackdir='center',  dotsize=0.25) + scale_fill_manual(values = c("DOWN"="red", "0"="white", "UP"="green", "LRT"="turquoise", "Wald"="plum1", "both"="orange"))
pp



############### thoracic heatmap

thor.03.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/03d/Thoracic_deseq/DE_results/DeSeq/Thoracic_03d_prop_uqua_DeSeq_all.tsv")
thor.07.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/07d/Thoracic_deseq/DE_results/DeSeq/Thoracic_07d_prop_uqua_DeSeq_all.tsv")
thor.14.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/14d/Thoracic_deseq/DE_results/DeSeq/Thoracic_14d_prop_uqua_DeSeq_all.tsv")
thor.56.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/56d/Thoracic_deseq/DE_results/DeSeq/Thoracic_56d_prop_uqua_DeSeq_all.tsv")

thor.03.deseq <- OrderByRownames(thor.03.deseq)
thor.07.deseq <- OrderByRownames(thor.07.deseq)
thor.14.deseq <- OrderByRownames(thor.14.deseq)
thor.56.deseq <- OrderByRownames(thor.56.deseq)

thor.deseq.binded.logs <- cbind(thor.03.deseq$log2FoldChange, thor.07.deseq$log2FoldChange, thor.14.deseq$log2FoldChange, thor.56.deseq$log2FoldChange)
rownames(thor.deseq.binded.logs) <- rownames(thor.03.deseq)
colnames(thor.deseq.binded.logs) <- c("03d", "07d", "14d", "56d")

thor.deseq.binded.logs.genes <- thor.deseq.binded.logs[which(rownames(thor.deseq.binded.logs) %in% qpcr.ens.gs.map$EnsemblID),]
thor.deseq.binded.logs.genes.c <- cbind(thor.deseq.binded.logs.genes, qpcr.ens.gs.map[which(qpcr.ens.gs.map$EnsemblID %in% rownames(thor.deseq.binded.logs.genes)),])
rownames(thor.deseq.binded.logs.genes.c) <- thor.deseq.binded.logs.genes.c$GeneName
thor.deseq.binded.logs.genes.c <- thor.deseq.binded.logs.genes.c[,-c(5,6)]

thor.deseq.fc <- 2^thor.deseq.binded.logs.genes.c
thor.qpcr.fc


## logfc heatmap
thor.qpcr.lfc <- log2(thor.qpcr.fc)
thor.deseq.lfc <- log2(thor.deseq.fc)

PlotQpcrRnaHeatmap(qpcr = thor.qpcr.lfc, rna = thor.deseq.lfc, heatmap.title = "thoracic PCR-DESeq2 lfc")
ComputeCorrelation(qpcr = thor.qpcr.lfc, rna = thor.deseq.lfc)



