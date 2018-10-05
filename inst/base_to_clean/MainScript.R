## MAIN SCRIPT
  source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/TimeCourseFunctions.R")
  source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/MainFunction.R")


setwd(dir = "/media/dario/dati/time_course_4/")

cervical.files <- list.files("Cervical/cervical_bam/")
cervical.bam <- cervical.files[grep(pattern = "*.bam", cervical.files)]

thoracic.files <- list.files("Thoracic/thoracic_bam/")
thoracic.bam <- thoracic.files[grep(pattern = "*.bam", thoracic.files)]
thoracic.bam <- thoracic.bam[-grep(pattern = "*.bai", thoracic.bam)]

rat5.ensembl.downloaded.gtf <- "/media/dario/dati/time_course_4/downloaded_references/ensembl_site/Rattus_norvegicus.Rnor_5.0.79_gene_biotype_protein_coding_header.gtf"
all.bam.folder <- "/media/dario/dati/time_course_4/cervical_thoracic/all_bams/"
cervical.thoracic.whole.bams <- list.files(all.bam.folder, pattern="bam$",full.names=TRUE)
cervical.thoracic.whole.counts.ens.gtf <- CountBamFilesFeatureCounts(bam.files = cervical.thoracic.whole.bams, annotation.file = rat5.ensembl.downloaded.gtf, output.folder = output.results.folder, prefix.output.file="whole_cervical_thoracic_18")


cervical.bam.folder <- "/media/dario/dati/time_course/Cervical/cervical_bam"
# cervical.annotation.file <- "/media/dario/dati/time_course/Cervical/cuffmerge/merged_assembly/merged.gtf"
cervical.design.file <- file.path("/media/dario/dati/time_course/Cervical/design_file/cervical_design_file.txt")

cervical.bam.files <- list.files(cervical.bam.folder, pattern="bam$",full.names=TRUE)

cervical.output.counts.folder <- file.path("/media/dario/dati/time_course/Cervical", "results", "counts")
cervical.output.de.folder <- file.path("/media/dario/dati/time_course/Cervical", "results", "de_genes", "DESeq2_TimeCourse")

# cervical.counts <- CountBamFilesFeatureCounts(bam.files = cervical.bam.files, annotation.file = cervical.annotation.file, counts.folder = cervical.output.counts.folder)
# colnames(cervical.counts)

#cervical.counts.ens.gtf <- CountBamFilesFeatureCounts(bam.files = cervical.bam.files, annotation.file = rat5.ensembl.downloaded.gtf, counts.folder = cervical.output.counts.folder, gtf.attr.type = "gene_id", prefix.output.file="Cervical")
dim(cervical.counts.ens.gtf)
# cervical.counts <- read.table(file.path(cervical.output.counts.folder, "Counts_FeatureCounts.txt"))

## automatizzare processo di rinomina e ordinamento delle colonne durante la fase di conteggio
cervical.renamed.counts <- cervical.counts.ens.gtf
colnames(cervical.renamed.counts) <- c( "cervical_14d_t1", "cervical_14d_t2", "cervical_14d_t3", "cervical_14d_t4", "cervical_14d_t5",
                                        "cervical_14d_u1", "cervical_14d_u2", "cervical_14d_u3",
                                        "cervical_03d_t1", "cervical_03d_t2", "cervical_03d_t3", "cervical_03d_t4", "cervical_03d_t5",
                                        "cervical_03d_u1", "cervical_03d_u2", "cervical_03d_u3",
                                        "cervical_56d_t1", "cervical_56d_t2", "cervical_56d_t3", "cervical_56d_t4", "cervical_56d_t5",
                                        "cervical_56d_u1", "cervical_56d_u2", "cervical_56d_u3",
                                        "cervical_07d_t1", "cervical_07d_t2", "cervical_07d_t3", "cervical_07d_t4", "cervical_07d_t5",
                                        "cervical_07d_u1", "cervical_07d_u2", "cervical_07d_u3"
)

# cervical.renamed.counts.reordered <- cervical.renamed.counts[ , c( 9:16, 25:32, 1:8, 17:24 )]
cervical.renamed.counts.reordered <- cervical.renamed.counts[ , order(colnames(cervical.renamed.counts))]
colnames(cervical.renamed.counts.reordered)
WriteDataFrameAsTsv(cervical.renamed.counts.reordered, file.name.path = file.path(cervical.output.counts.folder, "cervical_ensgtf_reordered_counts"))
# save(cervical.renamed.counts.reordered, file = "cervical_counts.RData")

load(file = "cervical_counts.RData")
cervical.output.folder <- file.path("/media/dario/dati/time_course/Cervical", "results")

cervical.design.matrix <- read.table(cervical.design.file, header = TRUE, row.names = 1)
cervical.de.genes.wilcoxon.uqua <- MainFunction(counts.dataframe = cervical.renamed.counts.reordered, design.matrix = cervical.design.matrix, output.results.folder = cervical.output.folder, prefix.output.file = "cervical", prefix.plot.label = "cervical", is.time.course = TRUE)
save(cervical.de.genes.wilcoxon.uqua, file="cervical.de.genes.wilcoxon.uqua.RData")
load(file="cervical.de.genes.wilcoxon.uqua.RData")
cervical.wilcoxon.unnormalized <- ReadDataFrameFromTsv("/media/dario/dati/time_course/Cervical/results/counts/FeatureCounts/wilcoxon/cervical__filtered_counts__Wilcoxon.tsv")
cervical.wilcoxon.uqua <- ReadDataFrameFromTsv("/media/dario/dati/time_course/Cervical/results/counts/FeatureCounts/cervical__filtered_counts__Wilcoxon__normalized__uqua.tsv")


cervical.output.plots.folder <- file.path(cervical.output.folder, "plots")
# PlotPCAPlotlyFunction(counts.data.frame = cervical.wilcoxon.unnormalized, design.matrix = cervical.design.matrix, colour.design.column.str = "Times", shape.design.column.str = "Conditions", output.path = cervical.output.plots.folder, prefix.plot = "cervical_unnormalized", title="Unnormalized PCA")
cervical.wilcoxon.estimated.neg.genes <- EstimateNegativeControlGenesForRUV(de.genes = cervical.de.genes.wilcoxon.uqua, n.tail.genes = 2000, counts.dataset = cervical.wilcoxon.unnormalized, n.genes.per.hist.break = 1)
cervical.de.genes.wilcoxon.ruvg <- MainFunction(counts.dataframe = cervical.renamed.counts.reordered, design.matrix = cervical.design.matrix, output.results.folder = cervical.output.folder, prefix.output.file = "cervical", prefix.plot.label = "cervical", is.time.course = TRUE, normalization.method = "ruvg", estimated.genes = cervical.wilcoxon.estimated.neg.genes)
cervical.wilcoxon.ruvg <- ReadDataFrameFromTsv("/media/dario/dati/time_course/Cervical/results/counts/FeatureCounts/wilcoxon/ruvg/cervical__filtered_counts__Wilcoxon__normalized__ruvg.tsv")
counts.data.frame = cervical.wilcoxon.ruvg; design.matrix = cervical.design.matrix; colour.design.column.str = "Times"; shape.design.column.str = "Conditions";
output.path = cervical.output.plots.folder; prefix.plot = "cervical_wilcoxon_ruvg"; title="Cervical wilcoxon Normalized RUVg PCA"



cervical.de.genes.proportion.uqua <- MainFunction(counts.dataframe = cervical.renamed.counts.reordered, design.matrix = cervical.design.matrix, output.results.folder = cervical.output.folder, prefix.output.file = "cervical", prefix.plot.label = "cervical", is.time.course = TRUE, filter.method = "Proportion")
cervical.proportion.unnormalized <- ReadDataFrameFromTsv("/media/dario/dati/time_course/Cervical/results/counts/FeatureCounts/proportion/cervical__filtered_counts__Proportion.tsv")
cervical.proportion.uqua <- ReadDataFrameFromTsv("/media/dario/dati/time_course/Cervical/results/counts/FeatureCounts/proportion/uqua/cervical__filtered_counts__Proportion__normalized__uqua.tsv")

cervical.proportion.estimated.neg.genes <- EstimateNegativeControlGenesForRUV(de.genes = cervical.de.genes.proportion.uqua, n.tail.genes = 2000, counts.dataset = cervical.proportion.unnormalized, n.genes.per.hist.break = 1)
cervical.de.genes.proportion.ruvg <- MainFunction(counts.dataframe = cervical.renamed.counts.reordered, design.matrix = cervical.design.matrix, output.results.folder = cervical.output.folder, prefix.output.file = "cervical", prefix.plot.label = "cervical", is.time.course = TRUE, normalization.method = "ruvg", estimated.genes = cervical.proportion.estimated.neg.genes, filter.method = "Proportion")

cervical.proportion.ruvg <- ReadDataFrameFromTsv("/media/dario/dati/time_course/Cervical/results/counts/FeatureCounts/proportion/ruvg/cervical__filtered_counts__Proportion__normalized__ruvg.tsv")
counts.data.frame = cervical.proportion.ruvg; design.matrix = cervical.design.matrix; colour.design.column.str = "Times"; shape.design.column.str = "Conditions";
output.path = cervical.output.plots.folder; prefix.plot = "cervical_proportion_ruvg"; title="Cervical Normalized RUVg PCA"

cervical.de.genes.proportion.uqua.sign <- SignificantDeGenesPAdj(cervical.de.genes.proportion.uqua)
WriteDataFrameAsTsv(cervical.de.genes.proportion.uqua.sign, file.name.path = "/media/dario/dati/time_course/Cervical/Cervical_results/de_genes/proportion/uqua/cervical__filtered_counts__Proportion__normalized__uqua_de_genes_LRT_significative_005")

cervical.de.genes.proportion.ruvg.sign <- SignificantDeGenesPAdj(cervical.de.genes.proportion.ruvg)
WriteDataFrameAsTsv(cervical.de.genes.proportion.ruvg.sign, file.name.path = "/media/dario/dati/time_course/Cervical/Cervical_results/de_genes/proportion/ruvg/cervical__filtered_counts__Proportion__normalized__ruvg_de_genes_LRT_significative_005")

# counts.dataframe = cervical.renamed.counts.reordered; design.matrix = cervical.design.matrix; output.results.folder = cervical.output.folder; prefix.output.file = "cervical"; prefix.plot.label = "cervical"; is.time.course = TRUE
# cervical.design.matrix <- read.table(cervical.design.file, header = TRUE, row.names = 1)
# cervical.filtered.counts.wilcoxon <- FilterLowCounts(counts.dataframe = cervical.renamed.counts.reordered, design.dataframe = cervical.design.matrix, is.normalized = FALSE, method.type = "Wilcoxon", cpm.cutoff = 1, cv.percentage = 1)
# WriteDataFrameAsTsv(cervical.filtered.counts.wilcoxon, file.name.path = file.path(cervical.output.counts.folder, "cervical_ensgtf_filtered_counts_wilcoxon_17555"))
# # filtered.counts.wilcoxon.reordered <- filtered.counts.wilcoxon[ , c( 9:16, 25:32, 1:8, 17:24 )]
#
#
# boxplot(log(cervical.renamed.counts.reordered+1))
# cervical.normalized.filtered.counts.wilcoxon <- NormalizeData(cervical.filtered.counts.wilcoxon, norm.type = "uqua")
# WriteDataFrameAsTsv(cervical.filtered.counts.wilcoxon, file.name.path = file.path(cervical.output.counts.folder, "cervical_ens_gtf_filtered_counts_wilcoxon_17371_normalized_upperquartile"))
# boxplot(log(cervical.normalized.filtered.counts.wilcoxon+1))
# colnames(cervical.normalized.filtered.counts.wilcoxon)
# cervical.normalized.filtered.counts.wilcoxon.reordered <- cervical.normalized.filtered.counts.wilcoxon[ , c( 9:16, 25:32, 1:8, 17:24 )]
# colnames(cervical.normalized.filtered.counts.wilcoxon.reordered)

# write.table(x = cervical.normalized.filtered.counts.wilcoxon.reordered, file = paste(cervical.output.counts.folder, "normalized_filtered_counts.tsv"), quote = FALSE, sep = "\t", col.names = NA)
#
# PlotPCAFunction(samples.dataframe = cervical.normalized.filtered.counts.wilcoxon, design.dataframe = cervical.design.matrix, plot.label = "Cervical PCA")
# PlotPCAFunctionBu(samples.dataframe = cervical.normalized.filtered.counts.wilcoxon, design.dataframe = cervical.design.matrix, plot.label = "Cervical PCA")
# PlotVariancePCA(samples.dataframe = cervical.normalized.filtered.counts.wilcoxon, plot.label = "Cervical screeplot")

# require("car")
# spm(normalized.filtered.counts.wilcoxon.reordered[,c(1:8)])
# which(is.na(normalized.filtered.counts.wilcoxon.reordered))
#
# cervical.de.genes <- DeSeqTimeCourse(counts.dataframe = round(cervical.normalized.filtered.counts.wilcoxon), design.dataframe = cervical.design.matrix)
# WriteDataFrameAsTsv(cervical.de.genes, file.name.path = file.path(cervical.output.de.folder, "cervical_ensgtf_de_genes_complete"))
# # which(is.na(cervical.de.genes))

# write.table(x = significant.cervical.de.genes, file = file.path(cervical.output.de.folder, "de_genes_complete.tsv"), quote = FALSE, sep = "\t", col.names = NA)
#
# ordered.padj.cervical.de.genes <- cervical.de.genes[order(cervical.de.genes$padj),]
# ordered.padj.cervical.de.genes.notna <- ordered.padj.cervical.de.genes[-which(is.na(ordered.padj.cervical.de.genes$padj)),]
# significant.cervical.de.genes <- ordered.padj.cervical.de.genes.notna[ordered.padj.cervical.de.genes.notna$padj<0.05,]
# WriteDataFrameAsTsv(significant.cervical.de.genes, file.name.path = file.path(cervical.output.de.folder, "cervical_ensgtf_de_genes_significant_005_925"))

# write.table(x = significant.cervical.de.genes, file = file.path(cervical.output.de.folder, "significant_de_genes_005_912.tsv"), quote = FALSE, sep = "\t", col.names = NA)
# dim(significant.cervical.de.genes)



# filtered.counts.proportion <- FilterLowCounts(counts.dataframe = cervical.renamed.counts.reordered, design.dataframe = cervical.design.matrix, is.normalized = FALSE, method.type = "Proportion", cpm.cutoff = 1, cv.percentage = 1)


 ##########################
cervical.deseq.dataset.counts <- DESeq2::DESeqDataSetFromMatrix(countData = round(cervical.normalized.filtered.counts.wilcoxon.reordered), colData = cervical.design.matrix, design = ~ Times + Conditions + Times:Conditions)

ENSRNOG00000037449


data<-plotCounts( dds=cervical.deseq.dataset.counts, gene = "Il33", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=cervical.deseq.dataset.counts, gene = "Zbp1", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=cervical.deseq.dataset.counts, gene = "Rspo1", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=cervical.deseq.dataset.counts, gene = "Pdlim1", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=cervical.deseq.dataset.counts, gene = "Gpnmb", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=cervical.deseq.dataset.counts, gene = "C1s", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()
######################################



############################################################################################################## thoracic
thoracic.files <- list.files("Thoracic/thoracic_bam/")
thoracic.bam <- thoracic.files[grep(pattern = "*.bam", thoracic.files)]
thoracic.bam <- thoracic.bam[-grep(pattern = "*.bai", thoracic.bam)]


thoracic.bam.folder <- "/media/dario/dati/time_course/Thoracic/thoracic_bam"
# thoracic.annotation.file <- "/media/dario/dati/time_course/Thoracic/cuffmerge/merged_assembly/merged.gtf"
thoracic.design.file <- file.path("/media/dario/dati/time_course/Thoracic/design_file/thoracic_design_file.txt")

thoracic.bam.files <- list.files(thoracic.bam.folder, pattern="bam$",full.names=TRUE)

thoracic.output.folder <- file.path("/media/dario/dati/time_course/Thoracic", "results")
thoracic.output.counts.folder <- file.path("/media/dario/dati/time_course/Thoracic", "results", "counts", "FeatureCounts")
thoracic.output.de.folder <- file.path("/media/dario/dati/time_course/Thoracic", "results", "de_genes")

dir.create(thoracic.output.counts.folder, recursive = TRUE)
dir.create(thoracic.output.de.folder, recursive = TRUE)

## feature counts
thoracic.ensgtf.counts <- CountBamFilesFeatureCounts(bam.files = thoracic.bam.files, annotation.file = rat5.ensembl.downloaded.gtf, counts.folder = thoracic.output.counts.folder, gtf.attr.type = "gene_id", prefix.output.file = "Thoracic")
WriteDataFrameAsTsv(thoracic.ensgtf.counts, file.name.path = file.path(thoracic.output.counts.folder, "thoracic_ensgtf_counts_FeatureCounts"))
# thoracic.ensgtf.counts<-read.table(file.path(thoracic.output.counts.folder, "thoracic_ensgtf_counts_FeatureCounts.tsv"), header = TRUE, sep = "\t", row.names = 1)
head(thoracic.ensgtf.counts)

thoracic.renamed.counts <- thoracic.ensgtf.counts
colnames(thoracic.renamed.counts) <- c( "thoracic_14d_t1", "thoracic_14d_t2", "thoracic_14d_t3", "thoracic_14d_t4", "thoracic_14d_t5",
                                        "thoracic_14d_u1", "thoracic_14d_u2", "thoracic_14d_u3",
                                        "thoracic_03d_t1", "thoracic_03d_t2", "thoracic_03d_t3", "thoracic_03d_t4", "thoracic_03d_t5",
                                        "thoracic_03d_u1", "thoracic_03d_u2", "thoracic_03d_u3",
                                        "thoracic_56d_t1", "thoracic_56d_t2", "thoracic_56d_t3", "thoracic_56d_t4", "thoracic_56d_t5",
                                        "thoracic_56d_u1",
                                        "thoracic_07d_t1", "thoracic_07d_t2", "thoracic_07d_t3", "thoracic_07d_t4", "thoracic_07d_t5",
                                        "thoracic_07d_u1", "thoracic_07d_u2", "thoracic_07d_u3"
)
# thoracic.renamed.counts.reordered <- thoracic.renamed.counts[ , c( 9:16, 23:30, 1:8, 17:22 )]
thoracic.renamed.counts.reordered <- thoracic.renamed.counts[ , order(colnames(thoracic.renamed.counts))]
colnames(thoracic.renamed.counts.reordered)
head(thoracic.renamed.counts.reordered)
WriteDataFrameAsTsv(thoracic.renamed.counts.reordered, file.name.path = file.path(thoracic.output.counts.folder, "thoracic_ensgtf_reordered_counts"))
# thoracic.renamed.counts.reordered<-read.table(file.path(thoracic.output.counts.folder, "thoracic_ensgtf_reordered_counts.tsv"), header = TRUE, sep = "\t", row.names = 1)

thoracic.design.matrix <- read.table(thoracic.design.file, header = TRUE, row.names = 1)
# thoracic.filtered.counts.wilcoxon <- FilterLowCounts(counts.dataframe = thoracic.renamed.counts.reordered, design.dataframe = thoracic.design.matrix, is.normalized = FALSE, method.type = "Wilcoxon", cpm.cutoff = 1, cv.percentage = 1)
# WriteDataFrameAsTsv(thoracic.filtered.counts.wilcoxon, file.name.path = file.path(thoracic.output.counts.folder, "thoracic_ensgtf_filtered_counts_wilcoxon_16990"))
# # thoracic.filtered.counts.proportion <- FilterLowCounts(counts.dataframe = thoracic.renamed.counts.reordered, design.dataframe = thoracic.design.matrix, is.normalized = FALSE, method.type = "Proportion", cpm.cutoff = 1, cv.percentage = 1)
#
#
# boxplot(log(thoracic.renamed.counts.reordered+1))
# thoracic.filtered.counts.wilcoxon.normalized <- NormalizeData(thoracic.filtered.counts.wilcoxon, norm.type = "uqua")
# WriteDataFrameAsTsv(thoracic.filtered.counts.wilcoxon.normalized, file.name.path = file.path(thoracic.output.counts.folder, "thoracic_ensgtf_filtered_counts_wilcoxon_16990_normalized_upperquartile"))
# boxplot(log(thoracic.filtered.counts.wilcoxon.normalized+1))
# colnames(thoracic.filtered.counts.wilcoxon.normalized)
# thoracic.filtered.counts.wilcoxon.normalized.reordered <- thoracic.filtered.counts.wilcoxon.normalized[ , c( 9:16, 23:30, 1:8, 17:22 )]
# colnames(thoracic.filtered.counts.wilcoxon.normalized.reordered)

# write.table(x = thoracic.filtered.counts.wilcoxon.normalized.reordered, file = paste(thoracic.output.counts.folder, "normalized_filtered_counts.tsv"), quote = FALSE, sep = "\t", col.names = NA)

# PlotPCAFunction(samples.dataframe = thoracic.filtered.counts.wilcoxon.normalized, design.dataframe = thoracic.design.matrix, plot.label = "Thoracic PCA")
# PlotPCAFunctionBu(samples.dataframe = thoracic.filtered.counts.wilcoxon.normalized, design.dataframe = thoracic.design.matrix, plot.label = "Thoracic PCA")
# PlotVariancePCA(samples.dataframe = thoracic.filtered.counts.wilcoxon.normalized, plot.label = "Thoracic screeplot")

# thoracic.de.genes <- DeSeqTimeCourse(counts.dataframe = round(thoracic.filtered.counts.wilcoxon.normalized), design.dataframe = thoracic.design.matrix)
# dim(thoracic.de.genes)
# WriteDataFrameAsTsv(thoracic.filtered.counts.wilcoxon, file.name.path = file.path(thoracic.output.de.folder, "thoracic_ensgtf_de_genes_complete"))
# # write.table(x = thoracic.de.genes, file = file.path(thoracic.output.de.folder, "de_genes_complete.tsv"), quote = FALSE, sep = "\t", col.names = NA)

# ordered.padj.thoracic.de.genes <- thoracic.de.genes[order(thoracic.de.genes$padj),]
# ordered.padj.thoracic.de.genes.notna <- ordered.padj.thoracic.de.genes[-which(is.na(ordered.padj.thoracic.de.genes$padj)),]
# significant.thoracic.de.genes <- ordered.padj.thoracic.de.genes.notna[ordered.padj.thoracic.de.genes.notna$padj<0.05,]
#
# dim(significant.thoracic.de.genes)
# WriteDataFrameAsTsv(thoracic.filtered.counts.wilcoxon, file.name.path = file.path(thoracic.output.de.folder, "thoracic_ensgtf_de_genes_005_321"))
# # write.table(x = significant.thoracic.de.genes, file = file.path(thoracic.output.de.folder, "significant_de_genes_005_322.tsv"), quote = FALSE, sep = "\t", col.names = NA)
#
#



# thoracic.deseq.dataset.counts <- DESeq2::DESeqDataSetFromMatrix(countData = round(thoracic.filtered.counts.wilcoxon.normalized.reordered), colData = thoracic.design.matrix, design = ~ Times + Conditions + Times:Conditions)


# SubstituteTimeLabels <- function(data, ordered.time.labels=NULL, ordered.time.new.labels=NULL) {
#   ## da migliorare con i parametri passati
#   levels(data$Times) <- c(levels(data$Times), "A_3d", "B_7d", "C_14d", "D_56d")
#   data$Times[data$Times =="3d"] <- "A_3d"
#   data$Times[data$Times =="7d"] <- "B_7d"
#   data$Times[data$Times =="14d"] <- "C_14d"
#   data$Times[data$Times =="56d"] <- "D_56d"
#   return(data)
# }

require(DESeq2)
thoracic.deseq.dataset.counts <- DESeq2::DESeqDataSetFromMatrix(countData = round(filtered.counts), colData = design.matrix, design = ~ Times + Conditions + Times:Conditions)

data <- plotCounts( dds=thoracic.deseq.dataset.counts, gene = "ENSRNOG00000007600", intgroup = c("Times", "Conditions"), returnData = TRUE )
# data <- SubstituteTimeLabels(data)
ggp <- ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+ggtitle(paste0("Thoracic ENSRNOG00000007600 gene overtime"))#+scale_y_log10()
ggplotly(ggp)

data<-plotCounts( dds=thoracic.deseq.dataset.counts, gene = "Zbp1", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=thoracic.deseq.dataset.counts, gene = "Rspo1", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=thoracic.deseq.dataset.counts, gene = "Pdlim1", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=thoracic.deseq.dataset.counts, gene = "Gpnmb", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=thoracic.deseq.dataset.counts, gene = "C1s", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=thoracic.deseq.dataset.counts, gene = "Mapk4", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()


data<-plotCounts( dds=thoracic.deseq.dataset.counts, gene = "Zfp865", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

data<-plotCounts( dds=thoracic.deseq.dataset.counts, gene = "Epn1", intgroup = c("Times", "Conditions"), returnData = TRUE )
data <- SubstituteTimeLabels(data)
ggplot(data, mapping = aes(x=Times, y=count, color=Conditions, group=Conditions))+geom_point()+stat_smooth(se=F, method="loess")+scale_y_log10()

############################################# scorporating thoracic by timesteps
# thoracic.counts.03 <- thoracic.renamed.counts.reordered[,c(1:8)]
# thoracic.counts.07 <- thoracic.renamed.counts.reordered[,c(9:16)]
# thoracic.counts.14 <- thoracic.renamed.counts.reordered[,c(17:24)]
# thoracic.counts.56 <- thoracic.renamed.counts.reordered[,c(25:30)]
#
# thoracic.design.matrix.03 <- thoracic.design.matrix[c(1:8) ,]
# thoracic.design.matrix.07 <- thoracic.design.matrix[c(9:16) ,]
# thoracic.design.matrix.14 <- thoracic.design.matrix[c(17:24) ,]
# thoracic.design.matrix.56 <- thoracic.design.matrix[c(25:30) ,]
#
# thoracic.output.folder.03 <- paste0(thoracic.output.folder,"_03")
# thoracic.de.genes.03 <- MainFunction(counts.dataframe = thoracic.counts.03, design.matrix = thoracic.design.matrix.03, output.results.folder = thoracic.output.folder.03, prefix.output.file = "thoracic_03", prefix.plot.label = "thoracic 03", is.time.course = FALSE)
# nn.genes=20
# thoracic.de.genes.03.tail <- TailDeGenesByPAdj(thoracic.de.genes.03, n.genes = nn.genes)
# ## su questi geni si deve fare una selezione delle distribuzioni in modo da ottenerle tutte o quasi diverse
#
# thoracic.dataset.set.03 <- newSeqExpressionSet( as.matrix(thoracic.counts.03),
#                                              phenoData= data.frame(
#                                                thoracic.design.matrix.03$Conditions, row.names = rownames(thoracic.design.matrix.03)
#                                              ))
#
#
# ruved.thoracic.set.03 <- RUVg(thoracic.dataset.set.03, thoracic.de.genes.03.tail, k=5)
# plotRLE(ruved.thoracic.set.03, outline=FALSE, ylim=c(-4,4), main=as.character(nn.genes))
# plotPCA(ruved.thoracic.set.03, cex=0.8)

############################
thoracic.renamed.counts <- thoracic.renamed.counts.reordered
thoracic.de.genes.wilcoxon.uqua <- MainFunction(counts.dataframe = thoracic.renamed.counts, design.matrix = thoracic.design.matrix, output.results.folder = thoracic.output.folder, prefix.output.file = "thoracic", prefix.plot.label = "thoracic", is.time.course = TRUE, filter.method = "Wilcoxon")
WriteDataFrameAsTsv(thoracic.de.genes.wilcoxon.uqua, file.name.path = file.path(thoracic.output.de.folder, "thoracic_de_genes_complete_wilcoxon_uqua"))

thoracic.de.genes.wilcoxon.ordered <- SortDeGenesByPAdj(thoracic.de.genes.wilcoxon.uqua)

thoracic.de.genes.proportion.uqua <- MainFunction(counts.dataframe = thoracic.renamed.counts, design.matrix = thoracic.design.matrix, output.results.folder = thoracic.output.folder, prefix.output.file = "thoracic", prefix.plot.label = "thoracic", is.time.course = TRUE, filter.method = "Proportion")
WriteDataFrameAsTsv(thoracic.de.genes.proportion.uqua, file.name.path = file.path(thoracic.output.de.folder, "thoracic_de_genes_complete_proportion_uqua"))
thoracic.de.genes.proportion.ordered <- SortDeGenesByPAdj(thoracic.de.genes.proportion.uqua)

thoracic.sign.de.genes.wilcoxon <- SignificantDeGenesPAdj(thoracic.de.genes.wilcoxon.ordered)
WriteDataFrameAsTsv(thoracic.sign.de.genes.wilcoxon, file.name.path = "/media/dario/dati/time_course/Thoracic/results/de_genes/wilcoxon/uqua/thoracic__filtered_counts__Wilcoxon__normalized__uqua_de_genes_LRT_significant_005.tsv")
thoracic.sign.de.genes.proportion <- SignificantDeGenesPAdj(thoracic.de.genes.proportion.ordered)
WriteDataFrameAsTsv(thoracic.sign.de.genes.proportion, file.name.path = "/media/dario/dati/time_course/Thoracic/results/de_genes/proportion/uqua/thoracic__filtered_counts__proportion__normalized__uqua_de_genes_LRT_significant_005.tsv")


thoracic.proportion.estimated.neg.genes <- EstimateNegativeControlGenesForRUV(de.genes = thoracic.de.genes.proportion.ordered, n.tail.genes = 2000, counts.dataset = thoracic.renamed.counts, n.genes.per.hist.break = 1)
# thoracic.proportion.unnormalized<-read.table("/media/dario/dati/time_course/Thoracic/results_backup_18_01_2016/counts/FeatureCounts/thoracic__filtered_counts__Proportion.tsv", sep="\t", row.names = 1, header = TRUE)
# ReadDataFrameFromTsv("/media/dario/dati/time_course/Thoracic/results/counts/FeatureCounts/thoracic_ensgtf_reordered_counts.tsv")
thoracic.de.genes.wilcoxon.ruvg <- MainFunction(counts.dataframe = thoracic.renamed.counts, design.matrix = thoracic.design.matrix, output.results.folder = thoracic.output.folder, prefix.output.file = "thoracic", prefix.plot.label = "thoracic", is.time.course = TRUE, filter.method = "Wilcoxon", normalization.method = "ruvg", estimated.genes = thoracic.proportion.estimated.neg.genes)

thoracic.de.genes.wilcoxon.ruvg.significant <- SignificantDeGenesPAdj(thoracic.de.genes.wilcoxon.ruvg)
WriteDataFrameAsTsv(thoracic.de.genes.wilcoxon.ruvg.significant, file.name.path = "/media/dario/dati/time_course/Thoracic/results/de_genes/wilcoxon/ruvg/thoracic__filtered_counts__Wilcoxon__normalized__ruvg_de_genes_LRT_significant_005.tsv")

thoracic.de.genes.proportion.ruvg <- MainFunction(counts.dataframe = thoracic.renamed.counts, design.matrix = thoracic.design.matrix, output.results.folder = thoracic.output.folder, prefix.output.file = "thoracic", prefix.plot.label = "thoracic", is.time.course = TRUE, filter.method = "Proportion", normalization.method = "ruvg", estimated.genes = thoracic.proportion.estimated.neg.genes)

thoracic.de.genes.proportion.ruvg.significant <- SignificantDeGenesPAdj(thoracic.de.genes.proportion.ruvg)
WriteDataFrameAsTsv(thoracic.de.genes.proportion.ruvg.significant, file.name.path = "/media/dario/dati/time_course/Thoracic/results/de_genes/proportion/ruvg/thoracic__filtered_counts__Proportion__normalized__ruvg_de_genes_LRT_significant_005.tsv")

# thoracic.proportion.counts.normalized.ruv <- NormalizeData(data.to.normalize = thoracic.proportion.unnormalized, norm.type = "ruvg", design.matrix = thoracic.design.matrix, estimated.genes = thoracic.proportion.estimated.neg.genes)
# PlotPCAPlotlyFunction(counts.data.frame = thoracic.proportion.counts.normalized.ruv, design.matrix = thoracic.design.matrix, colour.design.column.str = "Times", shape.design.column.str = "Conditions")

# de.genes = thoracic.de.genes.proportion.ordered; n.tail.genes = 2000; counts.dataset = thoracic.renamed.counts; n.genes.per.hist.break = 1

#
# dim(thoracic.sign.de.genes.wilcoxon)
# dim(thoracic.sign.de.genes.proportion)
# length(which(rownames(thoracic.sign.de.genes.proportion) %in% rownames(thoracic.sign.de.genes.wilcoxon)))
#
# nn.genes=2000
# thoracic.de.genes.tail.proportion <- TailDeGenesByPAdj(thoracic.de.genes.proportion.ordered, n.genes = nn.genes)
#
# thoracic.proportion.counts.tail <- thoracic.renamed.counts[which(rownames(thoracic.renamed.counts) %in% thoracic.de.genes.tail.proportion),]
#
# sd(thoracic.proportion.counts.tail[1,])
# mean(as.integer(thoracic.proportion.counts.tail[1,]))
#
# thoracic.proportion.mean.tail <- apply(thoracic.proportion.counts.tail, 1, mean)
# thoracic.proportion.sd.tail <- apply(thoracic.proportion.counts.tail, 1, sd)
#
# sort(thoracic.proportion.mean.tail)
#
# tt<-thoracic.proportion.mean.tail[which(thoracic.proportion.mean.tail <= summary(thoracic.proportion.mean.tail)[5])]
#
#
# tsd <- thoracic.proportion.sd.tail[which(names(thoracic.proportion.sd.tail) %in% names(tt))]
#
# tt
#
# hist(summary(t(as.matrix(thoracic.proportion.counts.tail))), breaks = 100, freq = TRUE, xlim = c(0, 4000))
#
# str(summary(t(as.matrix(thoracic.proportion.counts.tail))))



#
# est.genes <- SelectGenesFromHistBreaks(thoracic.proportion.counts.tail)
#
# require("RUVSeq")
# thoracic.dataset.set <- newSeqExpressionSet( as.matrix(thoracic.renamed.counts.reordered),
#                                              phenoData= data.frame(
#                                                thoracic.design.matrix$Conditions, row.names = rownames(thoracic.design.matrix)
#                                              ))
#
#
# ruved.thoracic.set <- RUVg(thoracic.dataset.set, est.genes, k=1)
# plotRLE(ruved.thoracic.set, outline=FALSE, ylim=c(-4,4))
# PlotPCAFunction(ruved.thoracic.set@assayData$normalizedCounts, design.dataframe = thoracic.design.matrix, scale = FALSE)
#


# nn.genes=2000
# thoracic.de.genes.tail.wilcoxon <- TailDeGenesByPAdj(thoracic.de.genes.wilcoxon.ordered, n.genes = nn.genes)
#
# thoracic.wilcoxon.counts.tail <- thoracic.renamed.counts[which(rownames(thoracic.renamed.counts) %in% thoracic.de.genes.tail.wilcoxon),]
# est.genes.wilcoxon <- SelectGenesFromHistBreaks(thoracic.wilcoxon.counts.tail, n.genes.per.break=1)
#
# require("RUVSeq")
# thoracic.dataset.set <- newSeqExpressionSet( as.matrix(thoracic.renamed.counts.reordered),
#                                              phenoData= data.frame(
#                                                thoracic.design.matrix$Conditions, row.names = rownames(thoracic.design.matrix)
#                                              ))
#
#
# ruved.thoracic.set <- RUVg(thoracic.dataset.set, est.genes.wilcoxon, k=1)

#
#
# plotRLE(ruved.thoracic.set, outline=FALSE, ylim=c(-4,4))
# PlotPCAFunction(ruved.thoracic.set@assayData$normalizedCounts, design.dataframe = thoracic.design.matrix, scale = FALSE)

thoracic.output.plots.folder <- file.path(thoracic.output.folder, "plots")

# wilcoxon - ruv
th.w.rg<-ReadDataFrameFromTsv("/media/dario/dati/time_course/Thoracic/results/counts/FeatureCounts/wilcoxon/ruvg/thoracic__filtered_counts__Wilcoxon__normalized__ruvg.tsv")

counts.data.frame = th.w.rg; design.matrix = thoracic.design.matrix; colour.design.column.str = "Times"; shape.design.column.str = "Conditions";
output.path = thoracic.output.plots.folder; prefix.plot = "thoracic_wilcoxon_ruvg"; title="Thoracic Normalized RUVg PCA"

# proportion - ruv
th.p.rg<-ReadDataFrameFromTsv("/media/dario/dati/time_course/Thoracic/results/counts/FeatureCounts/proportion/ruvg/thoracic__filtered_counts__Proportion__normalized__ruvg.tsv")

counts.data.frame = th.p.rg; design.matrix = thoracic.design.matrix; colour.design.column.str = "Times"; shape.design.column.str = "Conditions";
output.path = thoracic.output.plots.folder; prefix.plot = "thoracic_proportion_ruvg"; title="Thoracic Normalized RUVg PCA"


# wilcoxon - uqua
th.w.uqua<-ReadDataFrameFromTsv("/media/dario/dati/time_course/Thoracic/results/counts/FeatureCounts/thoracic__filtered_counts__Wilcoxon__normalized__uqua.tsv")

counts.data.frame = th.w.uqua; design.matrix = thoracic.design.matrix; colour.design.column.str = "Times"; shape.design.column.str = "Conditions";
output.path = thoracic.output.plots.folder; prefix.plot = "thoracic_wilcoxon_uqua"; title="Thoracic Normalized UQUA PCA"

# proportion - uqua
th.p.uqua<-ReadDataFrameFromTsv("/media/dario/dati/time_course/Thoracic/results/counts/FeatureCounts/proportion/uqua/thoracic__filtered_counts__Proportion__normalized__uqua.tsv")

counts.data.frame = th.p.uqua; design.matrix = thoracic.design.matrix; colour.design.column.str = "Times"; shape.design.column.str = "Conditions";
output.path = thoracic.output.plots.folder; prefix.plot = "thoracic_proportion_uqua"; title="Thoracic Normalized UQUA PCA"


