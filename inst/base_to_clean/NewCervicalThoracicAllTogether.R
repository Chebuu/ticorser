### all conditions 
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/includes.R")

cerv.counts <- ReadDataFrameFromTsv("Cervical/Cervical_results/counts/FeatureCounts/cervical_ensgtf_reordered_counts.tsv")
thor.counts <- ReadDataFrameFromTsv("Thoracic/Thoracic_results/counts/FeatureCounts/thoracic_ensgtf_reordered_counts.tsv")
all.counts <- cbind(cerv.counts, thor.counts)

all.design.file <- "cervical_thoracic/design_file/cervical_thoracic_all_design_file.txt"
cer.thor.des.mat <- ReadDataFrameFromTsv(all.design.file)

output.results.folder <- "cervical_thoracic/results_all"


NewMainFunction <- function(counts.dataframe, design.matrix, output.results.folder, prefix.output.file, prefix.plot.label, filter.method="Proportion", normalization.method= "uqua", is.time.course=TRUE, estimated.genes=NULL, ellipse.in.pca=FALSE) {
  
  output.counts.folder <- file.path(output.results.folder, "counts", "FeatureCounts")
  output.de.folder <- file.path(output.results.folder, "de_genes", "DESeq2")
  output.plots.folder <- file.path(output.results.folder, "plots")
  
  dir.create(output.counts.folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(output.de.folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(output.plots.folder, recursive = TRUE, showWarnings = FALSE)
  
  filename <- paste0(prefix.output.file)
  
  message("Filtering counts with ", filter.method, " method\n")
  filtered.counts <- FilterLowCounts(counts.dataframe = counts.dataframe, design.dataframe = design.matrix, is.normalized = FALSE, method.type = filter.method, cpm.cutoff = 1, cv.percentage = 1)
  filename <- UpdateFilename(filename, "_filtered_counts_", filter.method)
  WriteDataFrameAsTsv(filtered.counts, file.name.path = file.path(output.counts.folder, filename))
  PlotTimesBoxplot(data.frame.to.plot = filtered.counts, design.matrix = design.matrix, output.path = output.plots.folder, prefix.plot = prefix.plot.label, show.plot.flag = FALSE, save.plot = TRUE, plotly.flag = TRUE)
  PlotPCAPlotlyFunction(counts.data.frame = filtered.counts, design.matrix = design.matrix, scale = FALSE, plot.folder = output.plots.folder, prefix.plot = filename, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.name = "Cerv-Thor un-normalized PCA", ellipse.flag = ellipse.in.pca)
  
  message("Normalizing counts with ", normalization.method, "\n")
  normalized.filtered.counts <- NormalizeData(filtered.counts, norm.type = normalization.method, estimated.genes = estimated.genes, design.matrix = design.matrix)
  filename <- UpdateFilename(filename, "_normalized_", normalization.method)
  WriteDataFrameAsTsv(normalized.filtered.counts, file.name.path = file.path(output.counts.folder, filename))
  
  if(normalization.method=="ruvg") normalized.filtered.counts <- as.data.frame(normalized.filtered.counts)
  
  PlotTimesBoxplot(data.frame.to.plot = normalized.filtered.counts, design.matrix = design.matrix, output.path = output.plots.folder, prefix.plot = paste0(normalization.method, " boxplot"), show.plot.flag = FALSE, save.plot = TRUE, plotly.flag = TRUE)
  
  PlotPCAPlotlyFunction(counts.data.frame = normalized.filtered.counts, design.matrix = design.matrix, scale = FALSE, plot.folder = output.plots.folder, prefix.plot = filename, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.name=paste0(normalization.method," PCA"), ellipse.flag = ellipse.in.pca)
  
  
}


NewMainFunction(counts.dataframe = all.counts, design.matrix = cer.thor.des.mat, output.results.folder = output.results.folder, prefix.output.file = "cerv_thor_all_uqua", prefix.plot.label = "cerv_thor_all", filter.method = "Proportion", normalization.method = "uqua")

## cerv
cerv.thor.prop.uqua <- ReadDataFrameFromTsv("cervical_thoracic/results_all/counts/FeatureCounts/cerv_thor_all__filtered_counts__Proportion__normalized__uqua.tsv")
cerv.des.matrix <- ReadDataFrameFromTsv("Cervical/design_file/cervical_design_file.txt")
cerv.all.prop.uqua <- cerv.thor.prop.uqua[,(colnames(cerv.thor.prop.uqua) %in% rownames(cerv.des.matrix))]
colnames(cerv.all.prop.uqua)
cerv.de.genes <- DeSeqTimeCourse(counts.dataframe = round(cerv.all.prop.uqua), design.dataframe = cerv.des.matrix)
head(cerv.de.genes)
WriteDataFrameAsTsv(file.name.path = "cervical_thoracic/results_all/de_genes/DESeq2/cervical_all_prop_uqua_deseq_de_genes.tsv", data.frame.to.save = cerv.de.genes)
cerv.sign.de.genes <- SignificantDeGenesPAdj(cerv.de.genes)
WriteDataFrameAsTsv(file.name.path = "cervical_thoracic/results_all/de_genes/DESeq2/cervical_all_prop_uqua_deseq_de_genes_sign_005.tsv", data.frame.to.save = cerv.sign.de.genes)
cerv.est.genes <- EstimateNegativeControlGenesForRUV(de.genes = cerv.de.genes, n.tail.genes = 2000, counts.dataset = cerv.counts)

## thor
thor.des.matrix <- ReadDataFrameFromTsv("Thoracic/design_file/thoracic_design_file.txt")
thor.all.prop.uqua <- cerv.thor.prop.uqua[,(colnames(cerv.thor.prop.uqua) %in% rownames(thor.des.matrix))]
colnames(thor.all.prop.uqua)
thor.de.genes <- DeSeqTimeCourse(counts.dataframe = round(thor.all.prop.uqua), design.dataframe = thor.des.matrix)
WriteDataFrameAsTsv(file.name.path = "cervical_thoracic/results_all/de_genes/DESeq2/thoracic_all_prop_uqua_deseq_de_genes.tsv", data.frame.to.save = thor.de.genes)
thor.sign.de.genes <- SignificantDeGenesPAdj(thor.de.genes)
dim(thor.sign.de.genes)
WriteDataFrameAsTsv(file.name.path = "cervical_thoracic/results_all/de_genes/DESeq2/thoracic_all_prop_uqua_deseq_de_genes_sign_005.tsv", data.frame.to.save = thor.sign.de.genes)
thor.est.genes <- EstimateNegativeControlGenesForRUV(de.genes = thor.de.genes, n.tail.genes = 2000, counts.dataset = thor.counts)

## cerv-thor-t
cerv.thor.t.des.matrix <- ReadDataFrameFromTsv("cervical_thoracic/design_file/cervical_thoracic_design_file.txt")
cerv.thor.t.prop.uqua <- cerv.thor.prop.uqua[,(colnames(cerv.thor.prop.uqua) %in% rownames(cerv.thor.t.des.matrix))]
cerv.thor.t.de.genes <- DeSeqTimeCourse(counts.dataframe = round(cerv.thor.t.prop.uqua), design.dataframe = cerv.thor.t.des.matrix)
WriteDataFrameAsTsv(file.name.path = "cervical_thoracic/results_all/de_genes/DESeq2/cerv_thor_all_prop_uqua_deseq_de_genes.tsv", data.frame.to.save = cerv.thor.t.de.genes)
cerv.thor.t.sign.de.genes <- SignificantDeGenesPAdj(cerv.thor.t.de.genes)
dim(cerv.thor.t.sign.de.genes)
WriteDataFrameAsTsv(file.name.path = "cervical_thoracic/results_all/de_genes/DESeq2/cerv_thor_all_prop_uqua_deseq_de_genes_sign_005.tsv", data.frame.to.save = cerv.thor.t.sign.de.genes)
cerv.thor.t.counts <- all.counts[,(colnames(all.counts) %in% rownames(cerv.thor.t.des.matrix))]
cerv.thor.est.genes <- EstimateNegativeControlGenesForRUV(de.genes = cerv.thor.t.de.genes, n.tail.genes = 2000, counts.dataset = cerv.thor.t.counts)
####

Venn3de(x = cerv.est.genes, y = thor.est.genes, z = cerv.thor.est.genes, label1 = "Cervical", label2 = "Thoracic", label3 = "Cerv-Thor", title = "Venn of Negative Estimated Genes", intersection.flag = TRUE, intersection.exclusion.flag = FALSE, plot.dir = "cervical_thoracic/results_all/plots/")

cerv.est.genes.2 <- EstimateNegativeControlGenesForRUV(de.genes = cerv.de.genes, n.tail.genes = 2000, counts.dataset = cerv.counts, n.genes.per.hist.break = 2)
thor.est.genes.2 <- EstimateNegativeControlGenesForRUV(de.genes = thor.de.genes, n.tail.genes = 2000, counts.dataset = thor.counts, n.genes.per.hist.break = 2)
cerv.thor.est.genes.2 <- EstimateNegativeControlGenesForRUV(de.genes = cerv.thor.t.de.genes, n.tail.genes = 2000, counts.dataset = cerv.thor.t.counts, n.genes.per.hist.break = 2)
Venn3de(x = cerv.est.genes.2, y = thor.est.genes.2, z = cerv.thor.est.genes.2, label1 = "Cervical2", label2 = "Thoracic2", label3 = "Cerv-Thor2", title = "Venn of Negative Estimated Genes 2", intersection.flag = TRUE, intersection.exclusion.flag = FALSE, plot.dir = "cervical_thoracic/results_all/plots/")


cerv.est.genes.10 <- EstimateNegativeControlGenesForRUV(de.genes = cerv.de.genes, n.tail.genes = 2000, counts.dataset = cerv.counts, n.genes.per.hist.break = 10)
thor.est.genes.10 <- EstimateNegativeControlGenesForRUV(de.genes = thor.de.genes, n.tail.genes = 2000, counts.dataset = thor.counts, n.genes.per.hist.break = 10)
cerv.thor.est.genes.10 <- EstimateNegativeControlGenesForRUV(de.genes = cerv.thor.t.de.genes, n.tail.genes = 2000, counts.dataset = cerv.thor.t.counts, n.genes.per.hist.break = 10)
Venn3de(x = cerv.est.genes.10, y = thor.est.genes.10, z = cerv.thor.est.genes.10, label1 = "Cervical10", label2 = "Thoracic10", label3 = "Cerv-Thor10", title = "Venn of Negative Estimated Genes 10", intersection.flag = TRUE, intersection.exclusion.flag = FALSE, plot.dir = "cervical_thoracic/results_all/plots/")

cerv.est.genes.30 <- EstimateNegativeControlGenesForRUV(de.genes = cerv.de.genes, n.tail.genes = 2000, counts.dataset = cerv.counts, n.genes.per.hist.break = 30)
thor.est.genes.30 <- EstimateNegativeControlGenesForRUV(de.genes = thor.de.genes, n.tail.genes = 2000, counts.dataset = thor.counts, n.genes.per.hist.break = 30)
cerv.thor.est.genes.30 <- EstimateNegativeControlGenesForRUV(de.genes = cerv.thor.t.de.genes, n.tail.genes = 2000, counts.dataset = cerv.thor.t.counts, n.genes.per.hist.break = 30)
Venn3de(x = cerv.est.genes.30, y = thor.est.genes.30, z = cerv.thor.est.genes.30, label1 = "Cervical30", label2 = "Thoracic30", label3 = "Cerv-Thor30", title = "Venn of Negative Estimated Genes 30", intersection.flag = TRUE, intersection.exclusion.flag = FALSE, plot.dir = "cervical_thoracic/results_all/plots/")

neg.est.genes.inters <- ReadDataFrameFromTsv("cervical_thoracic/results_all/plots/Venn3/gene_lists/Cervical30_Thoracic30_Cerv-Thor30_genes_in_intersection.txt", row.names.col = NULL)

neg.est.genes.inters<-neg.est.genes.inters[-is.na(neg.est.genes.inters)]
## ruv normalization
NewMainFunction(counts.dataframe = all.counts, design.matrix = cer.thor.des.mat, output.results.folder = output.results.folder, prefix.output.file = "cerv_thor_all_ruvg", prefix.plot.label = "cerv_thor_all", filter.method = "Proportion", normalization.method = "ruvg", estimated.genes = neg.est.genes.inters)


cerv.est.genes.3000.20 <- EstimateNegativeControlGenesForRUV(de.genes = cerv.de.genes, n.tail.genes = 3000, counts.dataset = cerv.counts, n.genes.per.hist.break = 20)
thor.est.genes.3000.20 <- EstimateNegativeControlGenesForRUV(de.genes = thor.de.genes, n.tail.genes = 3000, counts.dataset = thor.counts, n.genes.per.hist.break = 20)
cerv.thor.est.genes.3000.20 <- EstimateNegativeControlGenesForRUV(de.genes = cerv.thor.t.de.genes, n.tail.genes = 3000, counts.dataset = cerv.thor.t.counts, n.genes.per.hist.break = 20)
Venn3de(x = cerv.est.genes.3000.20, y = thor.est.genes.3000.20, z = cerv.thor.est.genes.3000.20, label1 = "Cervical3k.20", label2 = "Thoracic3k.20", label3 = "Cerv-Thor3k.20", title = "Venn of Negative Estimated Genes 3000/20", intersection.flag = TRUE, intersection.exclusion.flag = FALSE, plot.dir = "cervical_thoracic/results_all/plots/")


neg.est.genes.inters.3k.20 <- ReadDataFrameFromTsv("cervical_thoracic/results_all/plots/Venn3/gene_lists/Cervical3k.20_Thoracic3k.20_Cerv-Thor3k.20_genes_in_intersection.txt", row.names.col = NULL)
neg.est.genes.inters.3k.20<-neg.est.genes.inters.3k.20[-is.na(neg.est.genes.inters.3k.20)]
## ruv normalization
NewMainFunction(counts.dataframe = all.counts, design.matrix = cer.thor.des.mat, output.results.folder = output.results.folder, prefix.output.file = "cerv_thor_all_ruvg", prefix.plot.label = "cerv_thor_all", filter.method = "Proportion", normalization.method = "ruvg", estimated.genes = neg.est.genes.inters.3k.20)



cerv.est.genes.10k.5 <- EstimateNegativeControlGenesForRUV(de.genes = cerv.de.genes, n.tail.genes = 10000, counts.dataset = cerv.counts, n.genes.per.hist.break = 5)
thor.est.genes.10k.5 <- EstimateNegativeControlGenesForRUV(de.genes = thor.de.genes, n.tail.genes = 10000, counts.dataset = thor.counts, n.genes.per.hist.break = 5)
cerv.thor.est.genes.10k.5 <- EstimateNegativeControlGenesForRUV(de.genes = cerv.thor.t.de.genes, n.tail.genes = 10000, counts.dataset = cerv.thor.t.counts, n.genes.per.hist.break = 5)
Venn3de(x = cerv.est.genes.10k.5, y = thor.est.genes.10k.5, z = cerv.thor.est.genes.10k.5, label1 = "Cervical10k.5", label2 = "Thoracic10k.5", label3 = "Cerv-Thor10k.5", title = "Venn of Negative Estimated Genes 3000/20", intersection.flag = TRUE, intersection.exclusion.flag = FALSE, plot.dir = "cervical_thoracic/results_all/plots/")
