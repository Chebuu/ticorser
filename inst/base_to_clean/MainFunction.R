

MainFunction <- function(counts.dataframe, design.matrix, output.results.folder, prefix.output.file, prefix.plot.label, filter.method="Wilcoxon", normalization.method= "uqua", is.time.course=TRUE, estimated.genes=NULL) {
  
  
  output.counts.folder <- file.path(output.results.folder, "counts", "FeatureCounts")
  output.de.folder <- file.path(output.results.folder, "de_genes", "DESeq2")
  output.plots.folder <- file.path(output.results.folder, "plots")
  
  dir.create(output.counts.folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(output.de.folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(output.plots.folder, recursive = TRUE, showWarnings = FALSE)
  
  filename <- paste0(prefix.output.file)
  
  cat("Filtering counts with ", filter.method,"method\n")
  filtered.counts <- FilterLowCounts(counts.dataframe = counts.dataframe, design.dataframe = design.matrix, is.normalized = FALSE, method.type = filter.method, cpm.cutoff = 1, cv.percentage = 1)
  filename <- UpdateFilename(filename, "_filtered_counts_", filter.method)
  WriteDataFrameAsTsv(filtered.counts, file.name.path = file.path(output.counts.folder, filename))

  jpeg(file.path(output.plots.folder, paste0(prefix.plot.label, "_", filter.method,"_unnormalized_boxplot.jpg")))
  boxplot(log(filtered.counts+1), main=paste0(prefix.plot.label, " unnormalized counts"), las=2)
  dev.off()
  
  # PlotPCAPlotlyFunction(counts.data.frame = filtered.counts, design.matrix = design.matrix, colour.design.column.str = "Times", shape.design.column.str = "Conditions", output.path = output.plots.folder, prefix.plot = prefix.plot.label, title="Unnormalized PCA")
  # counts.data.frame = filtered.counts; design.matrix = design.matrix; colour.design.column.str = "Times"; shape.design.column.str = "Conditions"; output.path = output.plots.folder; prefix.plot = prefix.plot.label; title="Unnormalized PCA"
  
  cat("Normalizing counts with ", normalization.method, "\n")
  normalized.filtered.counts <- NormalizeData(filtered.counts, norm.type = normalization.method, estimated.genes = estimated.genes, design.matrix = design.matrix)
  filename <- UpdateFilename(filename, "_normalized_", normalization.method)
  WriteDataFrameAsTsv(normalized.filtered.counts, file.name.path = file.path(output.counts.folder, filename))
  
  # PlotPCAPlotlyFunction(counts.data.frame = normalized.filtered.counts, design.matrix = design.matrix, colour.design.column.str = "Times", shape.design.column.str = "Conditions", output.path = output.plots.folder, prefix.plot = prefix.plot.label, title=paste0("Normalized ", normalization.method," PCA"))
  
  jpeg(file.path(output.plots.folder, paste0(prefix.plot.label, "_", filter.method, "_",  normalization.method,"_boxplot.jpg")))
  boxplot(log(normalized.filtered.counts+1), main=paste0(prefix.plot.label, " normalized ", normalization.method," counts"), las=2)
  dev.off()
  ####################################
  # names <- rownames(thoracic.design.matrix)
  # fill <- thoracic.design.matrix$Times
  # df <- as.data.frame(stack(log(ruved.thoracic.set@assayData$normalizedCounts+1)))
  # ggplot(df, aes(col, )) + geom_boxplot()
  
  # theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=8)) + guides(fill=guide_legend(title=NULL)) + ggtitle(title)
  
  ##########################
  # PlotPCAFunction(samples.dataframe = normalized.filtered.counts, design.dataframe = design.matrix, plot.label = paste0(prefix.plot.label, " normalized ", normalization.method," PCA"))
# PlotPCAFunctionBu(samples.dataframe = cervical.normalized.filtered.counts.wilcoxon, design.dataframe = design.matrix, plot.label = "Cervical PCA")
# PlotVariancePCA(samples.dataframe = cervical.normalized.filtered.counts.wilcoxon, plot.label = "Cervical screeplot")
# 
# # require("car")
# # spm(normalized.filtered.counts.wilcoxon.reordered[,c(1:8)])
# # which(is.na(normalized.filtered.counts.wilcoxon.reordered))
  
  if(is.time.course) { 
    test.type="LRT"
    de.genes <- DeSeqTimeCourse(counts.dataframe = round(normalized.filtered.counts), design.dataframe = design.matrix)
    filename <- UpdateFilename(filename, "de_genes_LRT")
  } else {
    test.type="Wald"
    de.genes <- DeSeqTimeCourse(counts.dataframe = round(normalized.filtered.counts), design.dataframe = design.matrix, test.type="Wald")
    filename <- UpdateFilename(filename, "de_genes_wald")
  }
  
  WriteDataFrameAsTsv(de.genes, file.name.path = file.path(output.de.folder, filename))
  
  return(de.genes)
# 
# WriteDataFrameAsTsv(cervical.de.genes, file.name.path = file.path(cervical.output.de.folder, "cervical_ensgtf_de_genes_complete"))
# # which(is.na(cervical.de.genes))
# 
# # write.table(x = significant.cervical.de.genes, file = file.path(cervical.output.de.folder, "de_genes_complete.tsv"), quote = FALSE, sep = "\t", col.names = NA)
# 
# ordered.padj.cervical.de.genes <- cervical.de.genes[order(cervical.de.genes$padj),]
# ordered.padj.cervical.de.genes.notna <- ordered.padj.cervical.de.genes[-which(is.na(ordered.padj.cervical.de.genes$padj)),]
# significant.cervical.de.genes <- ordered.padj.cervical.de.genes.notna[ordered.padj.cervical.de.genes.notna$padj<0.05,]
# WriteDataFrameAsTsv(significant.cervical.de.genes, file.name.path = file.path(cervical.output.de.folder, "cervical_ensgtf_de_genes_significant_005_925"))
# 
# # write.table(x = significant.cervical.de.genes, file = file.path(cervical.output.de.folder, "significant_de_genes_005_912.tsv"), quote = FALSE, sep = "\t", col.names = NA)
# dim(significant.cervical.de.genes)

}


# MainFunction(counts.dataframe = thoracic.renamed.counts.reordered, design.matrix = thoracic.design.matrix, output.results.folder = thoracic.output.folder, prefix.output.file = "thoracic", prefix.plot.label = "thoracic")
