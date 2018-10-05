## thoracic first 3 time points

thoracic.complete.counts <- ReadDataFrameFromTsv("Thoracic/results/counts/FeatureCounts/thoracic_ensgtf_reordered_counts.tsv")
colnames(thoracic.complete.counts)
thoracic.3.14.counts <- thoracic.complete.counts[, c(1:24)]
colnames(thoracic.3.14.counts)

thoracic.design.file <- file.path("/media/dario/dati/time_course/Thoracic/design_file/thoracic_design_file.txt")
thoracic.design.matrix <- ReadDataFrameFromTsv(thoracic.design.file)
thoracic.3.14.design.matrix <- thoracic.design.matrix[-which(thoracic.design.matrix$Times %in% "56d"),]

thoracic.output.folder.3.14 <- file.path("/media/dario/dati/time_course/Thoracic", "results_3_14")

thoracic.3.14.de.genes.proportion.uqua <- MainFunction(counts.dataframe = thoracic.3.14.counts, design.matrix = thoracic.3.14.design.matrix, output.results.folder = thoracic.output.folder.3.14, prefix.output.file = "thoracic_3_14", prefix.plot.label = "thoracic_3_14", is.time.course = TRUE, filter.method = "Proportion")

thoracic.3.14.de.genes.proportion.uqua.sign.005 <- SignificantDeGenesPAdj(thoracic.3.14.de.genes.proportion.uqua)

WriteDataFrameAsTsv(thoracic.3.14.de.genes.proportion.uqua.sign.005, file.name.path = "/media/dario/dati/time_course/Thoracic/results_3_14/de_genes/thoracic_3_14_filtered_counts__Proportion__normalized__uqua_de_genes_LRT_significant_005.tsv")

thoracic.3.14.proportion.estimated.neg.genes <- EstimateNegativeControlGenesForRUV(de.genes = thoracic.3.14.de.genes.proportion.uqua, n.tail.genes = 2000, counts.dataset = thoracic.3.14.counts, n.genes.per.hist.break = 1)

thoracic.3.14.de.genes.proportion.ruvg <- MainFunction(counts.dataframe = thoracic.3.14.counts, design.matrix = thoracic.3.14.design.matrix, output.results.folder = thoracic.output.folder.3.14, prefix.output.file = "thoracic_3_14", prefix.plot.label = "thoracic_3_14", is.time.course = TRUE, filter.method = "Proportion", normalization.method = "ruvg", estimated.genes = thoracic.3.14.proportion.estimated.neg.genes)
thoracic.3.14.de.genes.proportion.ruvg.sign.005 <- SignificantDeGenesPAdj(thoracic.3.14.de.genes.proportion.ruvg)
WriteDataFrameAsTsv(thoracic.3.14.de.genes.proportion.ruvg.sign.005, file.name.path = "/media/dario/dati/time_course/Thoracic/results_3_14/de_genes/thoracic_3_14_filtered_counts__Proportion__normalized__ruvg_de_genes_LRT_significant_005.tsv")
# 
# thoracic.3.14.de.genes.wilcoxon.uqua <- MainFunction(counts.dataframe = thoracic.3.14.counts, design.matrix = thoracic.3.14.design.matrix, output.results.folder = thoracic.output.folder.3.14, prefix.output.file = "thoracic_3_14", prefix.plot.label = "thoracic_3_14", is.time.course = TRUE, filter.method = "Wilcoxon")
# thoracic.3.14.de.genes.wilcoxon.uqua.sign.005 <- SignificantDeGenesPAdj(thoracic.3.14.de.genes.wilcoxon.uqua)
