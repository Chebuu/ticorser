## cervical first 3 time points

cervical.complete.counts <- ReadDataFrameFromTsv("Cervical/results/counts/FeatureCounts/cervical_ensgtf_reordered_counts.tsv")
colnames(cervical.complete.counts)
cervical.3.14.counts <- cervical.complete.counts[, c(1:24)]
colnames(cervical.3.14.counts)

cervical.design.file <- file.path("/media/dario/dati/time_course/Cervical/design_file/cervical_design_file.txt")
cervical.design.matrix <- ReadDataFrameFromTsv(cervical.design.file)
cervical.3.14.design.matrix <- cervical.design.matrix[-which(cervical.design.matrix$Times %in% "56d"),]

cervical.output.folder.3.14 <- file.path("/media/dario/dati/time_course/Cervical", "results_3_14")

cervical.3.14.de.genes.proportion.uqua <- MainFunction(counts.dataframe = cervical.3.14.counts, design.matrix = cervical.3.14.design.matrix, output.results.folder = cervical.output.folder.3.14, prefix.output.file = "cervical_3_14", prefix.plot.label = "cervical_3_14", is.time.course = TRUE, filter.method = "Proportion")

cervical.3.14.de.genes.proportion.uqua.sign.005 <- SignificantDeGenesPAdj(cervical.3.14.de.genes.proportion.uqua)

WriteDataFrameAsTsv(cervical.3.14.de.genes.proportion.uqua.sign.005, file.name.path = "/media/dario/dati/time_course/Cervical/results_3_14/de_genes/cervical_3_14_filtered_counts__Proportion__normalized__uqua_de_genes_LRT_significant_005.tsv")

cervical.3.14.proportion.estimated.neg.genes <- EstimateNegativeControlGenesForRUV(de.genes = cervical.3.14.de.genes.proportion.uqua, n.tail.genes = 2000, counts.dataset = cervical.3.14.counts, n.genes.per.hist.break = 1)

cervical.3.14.de.genes.proportion.ruvg <- MainFunction(counts.dataframe = cervical.3.14.counts, design.matrix = cervical.3.14.design.matrix, output.results.folder = cervical.output.folder.3.14, prefix.output.file = "cervical_3_14", prefix.plot.label = "cervical_3_14", is.time.course = TRUE, filter.method = "Proportion", normalization.method = "ruvg", estimated.genes = cervical.3.14.proportion.estimated.neg.genes)
cervical.3.14.de.genes.proportion.ruvg.sign.005 <- SignificantDeGenesPAdj(cervical.3.14.de.genes.proportion.ruvg)
WriteDataFrameAsTsv(cervical.3.14.de.genes.proportion.ruvg.sign.005, file.name.path = "/media/dario/dati/time_course/Cervical/results_3_14/de_genes/cervical_3_14_filtered_counts__Proportion__normalized__ruvg_de_genes_LRT_significant_005.tsv")
# 
# cervical.3.14.de.genes.wilcoxon.uqua <- MainFunction(counts.dataframe = cervical.3.14.counts, design.matrix = cervical.3.14.design.matrix, output.results.folder = cervical.output.folder.3.14, prefix.output.file = "cervical_3_14", prefix.plot.label = "cervical_3_14", is.time.course = TRUE, filter.method = "Wilcoxon")
# cervical.3.14.de.genes.wilcoxon.uqua.sign.005 <- SignificantDeGenesPAdj(cervical.3.14.de.genes.wilcoxon.uqua)
