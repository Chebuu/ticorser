## ### all conditions 
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/includes.R")

cerv.counts <- ReadDataFrameFromTsv("Cervical/Cervical_results/counts/FeatureCounts/cervical_ensgtf_reordered_counts.tsv")
thor.counts <- ReadDataFrameFromTsv("Thoracic/Thoracic_results/counts/FeatureCounts/thoracic_ensgtf_reordered_counts.tsv")
cerv.thor.counts <- cbind(cerv.counts, thor.counts)

all.design.file <- "cervical_thoracic/design_file/cervical_thoracic_all_design_file.txt"
cer.thor.des.mat <- ReadDataFrameFromTsv(all.design.file)
output.results.folder <- "cervical_thoracic_whole/results_final"

# cerv.thor.proj <<- Project$new(name="cervical_thoracic_whole", main.working.path=getwd())
# 
# working.project <<- cerv.thor.proj


## wilcoxon

# 
# wilc.whole.data.frames.fqua <- FilterAndNormalizeCounts(counts.dataframe = cerv.thor.counts, design.matrix = cer.thor.des.mat, output.results.folder = output.results.folder
#                                                    , prefix.output.file = "cervical_thoracic_whole", prefix.plot.label = "cerv_thor_all"
#                                                    , filter.method = "Wilcoxon", normalization.method = "fqua", ellipse.in.pca = "both")
# 
# wilc.whole.data.frames.uqua <- FilterAndNormalizeCounts(counts.dataframe = cerv.thor.counts, design.matrix = cer.thor.des.mat, output.results.folder = output.results.folder
#                                                    , prefix.output.file = "cervical_thoracic_whole", prefix.plot.label = "cerv_thor_all"
#                                                    , filter.method = "Wilcoxon", normalization.method = "uqua", ellipse.in.pca = "both")

## proportion
whole.data.frames.fqua <- FilterAndNormalizeCounts(counts.dataframe = cerv.thor.counts, design.matrix = cer.thor.des.mat, output.results.folder = output.results.folder
                                                   , prefix.output.file = "cervical_thoracic_whole", prefix.plot.label = "cerv_thor_all"
                                                   , filter.method = "Proportion", normalization.method = "fqua", ellipse.in.pca = "both")

whole.data.frames.uqua <- FilterAndNormalizeCounts(counts.dataframe = cerv.thor.counts, design.matrix = cer.thor.des.mat, output.results.folder = output.results.folder
                                                   , prefix.output.file = "cervical_thoracic_whole", prefix.plot.label = "cerv_thor_all"
                                                   , filter.method = "Proportion", normalization.method = "uqua", ellipse.in.pca = "both")

normalized.counts.prop.fqua <- whole.data.frames.fqua$normalized.counts.df
normalized.counts.prop.uqua <- whole.data.frames.uqua$normalized.counts.df
# normalized.counts.prop.uqua <- ReadDataFrameFromTsv("cervical_thoracic_whole/results_final/counts/FeatureCounts/Proportion/uqua/cervical_thoracic_whole_Proportion_uqua.tsv")

cerv.thor.prop.counts <- ReadDataFrameFromTsv("cervical_thoracic_whole/results_final/counts/FeatureCounts/Proportion/cervical_thoracic_whole_Proportion.tsv")

## cerv
cervical.design.matrix <- ReadDataFrameFromTsv(file.name.path = "Cervical/design_file/cervical_design_file.txt")

# output.counts.folder <- UpdateFolderPath(output.results.folder, "4Points_Cervical/counts/FeatureCounts/Proportion/fqua")
cerv.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Cervical")
# # plot.folder <- UpdateFolderPath(output.results.folder, "4Points_Cervical") 
# 
# cerv.de.fqua.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.fqua), design.matrix = cervical.design.matrix, 
#                                       de.test = "DeSeqTime", results.folder =  cerv.output.results.folder,
#                                       prefix.label = "4Points_Cervical prop fqua",
#                                       normalize.data.flag = TRUE, normalization.method = "uqua")
# 
#  
# 
# cerv.de.uqua.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
#                                       de.test = "DeSeqTime", results.folder =  cerv.output.results.folder,
#                                       prefix.label = "4Points_Cervical prop uqua",
#                                       normalize.data.flag = TRUE, normalization.method = "uqua")
# 
# cerv.de.uqua.fqua.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
#                                            de.test = "DeSeqTime", results.folder =  cerv.output.results.folder,
#                                            prefix.label = "4Points_Cervical prop uqua",
#                                            normalize.data.flag = TRUE, normalization.method = "fqua")
# 

cerv.de.uqua.notnorm2.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
                                                    de.test = "DeSeqTime", results.folder =  cerv.output.results.folder,
                                                    prefix.label = "4Points_Cervical prop uqua",
                                                    normalize.data.flag = FALSE)

##analisi cervical su matrice dei counts di input non normalizzata filtrata
# 
# cerv.de.unorm.uqua <- PerformDEAnalysis(whole.counts = cerv.thor.prop.counts, design.matrix = cervical.design.matrix, 
#                                         de.test = "DeSeqTime", results.folder =  cerv.output.results.folder,
#                                         prefix.label = "4Points_Cervical prop",
#                                         normalize.data.flag = TRUE, normalization.method = "uqua")
# cerv.prop.counts <- cerv.thor.prop.counts[, which(colnames(cerv.thor.prop.counts)%in%rownames(cervical.design.matrix))]
# cerv.uqua.est.neg.genes <- EstimateNegativeControlGenesForRUV(de.genes = cerv.de.unorm.uqua$de.not.na, n.tail.genes = 2000, counts.dataset = cerv.prop.counts, n.genes.per.hist.break = 1)
# 
# 
# cerv.de.unorm.ruvg <- PerformDEAnalysis(whole.counts = cerv.thor.prop.counts, design.matrix = cervical.design.matrix, 
#                                         de.test = "DeSeqTime", results.folder =  cerv.output.results.folder,
#                                         prefix.label = "4Points_Cervical prop",
#                                         normalize.data.flag = TRUE, normalization.method = "ruvg", negative.genes.list = cerv.uqua.est.neg.genes)
# 
# venn.cervical <- Venn3de(x = rownames(cerv.de.uqua.notnorm2.des4time$SIGN), y = rownames(cerv.de.unorm.uqua$SIGN), z = rownames(cerv.de.unorm.ruvg$SIGN), label1 = "cerv-whole-uqua", label2 = "cerv-singl-uqua", label3 = "cerv-singl-ruvg", title = "cervical de venn", intersection.flag = FALSE, intersection.exclusion.flag = FALSE, plot.dir = "./")

## thor
thor.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Thoracic")
thoracic.design.matrix <- ReadDataFrameFromTsv(file.name.path = "Thoracic/design_file/thoracic_design_file.txt")
# 
# thor.de.fqua.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.fqua), design.matrix = thoracic.design.matrix, 
#                                       de.test = "DeSeqTime", results.folder =  thor.output.results.folder,
#                                       prefix.label = "4Points_Thoracic prop fqua",
#                                       normalize.data.flag = TRUE, normalization.method = "uqua")
# 
# thor.de.uqua.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = thoracic.design.matrix, 
#                                            de.test = "DeSeqTime", results.folder =  thor.output.results.folder,
#                                            prefix.label = "4Points_Thoracic prop uqua",
#                                            normalize.data.flag = TRUE, normalization.method = "uqua")
# 
# 
# thor.de.uqua.fqua.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = thoracic.design.matrix, 
#                                            de.test = "DeSeqTime", results.folder =  thor.output.results.folder,
#                                            prefix.label = "4Points_Thoracic prop uqua",
#                                            normalize.data.flag = TRUE, normalization.method = "fqua")

thor.de.uqua.notnorm2.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = thoracic.design.matrix, 
                                                    de.test = "DeSeqTime", results.folder =  thor.output.results.folder,
                                                    prefix.label = "4Points_Thoracic prop uqua",
                                                    normalize.data.flag = FALSE)
# ### analisi thor su matrice whole non normalizzata
# thor.de.unorm.uqua <- PerformDEAnalysis(whole.counts = cerv.thor.prop.counts, design.matrix = thoracic.design.matrix, 
#                                                     de.test = "DeSeqTime", results.folder =  thor.output.results.folder,
#                                                     prefix.label = "4Points_Thoracic prop",
#                                                     normalize.data.flag = TRUE, normalization.method = "uqua")
# 
# thor.prop.counts <- cerv.thor.prop.counts[, which(colnames(cerv.thor.prop.counts)%in%rownames(thoracic.design.matrix))]
# 
# thor.estimated.genes.list <- EstimateNegativeControlGenesForRUV(de.genes = thor.de.unorm.uqua$de.not.na, n.tail.genes = 2000, counts.dataset = thor.prop.counts, n.genes.per.hist.break = 1)
# 
# thor.de.unorm.ruvg <- PerformDEAnalysis(whole.counts = cerv.thor.prop.counts, design.matrix = thoracic.design.matrix, 
#                                         de.test = "DeSeqTime", results.folder =  thor.output.results.folder,
#                                         prefix.label = "4Points_Thoracic prop",
#                                         normalize.data.flag = TRUE, normalization.method = "ruvg", negative.genes.list = thor.estimated.genes.list)
# 
# venn.thoracic <- Venn3de(x = rownames(thor.de.uqua.notnorm2.des4time$SIGN), y = rownames(thor.de.unorm.uqua$SIGN), z = rownames(thor.de.unorm.ruvg$SIGN), label1 = "thor-whole-uqua", label2 = "thor-singl-uqua", label3 = "thor-singl-ruvg", title = "Thoracic de venn", intersection.flag = FALSE, intersection.exclusion.flag = FALSE, plot.dir = "./")

## cerv-thor
cerv.thor.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic/design_file/cervical_thoracic_design_file.txt")
cerv.thor.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Cerv-Thor")

# 
# cerv.thor.de.fqua.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.fqua), design.matrix = cerv.thor.design.matrix, 
#                                            de.test = "DeSeqTime", results.folder =  cerv.thor.output.results.folder,
#                                            prefix.label = "4Points_Cerv_Thor prop fqua",
#                                            normalize.data.flag = TRUE, normalization.method = "uqua")
# 
# 
# 
# cerv.thor.de.uqua.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.design.matrix, 
#                                            de.test = "DeSeqTime", results.folder =  cerv.thor.output.results.folder,
#                                            prefix.label = "4Points_Cerv_Thor prop uqua",
#                                            normalize.data.flag = TRUE, normalization.method = "uqua")
# 
# cerv.thor.de.uqua.fqua.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.design.matrix, 
#                                                 de.test = "DeSeqTime", results.folder =  cerv.thor.output.results.folder,
#                                                 prefix.label = "4Points_Cerv_Thor prop uqua",
#                                                 normalize.data.flag = TRUE, normalization.method = "fqua")


cerv.thor.de.uqua.notnorm2.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.design.matrix, 
                                                         de.test = "DeSeqTime", results.folder =  cerv.thor.output.results.folder,
                                                         prefix.label = "4Points_Cerv_Thor prop uqua",
                                                         normalize.data.flag = FALSE)

# ### analisi thor su matrice whole non normalizzata
# cerv.thor.t.prop.counts <- cerv.thor.prop.counts[, which(colnames(cerv.thor.prop.counts)%in%rownames(cerv.thor.design.matrix))]
# 
# cerv.thor.de.unorm.uqua <- PerformDEAnalysis(whole.counts = cerv.thor.prop.counts, design.matrix = cerv.thor.design.matrix, 
#                                         de.test = "DeSeqTime", results.folder =  cerv.thor.output.results.folder,
#                                         prefix.label = "4Points_Thoracic prop",
#                                         normalize.data.flag = TRUE, normalization.method = "uqua")
# 
# cerv.thor.estimated.genes.list <- EstimateNegativeControlGenesForRUV(de.genes = cerv.thor.de.unorm.uqua$de.not.na, n.tail.genes = 2000, counts.dataset = cerv.thor.t.prop.counts, n.genes.per.hist.break = 1)
# 
# 
# cerv.thor.de.unorm.ruvg <- PerformDEAnalysis(whole.counts = cerv.thor.prop.counts, design.matrix = cerv.thor.design.matrix, 
#                                                          de.test = "DeSeqTime", results.folder =  cerv.thor.output.results.folder,
#                                                          prefix.label = "4Points_Cerv_Thor prop",
#                                                          normalize.data.flag = TRUE, normalization.method = "ruvg", negative.genes.list = cerv.thor.estimated.genes.list)
# 
# venn.thoracic <- Venn3de(x = rownames(cerv.thor.de.uqua.notnorm2.des4time$SIGN), y = rownames(cerv.thor.de.unorm.uqua$SIGN), z = rownames(cerv.thor.de.unorm.ruvg$SIGN), label1 = "cerv-thor-w-uqua", label2 = "cerv-thor-s-uqua", label3 = "cerv-thor-s-ruvg", title = "cerv-thor de venn", intersection.flag = FALSE, intersection.exclusion.flag = FALSE, plot.dir = "./")

# ############# comparing data
# cerv.ruved <- ReadDataFrameFromTsv("cervical_thoracic_whole/results_final/4Points_Cervical/counts/FeatureCounts/ruvg/4Points_Cervical_prop_ruvg.tsv")
# thor.ruved <- ReadDataFrameFromTsv("cervical_thoracic_whole/results_final/4Points_Thoracic/counts/FeatureCounts/ruvg/4Points_Thoracic_prop_ruvg.tsv")
# 
# cerv.thor.ruved <- cbind(cerv.ruved, thor.ruved)
# 
# cerv.thor.ruved.design <- rbind(cervical.design.matrix, thoracic.design.matrix)
# 
# PlotTimesBoxplot(data.frame.to.plot = cerv.thor.ruved, design.matrix = cerv.thor.design.matrix, output.path = "./", prefix.plot = "cerv-thor ruved", save.plot = TRUE,show.plot.flag = FALSE, plotly.flag = FALSE)
# 

