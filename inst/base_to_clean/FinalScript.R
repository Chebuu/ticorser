
# capture all the output to a file.
zz <- file("console1.Rout", open = "wt")
sink(zz)
sink(zz, type = "message")
# try(log("a"))

## ### all conditions 
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/includes.R")

rat5.ensembl.downloaded.gtf <- "/media/dario/dati/time_course/downloaded_references/ensembl_site/Rattus_norvegicus.Rnor_5.0.79_gene_biotype_protein_coding_canonical_chrs_NO_MT_withHeader.gtf"
output.results.folder <- "cervical_thoracic_whole/results_untreated"

all.bam.folder <- "/media/dario/dati/time_course/cervical_thoracic/all_bams/"
cervical.thoracic.whole.bams <- list.files(all.bam.folder, pattern="bam$",full.names=TRUE)

cervical.thoracic.whole.counts.ens.gtf <- CountBamFilesFeatureCounts(bam.files = cervical.thoracic.whole.bams, annotation.file = rat5.ensembl.downloaded.gtf, output.folder = output.results.folder, prefix.output.file="whole_cervical_thoracic")

cervical.thoracic.whole.renamed.counts <- cervical.thoracic.whole.counts.ens.gtf
colnames(cervical.thoracic.whole.renamed.counts) <- c( "cervical_03d_t1", "cervical_03d_t2", "cervical_03d_t3", "cervical_03d_t4", "cervical_03d_t5",
                                        "cervical_03d_u1", "cervical_03d_u2", "cervical_03d_u3",
                                        
                                        "cervical_07d_t1", "cervical_07d_t2", "cervical_07d_t3", "cervical_07d_t4", "cervical_07d_t5",
                                        "cervical_07d_u1", "cervical_07d_u2", "cervical_07d_u3",
                                        
                                        "cervical_14d_t1", "cervical_14d_t2", "cervical_14d_t3", "cervical_14d_t4", "cervical_14d_t5",
                                        "cervical_14d_u1", "cervical_14d_u2", "cervical_14d_u3",
                                        
                                        "cervical_56d_t1", "cervical_56d_t2", "cervical_56d_t3", "cervical_56d_t4", "cervical_56d_t5",
                                        "cervical_56d_u1", "cervical_56d_u2", "cervical_56d_u3",
                                        
                                        "thoracic_03d_t1", "thoracic_03d_t2", "thoracic_03d_t3", "thoracic_03d_t4", "thoracic_03d_t5",
                                        "thoracic_03d_u1", "thoracic_03d_u2", "thoracic_03d_u3", 
                                        
                                        "thoracic_07d_t1", "thoracic_07d_t2", "thoracic_07d_t3", "thoracic_07d_t4", "thoracic_07d_t5",
                                        "thoracic_07d_u1", "thoracic_07d_u2", "thoracic_07d_u3",
                                        
                                        "thoracic_14d_t1", "thoracic_14d_t2", "thoracic_14d_t3", "thoracic_14d_t4", "thoracic_14d_t5",
                                        "thoracic_14d_u1", "thoracic_14d_u2", "thoracic_14d_u3", 
                                        
                                        "thoracic_56d_t1", "thoracic_56d_t2", "thoracic_56d_t3", "thoracic_56d_t4", "thoracic_56d_t5",
                                        "thoracic_56d_u1"
                                        )
                                        
WriteDataFrameAsTsv(data.frame.to.save = cervical.thoracic.whole.renamed.counts, file.name.path = "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_8/counts/FeatureCounts/whole_cervical_thoracic_Counts_FeatureCounts_renamed.txt", col.names = NA)
# cerv.counts <- ReadDataFrameFromTsv("Cervical/Cervical_results/counts/FeatureCounts/cervical_ensgtf_reordered_counts.tsv")
# thor.counts <- ReadDataFrameFromTsv("Thoracic/Thoracic_results/counts/FeatureCounts/thoracic_ensgtf_reordered_counts.tsv")
# cerv.thor.counts <- cbind(cerv.counts, thor.counts)

all.design.file <- "/media/dario/dati/time_course/cervical_thoracic/design_file/cervical_thoracic_all_design_file.txt"
cer.thor.des.mat <- ReadDataFrameFromTsv(all.design.file)


design.files.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/design_file"
# normalized.counts.prop.uqua<-ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/james_results/counts/FeatureCounts/Proportion/uqua/cervical_thoracic_whole_Proportion_uqua.tsv")

# cervical.thoracic.whole.renamed.counts <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_6/counts/FeatureCounts/whole_cervical_thoracic_Counts_FeatureCounts_renamed.txt.tsv")

## proportion
whole.data.frames.uqua <- FilterAndNormalizeCounts(counts.dataframe = cervical.thoracic.whole.renamed.counts, design.matrix = cer.thor.des.mat, output.results.folder = output.results.folder
                                              , prefix.output.file = "cervical_thoracic_whole", prefix.plot.label = "cerv_thor_all"
                                              , filter.method = "Proportion", normalization.method = "uqua", ellipse.in.pca = "both")

# counts.dataframe = cervical.thoracic.whole.renamed.counts; design.matrix = cer.thor.des.mat; output.results.folder = output.results.folder; prefix.output.file = "cervical_thoracic_whole"; prefix.plot.label = "cerv_thor_all"; filter.method = "Proportion"; normalization.method = "uqua"; ellipse.in.pca = "both"
# normalized.counts.prop.uqua<-ReadDataFrameFromTsv("cervical_thoracic_whole/results_final_10/counts/FeatureCounts/Proportion/uqua/cervical_thoracic_whole_Proportion_uqua.tsv")
normalized.counts.prop.uqua <- whole.data.frames.uqua$normalized.counts.df


ens.symb.biomart.map <- ReadDataFrameFromTsv("downloaded_references/biomart/rnor5.0_biomart_ensembls_associatedgenenames.txt", row.names.col = NULL)
head(ens.symb.biomart.map)

# 
# ############# plotting molecules across time
# genes <- c( "Lama1", "Lama2", "Lama4", "Lama5", "Lamb1", "Lamc1", "Cspg4", "Cspg5", "Acan", "Bcan", "Ncan", "Vcan") 
# 
# map.ens.gs <- ReadDataFrameFromTsv("/home/dario/Scrivania/Dropbox/Lavori/IAC/coding/time_course/conversion_maps/rnor5.0_biomart_ensembls_associatedgenenames.txt", row.names.col = NULL)
# 
# normalized.counts.ordered <- normalized.counts.prop.uqua[order(rownames(normalized.counts.prop.uqua)),]
# map.ens.gs.ord <- map.ens.gs[order(map.ens.gs$Ensembl.Gene.ID),]
# gene.names <- map.ens.gs.ord$Associated.Gene.Name[which( map.ens.gs.ord$Ensembl.Gene.ID %in% rownames(normalized.counts.ordered))]
# normalized.counts.ordered.gs <- normalized.counts.ordered
# normalized.counts.ordered.gs$gene.names <- gene.names
# normalized.counts.ordered.gs$gene.names[which(tolower(normalized.counts.ordered.gs$gene.names) %in% tolower(genes))]
# time.plot.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_1/plots/Proportion/uqua"
# 
# for (gene in genes) {
#   PlotCountsAlongTimes(normalized.counts = normalized.counts.ordered.gs, design.matrix = cervical.design.matrix, gene.name = gene, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = time.plot.folder, prefix.plot = "Cervical")
#   PlotCountsAlongTimes(normalized.counts = normalized.counts.ordered.gs, design.matrix = thoracic.design.matrix, gene.name = gene, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = time.plot.folder, prefix.plot = "Thoracic")
#   PlotCountsAlongTimes(normalized.counts = normalized.counts.ordered.gs, design.matrix = cerv.thor.design.matrix, gene.name = gene, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = time.plot.folder, prefix.plot = "Cerv-Thor")
# }

# normalized.counts = normalized.counts.ordered.gs; design.matrix = thoracic.design.matrix; gene.name = gene; show.plot.flag = FALSE; plotly.flag = TRUE; save.plot = TRUE; plot.folder = time.plot.folder; prefix.plot = "Thoracic"


# normalized.counts.prop.uqua <- ReadDataFrameFromTsv("cervical_thoracic_whole/results_final/counts/FeatureCounts/Proportion/uqua/cervical_thoracic_whole_Proportion_uqua.tsv")

# cerv.thor.prop.counts <- ReadDataFrameFromTsv("cervical_thoracic_whole/results_final/counts/FeatureCounts/Proportion/cervical_thoracic_whole_Proportion.tsv")

# ## cerv 4times
# cervical.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/design_file/4timepoints/cervical_design_file.txt")
# 
# cerv.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Cervical")
# 
# 
# cerv.de.uqua.notnorm2.des4time.final.formula3 <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
#                                                                    de.test = "DeSeqTime", results.folder =  cerv.output.results.folder,
#                                                                    prefix.label = "4Points_Cervical prop uqua",
#                                                                    normalize.data.flag = FALSE, threshold = 0.05, conversion.map=ens.symb.biomart.map)
# 
# cerv.de.uqua.notnorm2.des4time.final.formula2 <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
#                                            de.test = "DeSeqTime", results.folder =  cerv.output.results.folder,
#                                            prefix.label = "4Points_Cervical prop uqua",
#                                            normalize.data.flag = FALSE, threshold = 0.05, conversion.map=ens.symb.biomart.map)
# 
# cerv.de.uqua.notnorm2.des4time.mid.formula <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
#                                                     de.test = "DeSeqTime", results.folder =  cerv.output.results.folder,
#                                                     prefix.label = "4Points_Cervical prop uqua",
#                                                     normalize.data.flag = FALSE, threshold = 0.05, conversion.map=ens.symb.biomart.map)
# 
# cerv.de.uqua.notnorm2.des4time.first.formula <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
#                                                     de.test = "DeSeqTime", results.folder =  cerv.output.results.folder,
#                                                     prefix.label = "4Points_Cervical prop uqua",
#                                                     normalize.data.flag = FALSE, threshold = 0.05, conversion.map=ens.symb.biomart.map)
# 
# # whole.counts = round(normalized.counts.prop.uqua); design.matrix = cervical.design.matrix; de.test = "DeSeqTime"; results.folder =  cerv.output.results.folder;  prefix.label = "4Points_Cervical prop uqua";  normalize.data.flag = FALSE; threshold = 0.05
# 
# ## thor 4times
# thor.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Thoracic")
# thoracic.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/design_file/4timepoints/thoracic_design_file.txt")
# 
# thor.de.uqua.notnorm2.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = thoracic.design.matrix, 
#                                                 de.test = "DeSeqTime", results.folder =  thor.output.results.folder,
#                                                 prefix.label = "4Points_Thoracic prop uqua",
#                                                 normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map)
# 
# ## cerv-thor 4times
# cerv.thor.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/design_file/4timepoints/cervical_thoracic_design_file.txt")
# cerv.thor.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Cerv-Thor")
# 
# 
# cerv.thor.de.uqua.notnorm2.des4time <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.design.matrix, 
#                                                      de.test = "DeSeqTime", results.folder =  cerv.thor.output.results.folder,
#                                                      prefix.label = "4Points_Cerv_Thor prop uqua",
#                                                      normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map)
# ############## convert data



# ############## comparing data 4 times
# venn.output.plot.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_1/plots/Proportion/uqua"
# 
# Venn3de(x = rownames(cerv.de.uqua.notnorm2.des4time$SIGN), y = rownames(thor.de.uqua.notnorm2.des4time$SIGN), z = rownames(cerv.thor.de.uqua.notnorm2.des4time$SIGN), 
#         label1 = "cervical", label2 = "thoracic", label3 = "cerv-thor", title = "DE significant genes", 
#         intersection.flag = TRUE, intersection.exclusion.flag = TRUE, plot.dir = venn.output.plot.folder)
# 
# Venn3de(x = rownames(cerv.de.uqua.notnorm2.des4time$UP), y = rownames(thor.de.uqua.notnorm2.des4time$UP), z = rownames(cerv.thor.de.uqua.notnorm2.des4time$UP), 
#         label1 = "cervical_up", label2 = "thoracic_up", label3 = "cerv-thor_up", title = "DE significant UP genes", 
#         intersection.flag = TRUE, intersection.exclusion.flag = TRUE, plot.dir = venn.output.plot.folder)
# 
# 
# Venn3de(x = rownames(cerv.de.uqua.notnorm2.des4time$DOWN), y = rownames(thor.de.uqua.notnorm2.des4time$DOWN), z = rownames(cerv.thor.de.uqua.notnorm2.des4time$DOWN), 
#         label1 = "cervical_down", label2 = "thoracic_down", label3 = "cerv-thor_down", title = "DE significant DOWN genes", 
#         intersection.flag = TRUE, intersection.exclusion.flag = TRUE, plot.dir = venn.output.plot.folder)

### subacute

##normalized megamatrix
# normalized.counts.prop.uqua <- whole.data.frames.uqua$normalized.counts.df

# ## cerv 3times
# subacute.design.folder <- file.path(design.files.folder, "subacute")
# 
# subacute.cervical.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(subacute.design.folder, "cervical_subacute_design_file.csv"))
# 
# cerv.subacute.output.results.folder <- UpdateFolderPath(output.results.folder, "subacute_Cervical")
# 
# cervical.de.uqua.notnorm2.subacute <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = subacute.cervical.design.matrix, 
#                                                          de.test = "DeSeqTime", results.folder =  cerv.subacute.output.results.folder,
#                                                          prefix.label = "Cervical subacute prop uqua",
#                                                          normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map)
# 
# ## thoracic subacute
# subacute.thoracic.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(subacute.design.folder, "thoracic_subacute_design_file.csv"))
# 
# thor.subacute.output.results.folder <- UpdateFolderPath(output.results.folder, "subacute_Thoracic")
# 
# thoracic.de.uqua.notnorm2.subacute <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = subacute.thoracic.design.matrix, 
#                                                         de.test = "DeSeqTime", results.folder =  thor.subacute.output.results.folder,
#                                                         prefix.label = "Thoracic subacute prop uqua",
#                                                         normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map)
# 
# ## cerv-thor subacute
# subacute.cerv.thor.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(subacute.design.folder, "cerv-thor_subacute_design_file.csv"))
# 
# cerv.thor.subacute.output.results.folder <- UpdateFolderPath(output.results.folder, "subacute_Cerv-Thor")
# 
# cerv.thor.de.uqua.notnorm2.subacute <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = subacute.cerv.thor.design.matrix, 
#                                                         de.test = "DeSeqTime", results.folder =  cerv.thor.subacute.output.results.folder,
#                                                         prefix.label = "Cerv-Thor subacute prop uqua",
#                                                         normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map)

# ### comparing subacute
# venn.output.plot.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_1/plots/Proportion/uqua/"
# 
# Venn3de(x = rownames(cervical.de.uqua.notnorm2.subacute$SIGN), y = rownames(thoracic.de.uqua.notnorm2.subacute$SIGN), z = rownames(cerv.thor.de.uqua.notnorm2.subacute$SIGN), 
#         label1 = "cervical_subac", label2 = "thoracic_subac", label3 = "cerv-thor_subac", title = "DE Subacute significant genes", 
#         intersection.flag = TRUE, intersection.exclusion.flag = TRUE, plot.dir = venn.output.plot.folder)
# 
# Venn3de(x = rownames(cervical.de.uqua.notnorm2.subacute$UP), y = rownames(thoracic.de.uqua.notnorm2.subacute$UP), z = rownames(cerv.thor.de.uqua.notnorm2.subacute$UP), 
#         label1 = "cervical_up_subac", label2 = "thoracic_up_subac", label3 = "cerv-thor_up_subac", title = "DE Subacute significant UP genes", 
#         intersection.flag = TRUE, intersection.exclusion.flag = TRUE, plot.dir = venn.output.plot.folder)
# 
# 
# Venn3de(x = rownames(cervical.de.uqua.notnorm2.subacute$DOWN), y = rownames(thoracic.de.uqua.notnorm2.subacute$DOWN), z = rownames(cerv.thor.de.uqua.notnorm2.subacute$DOWN), 
#         label1 = "cervical_down_subac", label2 = "thoracic_down_subac", label3 = "cerv-thor_down_subac", title = "DE Subacute significant DOWN genes", 
#         intersection.flag = TRUE, intersection.exclusion.flag = TRUE, plot.dir = venn.output.plot.folder)

####### single timepoint 3 days
design.folder.03d <- file.path(design.files.folder, "time_by_time/03d")

## cerv
cervical.03d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.03d, "cervical_03d_design_file.csv"))

cerv.03d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "03d", "Cervical_deseq")

cervical.de.uqua.notnorm2.03d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.03d.design.matrix, 
                                                        de.test = "DeSeq", results.folder =  cerv.03d.output.results.folder,
                                                        prefix.label = "Cervical 03d prop uqua",
                                                        normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)
# whole.counts = round(normalized.counts.prop.uqua); design.matrix = cervical.03d.design.matrix; de.test = "DeSeq"; results.folder =  cerv.03d.output.results.folder;prefix.label = "Cervical 03d prop uqua";normalize.data.flag = FALSE; conversion.map=ens.symb.biomart.map
# 
# cerv.03d.output.results.folder95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "03d", "Cervical_noiseq_095")
# 
# cervical.de.uqua.notnorm2.03d.noiseqbio <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cervical.03d.design.matrix, 
#                                                         de.test = "NOISeqBio", results.folder =  cerv.03d.output.results.folder95,
#                                                         prefix.label = "Cervical 03d prop uqua",
#                                                         normalize.data.flag = FALSE, threshold=0.95, conversion.map=ens.symb.biomart.map)

cerv.03d.output.results.folder99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "03d", "Cervical_noiseq_099")

cervical.de.uqua.notnorm2.03d.noiseqbio99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cervical.03d.design.matrix, 
                                                             de.test = "NOISeqBio", results.folder =  cerv.03d.output.results.folder99,
                                                             prefix.label = "Cervical 03d prop uqua",
                                                             normalize.data.flag = FALSE, threshold=0.99, conversion.map=ens.symb.biomart.map)


## thor 
thoracic.03d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.03d, "thoracic_03d_design_file.csv"))

thor.03d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "03d", "Thoracic_deseq")

thoracic.de.uqua.notnorm2.03d <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = thoracic.03d.design.matrix, 
                                                   de.test = "DeSeq", results.folder =  thor.03d.output.results.folder,
                                                   prefix.label = "Thoracic 03d prop uqua",
                                                   normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)
# 
# thor.03d.output.results.folder.noiseq95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "03d", "Thoracic_noiseq_095")
# 
# thoracic.de.uqua.notnorm2.03d.noiseq95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = thoracic.03d.design.matrix, 
#                                                    de.test = "NOISeqBio", results.folder =  thor.03d.output.results.folder.noiseq95,
#                                                    prefix.label = "Thoracic 03d prop uqua",
#                                                    normalize.data.flag = FALSE, threshold = 0.95, conversion.map=ens.symb.biomart.map)

# 
# thor.03d.output.results.folder.noiseq99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "03d", "Thoracic_noiseq_099")
# 
# thoracic.de.uqua.notnorm2.03d.noiseq99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = thoracic.03d.design.matrix, 
#                                                             de.test = "NOISeqBio", results.folder =  thor.03d.output.results.folder.noiseq99,
#                                                             prefix.label = "Thoracic 03d prop uqua",
#                                                             normalize.data.flag = FALSE, threshold = 0.99, conversion.map=ens.symb.biomart.map)


## cerv-thor 
cerv.thor.03d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.03d, "cerv-thor_03d_design_file.csv"))

cerv.thor.03d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "03d", "Cerv-Thor_deseq")

cerv.thor.de.uqua.notnorm2.03d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.03d.design.matrix, 
                                                   de.test = "DeSeq", results.folder =  cerv.thor.03d.output.results.folder,
                                                   prefix.label = "Cerv-Thor 03d prop uqua",
                                                   normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)
## nois
# cerv.thor.03d.output.results.folder.noi95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "03d", "Cerv-Thor_noiseq_95")
# 
# cerv.thor.de.uqua.notnorm2.03d.noi95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cerv.thor.03d.design.matrix,
#                                                     de.test = "NOISeqBio", results.folder =  cerv.thor.03d.output.results.folder.noi95,
#                                                     prefix.label = "Cerv-Thor 03d prop uqua",
#                                                     normalize.data.flag = FALSE, threshold = 0.95, conversion.map=ens.symb.biomart.map)
# 
# cerv.thor.03d.output.results.folder.noi99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "03d", "Cerv-Thor_noiseq_99")
# 
# cerv.thor.de.uqua.notnorm2.03d.noi99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cerv.thor.03d.design.matrix, 
#                                                     de.test = "NOISeqBio", results.folder =  cerv.thor.03d.output.results.folder.noi99,
#                                                     prefix.label = "Cerv-Thor 03d prop uqua",
#                                                     normalize.data.flag = FALSE, threshold = 0.99, conversion.map=ens.symb.biomart.map)


## cerv-thor untr 

cerv.thor.03d.design.matrix.untr <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.03d, "cervical_thoracic_03d_untr_design_file.txt"))

cerv.thor.03d.output.results.folder.untr <- UpdateFolderPath(output.results.folder, "Time_by_Time", "03d", "Cerv-Thor_untr_deseq")

cerv.thor.de.uqua.notnorm2.03d.deseq.untr <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.03d.design.matrix.untr, 
                                                          de.test = "DeSeq", results.folder =  cerv.thor.03d.output.results.folder.untr,
                                                          prefix.label = "Cerv-Thor untr 03d prop uqua",
                                                          normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)
# ## comparing 03d
# venn.output.plot.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_1/plots/Proportion/uqua/"
# 
# VennDE3ListsSuited(cervical.de.uqua.notnorm2.03d, thoracic.de.uqua.notnorm2.03d, cerv.thor.de.uqua.notnorm2.03d, title = "03 days", venn.plot.dir = venn.output.plot.folder)
# 


############### single point 07d
design.folder.07d <- file.path(design.files.folder, "time_by_time/07d")

## cerv
cervical.07d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.07d, "cervical_07d_design_file.csv"))

cerv.07d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "07d", "Cervical_deseq")

cervical.de.uqua.notnorm2.07d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.07d.design.matrix, 
                                                   de.test = "DeSeq", results.folder =  cerv.07d.output.results.folder,
                                                   prefix.label = "Cervical 07d prop uqua",
                                                   normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)

# 
# cerv.07d.output.results.folder.noiseq95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "07d", "Cervical_noiseq_95")
# 
# cervical.de.uqua.notnorm2.07d.noi95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cervical.07d.design.matrix, 
#                                                    de.test = "NOISeqBio", results.folder =  cerv.07d.output.results.folder.noiseq95,
#                                                    prefix.label = "Cervical 07d prop uqua",
#                                                    normalize.data.flag = FALSE, threshold = 0.95, conversion.map=ens.symb.biomart.map)


# cerv.07d.output.results.folder.noiseq99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "07d", "Cervical_noiseq_99")
# 
# cervical.de.uqua.notnorm2.07d.noi99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cervical.07d.design.matrix, 
#                                                    de.test = "NOISeqBio", results.folder =  cerv.07d.output.results.folder.noiseq99,
#                                                    prefix.label = "Cervical 07d prop uqua",
#                                                    normalize.data.flag = FALSE, threshold = 0.99, conversion.map=ens.symb.biomart.map)


## thor 
thoracic.07d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.07d, "thoracic_07d_design_file.csv"))

thor.07d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "07d", "Thoracic_deseq")

thoracic.de.uqua.notnorm2.07d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = thoracic.07d.design.matrix, 
                                                   de.test = "DeSeq", results.folder =  thor.07d.output.results.folder,
                                                   prefix.label = "Thoracic 07d prop uqua",
                                                   normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)
# 
# 
# thor.07d.output.results.folder.noi95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "07d", "Thoracic_noiseq_95")
# 
# thoracic.de.uqua.notnorm2.07d.noi95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = thoracic.07d.design.matrix, 
#                                                    de.test = "NOISeqBio", results.folder =  thor.07d.output.results.folder.noi95,
#                                                    prefix.label = "Thoracic 07d prop uqua",
#                                                    normalize.data.flag = FALSE, threshold = 0.95, conversion.map=ens.symb.biomart.map)
# 
# thor.07d.output.results.folder.noi99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "07d", "Thoracic_noiseq_99")
# 
# thoracic.de.uqua.notnorm2.07d.noi99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = thoracic.07d.design.matrix, 
#                                                          de.test = "NOISeqBio", results.folder =  thor.07d.output.results.folder.noi99,
#                                                          prefix.label = "Thoracic 07d prop uqua",
#                                                          normalize.data.flag = FALSE, threshold = 0.99, conversion.map=ens.symb.biomart.map)



## cerv-thor 
cerv.thor.07d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.07d, "cerv-thor_07d_design_file.csv"))

cerv.thor.07d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "07d", "Cerv-Thor_deseq")

cerv.thor.de.uqua.notnorm2.07d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.07d.design.matrix, 
                                                    de.test = "DeSeq", results.folder =  cerv.thor.07d.output.results.folder,
                                                    prefix.label = "Cerv-Thor 07d prop uqua",
                                                    normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)
# 
# cerv.thor.07d.output.results.folder.noiseq95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "07d", "Cerv-Thor_noiseq_95")
# 
# cerv.thor.de.uqua.notnorm2.07d.noi95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cerv.thor.07d.design.matrix, 
#                                                     de.test = "NOISeqBio", results.folder =  cerv.thor.07d.output.results.folder.noiseq95,
#                                                     prefix.label = "Cerv-Thor 07d prop uqua",
#                                                     normalize.data.flag = FALSE, threshold = 0.95, conversion.map=ens.symb.biomart.map)
# 
# cerv.thor.07d.output.results.folder.noiseq99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "07d", "Cerv-Thor_noiseq_99")
# 
# cerv.thor.de.uqua.notnorm2.07d.noi99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cerv.thor.07d.design.matrix, 
#                                                           de.test = "NOISeqBio", results.folder =  cerv.thor.07d.output.results.folder.noiseq99,
#                                                           prefix.label = "Cerv-Thor 07d prop uqua",
#                                                           normalize.data.flag = FALSE, threshold = 0.99, conversion.map=ens.symb.biomart.map)


## cerv-thor untr 

cerv.thor.07d.design.matrix.untr <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.07d, "cervical_thoracic_07d_untr_design_file.txt"))

cerv.thor.07d.output.results.folder.untr <- UpdateFolderPath(output.results.folder, "Time_by_Time", "07d", "Cerv-Thor_untr_deseq")

cerv.thor.de.uqua.notnorm2.07d.deseq.untr <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.07d.design.matrix.untr, 
                                                               de.test = "DeSeq", results.folder =  cerv.thor.07d.output.results.folder.untr,
                                                               prefix.label = "Cerv-Thor untr 07d prop uqua",
                                                               normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)

# ## comparing 07d
# venn.output.plot.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_1/plots/Proportion/uqua/"
# 
# VennDE3ListsSuited(cervical.de.uqua.notnorm2.07d, thoracic.de.uqua.notnorm2.07d, cerv.thor.de.uqua.notnorm2.07d, title = "07 days", venn.plot.dir = venn.output.plot.folder)
# 
# 


############### single point 14d
design.folder.14d <- file.path(design.files.folder, "time_by_time/14d")

## cerv
cervical.14d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.14d, "cervical_14d_design_file.csv"))

cerv.14d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "14d", "Cervical_deseq")

cervical.de.uqua.notnorm2.14d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.14d.design.matrix, 
                                                   de.test = "DeSeq", results.folder =  cerv.14d.output.results.folder,
                                                   prefix.label = "Cervical 14d prop uqua",
                                                   normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)
# 
# cerv.14d.output.results.folder.noi95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "14d", "Cervical_noiseq95")
# 
# cervical.de.uqua.notnorm2.14d.noi95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cervical.14d.design.matrix, 
#                                                    de.test = "NOISeqBio", results.folder =  cerv.14d.output.results.folder.noi95,
#                                                    prefix.label = "Cervical 14d prop uqua",
#                                                    normalize.data.flag = FALSE, threshold = 0.95, conversion.map=ens.symb.biomart.map)
# 
# cerv.14d.output.results.folder.noi99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "14d", "Cervical_noiseq99")
# 
# cervical.de.uqua.notnorm2.14d.noi99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cervical.14d.design.matrix, 
#                                                          de.test = "NOISeqBio", results.folder =  cerv.14d.output.results.folder.noi99,
#                                                          prefix.label = "Cervical 14d prop uqua",
#                                                          normalize.data.flag = FALSE, threshold = 0.99, conversion.map=ens.symb.biomart.map)

## thor 
thoracic.14d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.14d, "thoracic_14d_design_file.csv"))

thor.14d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "14d", "Thoracic_deseq")

thoracic.de.uqua.notnorm2.14d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = thoracic.14d.design.matrix, 
                                                   de.test = "DeSeq", results.folder =  thor.14d.output.results.folder,
                                                   prefix.label = "Thoracic 14d prop uqua",
                                                   normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)
# 
# 
# thor.14d.output.results.folder.noi95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "14d", "Thoracic_noiseq_95")
# 
# thoracic.de.uqua.notnorm2.14d.noi95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = thoracic.14d.design.matrix, 
#                                                    de.test = "NOISeqBio", results.folder =  thor.14d.output.results.folder.noi95,
#                                                    prefix.label = "Thoracic 14d prop uqua",
#                                                    normalize.data.flag = FALSE, threshold = 0.95, conversion.map=ens.symb.biomart.map)
# 
# thor.14d.output.results.folder.noi99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "14d", "Thoracic_noiseq_99")
# 
# thoracic.de.uqua.notnorm2.14d.noi99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = thoracic.14d.design.matrix, 
#                                                    de.test = "NOISeqBio", results.folder =  thor.14d.output.results.folder.noi99,
#                                                    prefix.label = "Thoracic 14d prop uqua",
#                                                    normalize.data.flag = FALSE, threshold = 0.99, conversion.map=ens.symb.biomart.map)


## cerv-thor 
cerv.thor.14d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.14d, "cerv-thor_14d_design_file.csv"))

cerv.thor.14d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "14d", "Cerv-Thor_deseq")

cerv.thor.de.uqua.notnorm2.14d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.14d.design.matrix, 
                                                    de.test = "DeSeq", results.folder =  cerv.thor.14d.output.results.folder,
                                                    prefix.label = "Cerv-Thor 14d prop uqua",
                                                    normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)


# cerv.thor.14d.output.results.folder.noi95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "14d", "Cerv-Thor_noiseq_95")
# 
# cerv.thor.de.uqua.notnorm2.14d.noi95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cerv.thor.14d.design.matrix, 
#                                                     de.test = "NOISeqBio", results.folder =  cerv.thor.14d.output.results.folder.noi95,
#                                                     prefix.label = "Cerv-Thor 14d prop uqua",
#                                                     normalize.data.flag = FALSE, threshold = 0.95, conversion.map=ens.symb.biomart.map)
# 
# cerv.thor.14d.output.results.folder.noi99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "14d", "Cerv-Thor_noiseq_99")
# 
# cerv.thor.de.uqua.notnorm2.14d.noi99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cerv.thor.14d.design.matrix, 
#                                                           de.test = "NOISeqBio", results.folder =  cerv.thor.14d.output.results.folder.noi99,
#                                                           prefix.label = "Cerv-Thor 14d prop uqua",
#                                                           normalize.data.flag = FALSE, threshold = 0.99, conversion.map=ens.symb.biomart.map)


## cerv-thor untr 

cerv.thor.14d.design.matrix.untr <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.14d, "cervical_thoracic_14d_untr_design_file.txt"))

cerv.thor.14d.output.results.folder.untr <- UpdateFolderPath(output.results.folder, "Time_by_Time", "14d", "Cerv-Thor_untr_deseq")

cerv.thor.de.uqua.notnorm2.14d.deseq.untr <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.14d.design.matrix.untr, 
                                                               de.test = "DeSeq", results.folder =  cerv.thor.14d.output.results.folder.untr,
                                                               prefix.label = "Cerv-Thor untr 14d prop uqua",
                                                               normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)

# ## comparing 03d
# venn.output.plot.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_1/plots/Proportion/uqua/"
# 
# VennDE3ListsSuited(cervical.de.uqua.notnorm2.14d, thoracic.de.uqua.notnorm2.14d, cerv.thor.de.uqua.notnorm2.14d, title = "14 days", venn.plot.dir = venn.output.plot.folder)




############### single point 56d
design.folder.56d <- file.path(design.files.folder, "time_by_time/56d")

## cerv
cervical.56d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.56d, "cervical_56d_design_file.csv"))

cerv.56d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "56d", "Cervical_deseq")

cervical.de.uqua.notnorm2.56d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.56d.design.matrix, 
                                                   de.test = "DeSeq", results.folder =  cerv.56d.output.results.folder,
                                                   prefix.label = "Cervical 56d prop uqua",
                                                   normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)

# 
# cerv.56d.output.results.folder.noi95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "56d", "Cervical_noiseq_95")
# 
# cervical.de.uqua.notnorm2.56d.noi95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cervical.56d.design.matrix, 
#                                                    de.test = "NOISeqBio", results.folder =  cerv.56d.output.results.folder.noi95,
#                                                    prefix.label = "Cervical 56d prop uqua",
#                                                    normalize.data.flag = FALSE, threshold = 0.95, conversion.map=ens.symb.biomart.map)

# cerv.56d.output.results.folder.noi99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "56d", "Cervical_noiseq_99")
# 
# cervical.de.uqua.notnorm2.56d.noi99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cervical.56d.design.matrix, 
#                                                          de.test = "NOISeqBio", results.folder =  cerv.56d.output.results.folder.noi99,
#                                                          prefix.label = "Cervical 56d prop uqua",
#                                                          normalize.data.flag = FALSE, threshold = 0.99, conversion.map=ens.symb.biomart.map)



## thor 
thoracic.56d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.56d, "thoracic_56d_design_file.csv"))

thor.56d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "56d", "Thoracic_deseq")

thoracic.de.uqua.notnorm2.56d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = thoracic.56d.design.matrix, 
                                                   de.test = "DeSeq", results.folder =  thor.56d.output.results.folder,
                                                   prefix.label = "Thoracic 56d prop uqua",
                                                   normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)

# whole.counts = round(normalized.counts.prop.uqua); design.matrix = thoracic.56d.design.matrix; de.test = "DeSeq"; results.folder =  thor.56d.output.results.folder; prefix.label = "Thoracic 56d prop uqua" ; normalize.data.flag = FALSE
# 
# 
# thor.56d.output.results.folder.noi95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "56d", "Thoracic_noiseq_95")
# 
# thoracic.de.uqua.notnorm2.56d.noi95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = thoracic.56d.design.matrix, 
#                                                    de.test = "NOISeqBio", results.folder =  thor.56d.output.results.folder.noi95,
#                                                    prefix.label = "Thoracic 56d prop uqua",
#                                                    normalize.data.flag = FALSE, threshold = 0.95)
# 
# thor.56d.output.results.folder.noi99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "56d", "Thoracic_noiseq_99")
# 
# thoracic.de.uqua.notnorm2.56d.noi99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = thoracic.56d.design.matrix, 
#                                                          de.test = "NOISeqBio", results.folder =  thor.56d.output.results.folder.noi99,
#                                                          prefix.label = "Thoracic 56d prop uqua",
#                                                          normalize.data.flag = FALSE, threshold = 0.99)



## cerv-thor 
cerv.thor.56d.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.56d, "cerv-thor_56d_design_file.csv"))

cerv.thor.56d.output.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time", "56d", "Cerv-Thor_deseq")

cerv.thor.de.uqua.notnorm2.56d.deseq <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.56d.design.matrix, 
                                                    de.test = "DeSeq", results.folder =  cerv.thor.56d.output.results.folder,
                                                    prefix.label = "Cerv-Thor 56d prop uqua",
                                                    normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)

# 
# cerv.thor.56d.output.results.folder.noi95 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "56d", "Cerv-Thor_noiseqbio_95")
# 
# cerv.thor.de.uqua.notnorm2.56d.noi95 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cerv.thor.56d.design.matrix, 
#                                                     de.test = "NOISeqBio", results.folder =  cerv.thor.56d.output.results.folder.noi95,
#                                                     prefix.label = "Cerv-Thor 56d prop uqua",
#                                                     normalize.data.flag = FALSE, threshold = 0.95, conversion.map=ens.symb.biomart.map)
# 
# cerv.thor.56d.output.results.folder.noi99 <- UpdateFolderPath(output.results.folder, "Time_by_Time", "56d", "Cerv-Thor_noiseqbio_99")
# 
# cerv.thor.de.uqua.notnorm2.56d.noi99 <- PerformDEAnalysis(whole.counts = normalized.counts.prop.uqua, design.matrix = cerv.thor.56d.design.matrix, 
#                                                           de.test = "NOISeqBio", results.folder =  cerv.thor.56d.output.results.folder.noi99,
#                                                           prefix.label = "Cerv-Thor 56d prop uqua",
#                                                           normalize.data.flag = FALSE, threshold = 0.99, conversion.map=ens.symb.biomart.map)



## cerv-thor untr 

cerv.thor.56d.design.matrix.untr <- ReadDataFrameFromTsv(file.name.path = file.path(design.folder.56d, "cervical_thoracic_56d_untr_design_file.txt"))

cerv.thor.56d.output.results.folder.untr <- UpdateFolderPath(output.results.folder, "Time_by_Time", "56d", "Cerv-Thor_untr_deseq")

cerv.thor.de.uqua.notnorm2.56d.deseq.untr <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.56d.design.matrix.untr, 
                                                               de.test = "DeSeq", results.folder =  cerv.thor.56d.output.results.folder.untr,
                                                               prefix.label = "Cerv-Thor untr 56d prop uqua",
                                                               normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)


# 
# ## comparing 56d
# venn.output.plot.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_1/plots/Proportion/uqua/"
# 
# VennDE3ListsSuited(cervical.de.uqua.notnorm2.56d, thoracic.de.uqua.notnorm2.56d, cerv.thor.de.uqua.notnorm2.56d, title = "56 days", venn.plot.dir = venn.output.plot.folder)
# 
