

source("/src/TimeCourseFunctions.R")
source("./src/MainFunction.R")
source("./src/BoxplotFunctions.R")

setwd(dir = "/media/dario/dati/time_course_4/")

all.bam.folder <- "/media/dario/dati/time_course_4/cervical_thoracic/all_bams/"
# bam.files <- list.files(cerv.thor.bam.folder)
# all.bam <- bam.files[grep(pattern = "*.bam", bam.files)]
# just.bam <- all.bam[-grep(pattern = "*.bai", all.bam)]

bam.files <- list.files(all.bam.folder, pattern="bam$",full.names=TRUE)


rat5.ensembl.downloaded.gtf <- "/media/dario/dati/time_course_4/downloaded_references/illumina_igenome_site/Rattus_norvegicus_Ensembl_Rnor_5.0/Rattus_norvegicus/Ensembl/Rnor_5.0/Annotation/Genes/genes.gtf"
output.folder <- file.path("/media/dario/dati/time_course_4/cervical_thoracic/results_2018")
dir.create(output.folder)

counts.output.folder <- file.path(output.folder, "counts")

counts.ens.gtf <- countBamFilesFeatureCounts(bam.files = bam.files, annotation.file = rat5.ensembl.downloaded.gtf, counts.folder = counts.output.folder, gtf.attr.type = "gene_id", prefix.output.file="Cervical_Thoracic")
colnames(counts.ens.gtf)

counts.ens.gtf.renamed <- counts.ens.gtf
colnames(counts.ens.gtf.renamed) <- c( "cervical_03d_t1", "cervical_03d_t2", "cervical_03d_t3", "cervical_03d_t4", "cervical_03d_t5",
                                        "cervical_07d_t1", "cervical_07d_t2", "cervical_07d_t3", "cervical_07d_t4", "cervical_07d_t5",
                                        "cervical_14d_t1", "cervical_14d_t2", "cervical_14d_t3", "cervical_14d_t4", "cervical_14d_t5",
                                        "cervical_56d_t1", "cervical_56d_t2", "cervical_56d_t3", "cervical_56d_t4", "cervical_56d_t5",
                                        "thoracic_03d_t1", "thoracic_03d_t2", "thoracic_03d_t3", "thoracic_03d_t4", "thoracic_03d_t5",
                                        "thoracic_07d_t1", "thoracic_07d_t2", "thoracic_07d_t3", "thoracic_07d_t4", "thoracic_07d_t5",
                                        "thoracic_14d_t1", "thoracic_14d_t2", "thoracic_14d_t3", "thoracic_14d_t4", "thoracic_14d_t5",
                                        "thoracic_56d_t1", "thoracic_56d_t2", "thoracic_56d_t3", "thoracic_56d_t4", "thoracic_56d_t5"
)
colnames(counts.ens.gtf.renamed)
WriteDataFrameAsTsv(data.frame.to.save = counts.ens.gtf.renamed, file.name.path = file.path(counts.output.folder, "Cervical_Thoracic_Counts_Renamed_FeatureCounts"))

design.file <- "/media/dario/dati/time_course/cervical_thoracic/design_file/cervical_thoracic_design_file.txt"
design.matrix <- read.table(design.file, header = TRUE, row.names = 1)

de.genes.proportion.uqua <- MainFunction(counts.dataframe = counts.ens.gtf.renamed, design.matrix = design.matrix, output.results.folder = output.folder, prefix.output.file = "cervical_thoracic", prefix.plot.label = "cerv_thor", is.time.course = TRUE, filter.method = "Proportion")

unnormalized.counts <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic/results/counts/FeatureCounts/cervical_thoracic__filtered_counts__Proportion.tsv")
head(unnormalized.counts)
PlotTimesBoxplot(data.frame.to.plot = unnormalized.counts, design.matrix = design.matrix, output.path = file.path(output.folder, "plots"), prefix.plot = "cerv_thor un-normalized")

normalized.counts <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic/results/counts/FeatureCounts/cervical_thoracic__filtered_counts__Proportion__normalized__uqua.tsv")
head(normalized.counts)
PlotTimesBoxplot(data.frame.to.plot = normalized.counts, design.matrix = design.matrix, output.path = file.path(output.folder, "plots"), prefix.plot = "cerv_thor normalized uqua")

# PlotPCAFunction(samples.dataframe = normalized.counts, design.dataframe = design.matrix,scale = TRUE)

sign.de.genes.proportion.uqua <- SignificantDeGenesPAdj(de.genes.proportion.uqua)
WriteDataFrameAsTsv(sign.de.genes.proportion.uqua, file.name.path = file.path(output.folder, "de_genes", "DESeq2", "cervical_thoracic__filtered_counts__Proportion__normalized__uqua_de_genes_LRT_sign_005"))

proportion.estimated.neg.genes <- EstimateNegativeControlGenesForRUV(de.genes = de.genes.proportion.uqua, n.tail.genes = 2000, counts.dataset = unnormalized.counts, n.genes.per.hist.break = 1)
# thoracic.proportion.unnormalized<-read.table("/media/dario/dati/time_course/Thoracic/results_backup_18_01_2016/counts/FeatureCounts/thoracic__filtered_counts__Proportion.tsv", sep="\t", row.names = 1, header = TRUE)
# ReadDataFrameFromTsv("/media/dario/dati/time_course/Thoracic/results/counts/FeatureCounts/thoracic_ensgtf_reordered_counts.tsv")
de.genes.proportion.ruvg <- MainFunction(counts.dataframe = counts.ens.gtf.renamed, design.matrix = design.matrix, output.results.folder = output.folder, prefix.output.file = "cerv_thor_treat", prefix.plot.label = "cerv_thor_treat", is.time.course = TRUE, filter.method = "Proportion", normalization.method = "ruvg", estimated.genes = proportion.estimated.neg.genes)
sign.de.genes.proportion.ruvg <- SignificantDeGenesPAdj(de.genes.proportion.ruvg)
WriteDataFrameAsTsv(sign.de.genes.proportion.ruvg, file.name.path = file.path(output.folder, "de_genes", "DESeq2", "cervical_thoracic__filtered_counts__Proportion__normalized__ruvg_de_genes_LRT_sign_005"))


####################### UNTREATED
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/TimeCourseFunctions.R")
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/MainFunction.R")
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/BoxplotFunctions.R")

setwd(dir = "/media/dario/dati/time_course/")

cerv.thor.bam.folder.u <- "/media/dario/dati/time_course/cervical_thoracic/bam_u"
# bam.files <- list.files(cerv.thor.bam.folder)
# all.bam <- bam.files[grep(pattern = "*.bam", bam.files)]
# just.bam <- all.bam[-grep(pattern = "*.bai", all.bam)]

bam.files.u <- list.files(cerv.thor.bam.folder.u, pattern="bam$",full.names=TRUE)


rat5.ensembl.downloaded.gtf <- "/media/dario/dati/time_course/downloaded_references/illumina_igenome_site/Rattus_norvegicus_Ensembl_Rnor_5.0/Rattus_norvegicus/Ensembl/Rnor_5.0/Annotation/Genes/genes.gtf"
output.folder.u <- file.path("/media/dario/dati/time_course/cervical_thoracic/results_u")
dir.create(output.folder.u)

counts.output.folder.u <- file.path(output.folder.u, "counts")

counts.ens.gtf.u <- CountBamFilesFeatureCounts(bam.files = bam.files.u, annotation.file = rat5.ensembl.downloaded.gtf, counts.folder = counts.output.folder.u, gtf.attr.type = "gene_id", prefix.output.file="Cervical_Thoracic_Untr")
colnames(counts.ens.gtf.u)

counts.ens.gtf.renamed.u <- counts.ens.gtf.u
colnames(counts.ens.gtf.renamed.u) <- c(
                                         "cervical_03d_u1", "cervical_03d_u2", "cervical_03d_u3",
                                         "cervical_07d_u1", "cervical_07d_u2", "cervical_07d_u3",
                                         "cervical_14d_u1", "cervical_14d_u2", "cervical_14d_u3",
                                         "cervical_56d_u1", "cervical_56d_u2", "cervical_56d_u3",
                                         "thoracic_03d_u1", "thoracic_03d_u2", "thoracic_03d_u3",
                                         "thoracic_07d_u1", "thoracic_07d_u2", "thoracic_07d_u3",
                                         "thoracic_14d_u1", "thoracic_14d_u2", "thoracic_14d_u3",
                                         "thoracic_56d_u1"
)
colnames(counts.ens.gtf.renamed.u)
WriteDataFrameAsTsv(data.frame.to.save = counts.ens.gtf.renamed.u, file.name.path = file.path(counts.output.folder.u, "Cervical_Thoracic_Untr_Counts_Renamed_FeatureCounts"))

design.file.u <- "/media/dario/dati/time_course/cervical_thoracic/design_file/cervical_thoracic_untr_design_file.txt"
design.matrix.u <- read.table(design.file.u, header = TRUE, row.names = 1)

de.genes.proportion.uqua.u <- MainFunction(counts.dataframe = counts.ens.gtf.renamed.u, design.matrix = design.matrix.u, output.results.folder = output.folder.u, prefix.output.file = "cervical_thoracic_untr", prefix.plot.label = "cerv_thor_untr", is.time.course = TRUE, filter.method = "Proportion")

unnormalized.counts.u <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic/results_u/counts/FeatureCounts/cervical_thoracic_untr__filtered_counts__Proportion.tsv")
head(unnormalized.counts.u)
PlotTimesBoxplot(data.frame.to.plot = unnormalized.counts.u, design.matrix = design.matrix.u, output.path = file.path(output.folder.u, "plots"), prefix.plot = "cerv_thor_untr un-normalized")

normalized.counts.u <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic/results_u/counts/FeatureCounts/cervical_thoracic_untr__filtered_counts__Proportion__normalized__uqua.tsv")
head(normalized.counts.u)
PlotTimesBoxplot(data.frame.to.plot = normalized.counts.u, design.matrix = design.matrix.u, output.path = file.path(output.folder.u, "plots"), prefix.plot = "cerv_thor_untr normalized uqua")

# PlotPCAFunction(samples.dataframe = normalized.counts, design.dataframe = design.matrix,scale = TRUE)

sign.de.genes.proportion.uqua.u <- SignificantDeGenesPAdj(de.genes.proportion.uqua.u)
WriteDataFrameAsTsv(sign.de.genes.proportion.uqua.u, file.name.path = file.path(output.folder.u, "de_genes", "DESeq2", "cervical_thoracic_untr_filtered_counts__Proportion__normalized__uqua_de_genes_LRT_sign_005"))

proportion.estimated.neg.genes.u <- EstimateNegativeControlGenesForRUV(de.genes = de.genes.proportion.uqua.u, n.tail.genes = 2000, counts.dataset = unnormalized.counts.u, n.genes.per.hist.break = 1)
# thoracic.proportion.unnormalized<-read.table("/media/dario/dati/time_course/Thoracic/results_backup_18_01_2016/counts/FeatureCounts/thoracic__filtered_counts__Proportion.tsv", sep="\t", row.names = 1, header = TRUE)
# ReadDataFrameFromTsv("/media/dario/dati/time_course/Thoracic/results/counts/FeatureCounts/thoracic_ensgtf_reordered_counts.tsv")
de.genes.proportion.ruvg.u <- MainFunction(counts.dataframe = counts.ens.gtf.renamed.u, design.matrix = design.matrix.u, output.results.folder = output.folder.u, prefix.output.file = "cerv_thor_untr", prefix.plot.label = "cerv_thor_untr", is.time.course = TRUE, filter.method = "Proportion", normalization.method = "ruvg", estimated.genes = proportion.estimated.neg.genes.u)
sign.de.genes.proportion.ruvg.u <- SignificantDeGenesPAdj(de.genes.proportion.ruvg.u)
WriteDataFrameAsTsv(sign.de.genes.proportion.ruvg.u, file.name.path = file.path(output.folder.u, "de_genes", "DESeq2", "cervical_thoracic_untr_filtered_counts__Proportion__normalized__ruvg_de_genes_LRT_sign_005"))
