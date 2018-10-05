
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/includes.R")

output.results.folder <- "cervical_thoracic_whole/results_DESeqTests"

all.design.file <- "cervical_thoracic/design_file/cervical_thoracic_all_design_file.txt"
cer.thor.des.mat <- ReadDataFrameFromTsv(all.design.file)


design.files.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/design_file"

ens.symb.biomart.map <- ReadDataFrameFromTsv("downloaded_references/biomart/rnor5.0_biomart_ensembls_associatedgenenames.txt", row.names.col = NULL)
head(ens.symb.biomart.map)

cervical.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/design_file/4timepoints/cervical_design_file.txt")

cerv.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Cervical")

cerv.de.uqua.notnorm2.des4time.TC <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
                                                                   de.test = "DeSeqTime_TC", results.folder =  cerv.output.results.folder,
                                                                   prefix.label = "4Points_Cervical prop uqua",
                                                                   normalize.data.flag = FALSE, threshold = 0.05, conversion.map=ens.symb.biomart.map)

cerv.de.uqua.notnorm2.des4time.T <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
                                                                   de.test = "DeSeqTime_T", results.folder =  cerv.output.results.folder,
                                                                   prefix.label = "4Points_Cervical prop uqua",
                                                                   normalize.data.flag = FALSE, threshold = 0.05, conversion.map=ens.symb.biomart.map)

cerv.de.uqua.notnorm2.des4time.NoInt <- PerformDEAnalysis(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
                                                                de.test = "DeSeqTime_NoInteraction", results.folder =  cerv.output.results.folder,
                                                                prefix.label = "4Points_Cervical prop uqua",
                                                                normalize.data.flag = FALSE, threshold = 0.05, conversion.map=ens.symb.biomart.map)


cerv.de.uqua.notnorm2.des4time.TC$de.total[which(rownames(cerv.de.uqua.notnorm2.des4time.TC$de.total) %in% "ENSRNOG00000019689"),]

cerv.norm.counts <- round(normalized.counts.prop.uqua[,which(colnames(normalized.counts.prop.uqua) %in% rownames(cervical.design.matrix))])
cerv.norm.counts[which(rownames(cerv.norm.counts) %in% "ENSRNOG00000019689"),]

n.03.07 <- union(rownames(cervical.de.uqua.notnorm2.03d.deseq$SIGN), rownames(cervical.de.uqua.notnorm2.07d$SIGN))
length(n.03.07)
n.03.07.14 <- union(n.03.07, rownames(cervical.de.uqua.notnorm2.14d$SIGN))
length(n.03.07.14)
n.03.07.14.56 <- unique(union(n.03.07, rownames(cervical.de.uqua.notnorm2.56d$SIGN)))
length(n.03.07.14.56)


cervical.design.matrix
length(intersect(rownames(cerv.de.uqua.notnorm2.des4time.NoInt$SIGN), rownames(cerv.de.uqua.notnorm2.des4time.TC$SIGN)))


mcols(de.results)


length(rownames(cervical.de.uqua.notnorm2.03d.deseq$SIGN))
cond.tr.vs.untr<-results(de.results, name="Conditions_treated_vs_untreated", test="Wald")
cond.tr.vs.untr.sign <- subset(cond.tr.vs.untr, padj<0.05)
length(rownames(cond.tr.vs.untr.sign))

length(intersect( rownames(cervical.de.uqua.notnorm2.03d.deseq$SIGN) , rownames(cond.tr.vs.untr.sign)))
length(intersect( n.03.07.14.56 , rownames(cond.tr.vs.untr.sign)))





de.outliers <- rownames(readable.results)[which(is.na(readable.results$padj))]

ma.outliers <- colnames(cerv.p.masig.t$influ.info)

length(intersect(de.outliers, ma.outliers))















