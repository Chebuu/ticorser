

source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/includes.R")

date <- gsub(pattern = " ", replacement = "_", date())
date <- gsub(pattern = ":", replacement = "_", date)
output.results.folder <- paste0("cervical_thoracic_whole/", date)
all.design.file <- "/media/dario/dati/time_course/cervical_thoracic/design_file/cervical_thoracic_all_design_file.txt"
all.des.mat <- ReadDataFrameFromTsv(all.design.file)

design.files.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/design_file"

normalized.counts.prop.uqua <- ReadDataFrameFromTsv("cervical_thoracic_whole/results_final_10/counts/FeatureCounts/Proportion/uqua/cervical_thoracic_whole_Proportion_uqua.tsv")
ens.symb.biomart.map <- ReadDataFrameFromTsv("downloaded_references/biomart/rnor5.0_biomart_ensembls_associatedgenenames.txt", row.names.col = NULL)

## cerv 4times
cervical.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/design_file/4timepoints/cervical_design_file.txt")

cerv.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Cervical")


# whole.counts = round(normalized.counts.prop.uqua); design.matrix = cervical.design.matrix; de.test = "DeSeqTime_TC"; results.folder =  cerv.output.results.folder;  prefix.label = "4Points_Cervical prop uqua";  normalize.data.flag = FALSE; threshold = 0.05, conversion.map=ens.symb.biomart.map

cerv.de.uqua.notnorm2.des4time.lrt.wald <- PerformDEAnalysisLRT(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cervical.design.matrix, 
                                                                  de.test = "DeSeqTime_TC", results.folder =  cerv.output.results.folder,
                                                                  prefix.label = "4Points_Cervical prop uqua",
                                                                  normalize.data.flag = FALSE, threshold = 0.05, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)


beep()

## thor 4times
thor.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Thoracic")
thoracic.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/design_file/4timepoints/thoracic_design_file.txt")

thor.de.uqua.notnorm2.des4time.lrt.wald <- PerformDEAnalysisLRT(whole.counts = round(normalized.counts.prop.uqua), design.matrix = thoracic.design.matrix, 
                                                    de.test = "DeSeqTime_TC", results.folder =  thor.output.results.folder,
                                                    prefix.label = "4Points_Thoracic prop uqua",
                                                    normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)

## cerv-thor 4times
cerv.thor.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/design_file/4timepoints/cervical_thoracic_design_file.txt")
cerv.thor.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Cerv-Thor")


cerv.thor.de.uqua.notnorm2.des4time.lrt.wald <- PerformDEAnalysisLRT(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.design.matrix, 
                                                         de.test = "DeSeqTime_TC", results.folder =  cerv.thor.output.results.folder,
                                                         prefix.label = "4Points_Cerv_Thor prop uqua",
                                                         normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)



## cerv-thor 4times UNTREATED
cerv.thor.untr.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/design_file/4timepoints/cervical_thoracic_untr_design_file.txt")
cerv.thor.untr.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Cerv-Thor_Untr")


cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald <- PerformDEAnalysisLRT(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.untr.design.matrix, 
                                                                     de.test = "DeSeqTime_TC", results.folder =  cerv.thor.untr.output.results.folder,
                                                                     prefix.label = "4Points_Cerv_Thor_untr prop uqua",
                                                                     normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, enrich.results.flag = FALSE)


graphics.off()
design.files.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/design_file"
## cerv 3times
subacute.design.folder <- file.path(design.files.folder, "subacute")

subacute.cervical.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(subacute.design.folder, "cervical_subacute_design_file.csv"))

cerv.subacute.output.results.folder <- UpdateFolderPath(output.results.folder, "subacute_Cervical")

cervical.de.uqua.notnorm2.subacute.lrt.wald <- PerformDEAnalysisLRT(whole.counts = round(normalized.counts.prop.uqua), design.matrix = subacute.cervical.design.matrix, 
                                                        de.test = "DeSeqTime_TC", results.folder =  cerv.subacute.output.results.folder,
                                                        prefix.label = "Cervical subacute prop uqua",
                                                        normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map)


## thoracic subacute
subacute.thoracic.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(subacute.design.folder, "thoracic_subacute_design_file.csv"))

thor.subacute.output.results.folder <- UpdateFolderPath(output.results.folder, "subacute_Thoracic")

thoracic.de.uqua.notnorm2.subacute.lrt.wald <- PerformDEAnalysisLRT(whole.counts = round(normalized.counts.prop.uqua), design.matrix = subacute.thoracic.design.matrix, 
                                                        de.test = "DeSeqTime_TC", results.folder =  thor.subacute.output.results.folder,
                                                        prefix.label = "Thoracic subacute prop uqua",
                                                        normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map)

## cerv-thor subacute
subacute.cerv.thor.design.matrix <- ReadDataFrameFromTsv(file.name.path = file.path(subacute.design.folder, "cerv-thor_subacute_design_file.csv"))

cerv.thor.subacute.output.results.folder <- UpdateFolderPath(output.results.folder, "subacute_Cerv-Thor")

cerv.thor.de.uqua.notnorm2.subacute.lrt.wald <- PerformDEAnalysisLRT(whole.counts = round(normalized.counts.prop.uqua), design.matrix = subacute.cerv.thor.design.matrix, 
                                                         de.test = "DeSeqTime_TC", results.folder =  cerv.thor.subacute.output.results.folder,
                                                         prefix.label = "Cerv-Thor subacute prop uqua",
                                                         normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map)


## cerv-thor subacute UNTREATED
# output.results.folder <- "cervical_thoracic_whole/results_final_20"
cerv.thor.subacute.untr.design.matrix <- ReadDataFrameFromTsv(file.name.path = "/media/dario/dati/time_course/cervical_thoracic_whole/design_file/subacute/cervical_thoracic_untr_design_file.txt")
cerv.thor.subacute.untr.output.results.folder <- UpdateFolderPath(output.results.folder, "subacute_Cerv-Thor_Untr")


cerv.thor.untr.de.uqua.notnorm2.subacute.lrt.wald <- PerformDEAnalysisLRT(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.subacute.untr.design.matrix, 
                                                                          de.test = "DeSeqTime_TC", results.folder =  cerv.thor.subacute.untr.output.results.folder,
                                                                          prefix.label = "Cerv_Thor_untr subacute prop uqua",
                                                                          normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map)
