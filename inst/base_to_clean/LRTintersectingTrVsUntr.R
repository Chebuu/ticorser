cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald

cerv.thor.de.uqua.notnorm2.des4time.lrt.wald


intersecting.cerv.thor.tr.vs.untr.lrt <- intersect(rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN), rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN)) ##0

intersecting.cerv.thor.tr.vs.untr.wald <- intersect(rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN), rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN))

length(intersecting.cerv.thor.tr.vs.untr.wald) #486
length(rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN)) #1395
length(rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN)) # 2289

intersecting.cerv.thor.tr.minus.untr.wald <- setdiff(rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN), rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN))
length(intersecting.cerv.thor.tr.minus.untr.wald) #1803

functional.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_minus_Untr/functional/WaldLRT_Tr_minus_Untr/"
prefix.label = "4Points_Cerv_Thor_tr_vs_untr prop uqua WaldLRT_tr_minus_untr"
enrichSuitedSingleList(de.gene.list = intersecting.cerv.thor.tr.minus.untr.wald, functional.folder = functional.folder, filename = prefix.label)



functional.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_intersect_Untr/functional/WaldLRT_Tr_intersect_Untr/"
prefix.label = "4Points_Cerv_Thor_tr_vs_untr prop uqua WaldLRT_tr_intersect_untr"
enrichSuitedSingleList(de.gene.list = intersecting.cerv.thor.tr.vs.untr.wald, functional.folder = functional.folder, filename = prefix.label)



normalized.counts.prop.uqua.tr.minus.untr.intersection <- normalized.counts.prop.uqua[-which(rownames(normalized.counts.prop.uqua) %in% intersecting.cerv.thor.tr.vs.untr.wald), ]
dim(normalized.counts.prop.uqua.tr.minus.untr.intersection)
dim(normalized.counts.prop.uqua)

all.design.file <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/design_file/whole/cervical_thoracic_all_design_file.txt")
subsetted.counts <- normalized.counts.prop.uqua[ which(rownames(normalized.counts.prop.uqua) %in% intersecting.cerv.thor.tr.minus.untr.wald),]
subsetted.counts <- AttachConvertedColumn(subsetted.counts, ens.symb.biomart.map)
# normalized.counts = subsetted.counts; design.matrix = all.design.file; gene.name = rownames(subsetted.counts)[1]; gene.name.column.name = "symbol"; show.plot.flag = FALSE; plotly.flag = TRUE; save.plot = TRUE; plot.folder = "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_minus_Untr/profile_plotting"
for(i in 1:dim(subsetted.counts)[1]) {
  tryCatch({
    PlotCountsAlongTimes(normalized.counts = subsetted.counts, design.matrix = all.design.file, gene.name = subsetted.counts$symbol[i], gene.name.column.name = "symbol", show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = "cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_minus_Untr/profile_plotting")
  }, error=function(e) {
    print(e)
  })
  
}

subsetted.counts <- normalized.counts.prop.uqua[ which(rownames(normalized.counts.prop.uqua) %in% intersecting.cerv.thor.tr.vs.untr.wald), ]
subsetted.counts <- AttachConvertedColumn(subsetted.counts, ens.symb.biomart.map)

for(i in 1:dim(subsetted.counts)[1]) {
  tryCatch({
    PlotCountsAlongTimes(normalized.counts = subsetted.counts, design.matrix = all.design.file, gene.name = subsetted.counts$symbol[i], gene.name.column.name = "symbol", show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = "cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_minus_Untr/intersection_genes", prefix.plot = "intersection gene")
  }, error=function(e) {
    print(e)
  })
  
}


## cerv-thor 4times UNTREATED 0.1
cerv.thor.untr.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/design_file/4timepoints/cervical_thoracic_untr_design_file.txt")
cerv.thor.untr.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Cerv-Thor_Untr_thr01")


cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01 <- PerformDEAnalysisLRT(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.untr.design.matrix, 
                                                                          de.test = "DeSeqTime_TC", results.folder =  cerv.thor.untr.output.results.folder,
                                                                          prefix.label = "4Points_Cerv_Thor_untr prop uqua",
                                                                          normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, threshold = 0.1, enrich.results.flag = FALSE)

dim(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01$waldLRT$SIGN) ## 2085

intersecting.cerv.thor.tr.vs.untr.lrt <- intersect(rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01$LRT$SIGN), rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN)) ##0

intersecting.cerv.thor.tr.vs.untr.wald.thr01 <- intersect(rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01$waldLRT$SIGN), rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN))
length(intersecting.cerv.thor.tr.vs.untr.wald.thr01) ## 645

intersecting.cerv.thor.tr.minus.untr.wald.thr01 <- setdiff(rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN), rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01$waldLRT$SIGN))
length(intersecting.cerv.thor.tr.minus.untr.wald.thr01)

subsetted.counts <- normalized.counts.prop.uqua[ which(rownames(normalized.counts.prop.uqua) %in% intersecting.cerv.thor.tr.minus.untr.wald.thr01), ]
subsetted.counts <- AttachConvertedColumn(subsetted.counts, ens.symb.biomart.map)

for(i in 1:dim(subsetted.counts)[1]) {
  tryCatch({
    PlotCountsAlongTimes(normalized.counts = subsetted.counts, design.matrix = all.design.file, gene.name = subsetted.counts$symbol[i], gene.name.column.name = "symbol", show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = "cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_intersect_Untr_thr01/intersection_genes", prefix.plot = "intersection gene")
  }, error=function(e) {
    print(e)
  })
  
}
# 
# functional.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_intersect_Untr_thr01/functional/WaldLRT_Tr_intersect_Untr/"
# prefix.label = "4Points_Cerv_Thor_tr_vs_untr prop uqua WaldLRT_tr_intersect_untr_thr01"
# enrichSuitedSingleList(de.gene.list = intersecting.cerv.thor.tr.minus.untr.wald.thr01, functional.folder = functional.folder, filename = prefix.label)





## cerv-thor 4times UNTREATED 0.2
cerv.thor.untr.design.matrix <- ReadDataFrameFromTsv(file.name.path = "cervical_thoracic_whole/design_file/4timepoints/cervical_thoracic_untr_design_file.txt")
cerv.thor.untr.output.results.folder <- UpdateFolderPath(output.results.folder, "4Points_Cerv-Thor_Untr_thr02")


cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr02 <- PerformDEAnalysisLRT(whole.counts = round(normalized.counts.prop.uqua), design.matrix = cerv.thor.untr.design.matrix, 
                                                                                de.test = "DeSeqTime_TC", results.folder =  cerv.thor.untr.output.results.folder,
                                                                                prefix.label = "4Points_Cerv_Thor_untr prop uqua",
                                                                                normalize.data.flag = FALSE, conversion.map=ens.symb.biomart.map, threshold = 0.2, enrich.results.flag = FALSE)

dim(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr02$waldLRT$SIGN) ## 3184

intersecting.cerv.thor.tr.vs.untr.lrt <- intersect(rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr02$LRT$SIGN), rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN)) ##0

intersecting.cerv.thor.tr.vs.untr.wald.thr02 <- intersect(rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr02$waldLRT$SIGN), rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN))
length(intersecting.cerv.thor.tr.vs.untr.wald.thr02) ## 856

intersecting.cerv.thor.tr.minus.untr.wald.thr02 <- setdiff(rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN), rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr02$waldLRT$SIGN))
length(intersecting.cerv.thor.tr.minus.untr.wald.thr02)

subsetted.counts <- normalized.counts.prop.uqua[ which(rownames(normalized.counts.prop.uqua) %in% intersecting.cerv.thor.tr.minus.untr.wald.thr02), ]
subsetted.counts <- AttachConvertedColumn(subsetted.counts, ens.symb.biomart.map)

for(i in 1:dim(subsetted.counts)[1]) {
  tryCatch({
    PlotCountsAlongTimes(normalized.counts = subsetted.counts, design.matrix = all.design.file, gene.name = subsetted.counts$symbol[i], gene.name.column.name = "symbol", show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = "cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_intersect_Untr_thr02/intersection_genes", prefix.plot = "intersection gene")
  }, error=function(e) {
    print(e)
  })
  
}

functional.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_intersect_Untr_thr02/functional/WaldLRT_Tr_intersect_Untr/"
prefix.label = "4Points_Cerv_Thor_tr_vs_untr prop uqua WaldLRT_tr_intersect_untr_thr01"
enrichSuitedSingleList(de.gene.list = intersecting.cerv.thor.tr.minus.untr.wald.thr01, functional.folder = functional.folder, filename = prefix.label)




###

cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.pval.minorthr005 <- cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01$waldLRT$de.total[cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01$waldLRT$de.total$pvalue < 0.05,]

cerv.thor.de.uqua.notnorm2.des4time.wald.minus.pval.005 <- setdiff(rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN), rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.pval.minorthr005)) 
length(cerv.thor.de.uqua.notnorm2.des4time.wald.minus.pval.005)## 1416



subsetted.counts <- normalized.counts.prop.uqua[ which(rownames(normalized.counts.prop.uqua) %in% cerv.thor.de.uqua.notnorm2.des4time.lrt.wald.minus.pval.005), ]
subsetted.counts <- AttachConvertedColumn(subsetted.counts, ens.symb.biomart.map)

for(i in 1:dim(subsetted.counts)[1]) {
  tryCatch({
    PlotCountsAlongTimes(normalized.counts = subsetted.counts, design.matrix = all.design.file, gene.name = subsetted.counts$symbol[i], gene.name.column.name = "symbol", show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = "cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_intersect_Untr_minus_pval_005/intersection_genes", prefix.plot = "intersection gene")
  }, error=function(e) {
    print(e)
  })
  
}


functional.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/4Points_Cerv-Thor_Tr_intersect_Untr_pval_005/functional/WaldLRT_Tr_intersect_Untr/"
prefix.label = "4Points_Cerv_Thor_tr_vs_untr prop uqua WaldLRT_tr_intersect_untr_pval_005"
enrichSuitedSingleList(de.gene.list = cerv.thor.de.uqua.notnorm2.des4time.lrt.wald.minus.pval.005, functional.folder = functional.folder, filename = prefix.label)

############## subacute intersections


cerv.thor.de.uqua.notnorm2.subacute.lrt.wald
cerv.thor.untr.de.uqua.notnorm2.subacute.lrt.wald

cerv.thor.untr.de.uqua.notnorm2.subacute.lrt.wald.pval.minorthr005 <- cerv.thor.untr.de.uqua.notnorm2.subacute.lrt.wald$waldLRT$de.total[cerv.thor.untr.de.uqua.notnorm2.subacute.lrt.wald$waldLRT$de.total$pvalue < 0.05,]
dim(cerv.thor.untr.de.uqua.notnorm2.subacute.lrt.wald.pval.minorthr005)##3345

cerv.thor.de.uqua.notnorm2.subacute.lrt.wald.minus.pval.005 <- setdiff(rownames(cerv.thor.de.uqua.notnorm2.subacute.lrt.wald$waldLRT$SIGN), rownames(cerv.thor.untr.de.uqua.notnorm2.subacute.lrt.wald.pval.minorthr005)) 
length(cerv.thor.de.uqua.notnorm2.subacute.lrt.wald.minus.pval.005)## 1147


subsetted.counts <- normalized.counts.prop.uqua[ which(rownames(normalized.counts.prop.uqua) %in% cerv.thor.de.uqua.notnorm2.subacute.lrt.wald.minus.pval.005), ]
subsetted.counts <- AttachConvertedColumn(subsetted.counts, ens.symb.biomart.map)

for(i in 1:dim(subsetted.counts)[1]) {
  tryCatch({
    PlotCountsAlongTimes(normalized.counts = subsetted.counts, design.matrix = all.design.file, gene.name = subsetted.counts$symbol[i], gene.name.column.name = "symbol", show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = "cervical_thoracic_whole/results_final_20/Subacute_Cerv-Thor_Tr_intersect_Untr_minus_pval_005/intersection_genes", prefix.plot = "intersection gene")
  }, error=function(e) {
    print(e)
  })
  
}


functional.folder <- "/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/subacute_Cerv-Thor_Tr_intersect_Untr_pval_005/functional/WaldLRT_Tr_intersect_Untr/"
prefix.label = "Subacute_Cerv_Thor_tr_vs_untr prop uqua WaldLRT_tr_intersect_untr_pval_005"
enrichSuitedSingleList(de.gene.list = cerv.thor.de.uqua.notnorm2.subacute.lrt.wald.minus.pval.005, functional.folder = functional.folder, filename = prefix.label)





######################## final union/intersection

cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.pval.minorthr005 <- cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01$waldLRT$de.total[cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01$waldLRT$de.total$pvalue < 0.05, ]
dim(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.pval.minorthr005)

sum(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01$LRT$de.total$padj < 0.05, na.rm = TRUE)

cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.pval.minorthr005.genes <- union(rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.pval.minorthr005), rownames(cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.thr01$LRT$SIGN))

cerv.thor.de.uqua.notnorm2.des4time.lrt.wald.union.genes <- union(rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN), rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN))

cerv.thor.de.uqua.notnorm2.des4time.lrt.wald.union.genes.minus.untr <-  setdiff(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald.union.genes, cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald.pval.minorthr005.genes)


output.results.folder <- UpdateFolderPath(output.results.folder, "final_intersection")
Venn3de(cervical.4times.union, thoracic.4times.union, cerv.thor.de.uqua.notnorm2.des4time.lrt.wald.union.genes.minus.untr, label1 = "Cervical_union", label2 = "Thoracic_union", label3 = "Cerv-Thor_tr_minus_untr", intersection.flag = TRUE, title = "DE Genes", plot.dir = output.results.folder, enrich.lists.flag = TRUE, intersection.exclusion.flag = TRUE)










