### intersections 03d

# 
# FilterOutUntreatedGenes <- function(dataset.treated, dataset.untreated, prefix, output.dir=NULL) {
#   dataset.untreated.pval <- dataset.untreated[dataset.untreated$pvalue < 0.05,]
#   
#   dataset.untreated.padj <- dataset.untreated[dataset.untreated$padj < 0.05,]
#   
#   dataset.treated.minus.untreted.pval <- setdiff(rownames(dataset.treated), rownames(dataset.untreated.pval))
#   dataset.treated.minus.untreted.padj <- setdiff(rownames(dataset.treated), rownames(dataset.untreated.padj))
#   
#   dataset.treated.intersect.untreted.pval <- intersect(rownames(dataset.treated), rownames(dataset.untreated.pval))
#   dataset.treated.intersect.untreted.padj <- intersect(rownames(dataset.treated), rownames(dataset.untreated.padj))
#   
#   
#   dimensions.dataframe <- cbind(rep(dim(dataset.treated)[1], times=2), 
#                                 rbind(dim(dataset.untreated.pval)[1],  dim(dataset.untreated.padj)[1]), 
#                                 rbind(length(dataset.treated.minus.untreted.pval), length(dataset.treated.minus.untreted.padj)), 
#                                 rbind(length(dataset.treated.intersect.untreted.pval), length(dataset.treated.intersect.untreted.padj))
#                           )
#   
#   colnames(dimensions.dataframe) <- c("dim_tr_sign", "dim_untr_sign", "dim_res_filt", "dim_filt_out")
#   rownames(dimensions.dataframe) <- c("pval", "padj")
#   
#   print(prefix)
#   print(dimensions.dataframe)
#   return( list("dataset.tr.minus.untr.pval"= dataset.treated.minus.untreted.pval, "dataset.tr.minus.untr.padj"= dataset.treated.minus.untreted.padj, "dataset.treated.intersect.untreted.pval"=dataset.treated.intersect.untreted.pval, "dataset.treated.intersect.untreted.padj"=dataset.treated.intersect.untreted.padj, "dimensions.dataframe"=dimensions.dataframe) )
# }

# 
# filtered.datasets.03d.wald.SIGN <- FilterOutUntreatedGenesAllThresholds(dataset.treated = cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN, dataset.untreated = cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$de.not.na, prefix="SIGN", 
#                                                                         output.dir = output.filtering.folder)
# 
# 
# cerv.thor.03.filtered.sign <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.03d.deseq$SIGN, dataset.untreated = cerv.thor.de.uqua.notnorm2.03d.deseq.untr$de.not.na, prefix = "03d SIGN")
# cerv.thor.07.filtered.sign <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.07d.deseq$SIGN, dataset.untreated = cerv.thor.de.uqua.notnorm2.07d.deseq.untr$de.not.na, prefix = "07d SIGN")
# cerv.thor.14.filtered.sign <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.14d.deseq$SIGN, dataset.untreated = cerv.thor.de.uqua.notnorm2.14d.deseq.untr$de.not.na, prefix = "14d SIGN")
# cerv.thor.56.filtered.sign <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.56d.deseq$SIGN, dataset.untreated = cerv.thor.de.uqua.notnorm2.56d.deseq.untr$de.not.na, prefix = "56d SIGN")
# 
# 
# cerv.thor.03.filtered.up <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.03d.deseq$UP, dataset.untreated = cerv.thor.de.uqua.notnorm2.03d.deseq.untr$de.not.na, prefix = "03d up")
# cerv.thor.07.filtered.up <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.07d.deseq$UP, dataset.untreated = cerv.thor.de.uqua.notnorm2.07d.deseq.untr$de.not.na, prefix = "07d up")
# cerv.thor.14.filtered.up <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.14d.deseq$UP, dataset.untreated = cerv.thor.de.uqua.notnorm2.14d.deseq.untr$de.not.na, prefix = "14d up")
# cerv.thor.56.filtered.up <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.56d.deseq$UP, dataset.untreated = cerv.thor.de.uqua.notnorm2.56d.deseq.untr$de.not.na, prefix = "56d up")
# 
# 
# cerv.thor.03.filtered.down <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.03d.deseq$DOWN, dataset.untreated = cerv.thor.de.uqua.notnorm2.03d.deseq.untr$de.not.na, prefix = "03d down")
# cerv.thor.07.filtered.down <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.07d.deseq$DOWN, dataset.untreated = cerv.thor.de.uqua.notnorm2.07d.deseq.untr$de.not.na, prefix = "07d down")
# cerv.thor.14.filtered.down <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.14d.deseq$DOWN, dataset.untreated = cerv.thor.de.uqua.notnorm2.14d.deseq.untr$de.not.na, prefix = "14d down")
# cerv.thor.56.filtered.down <- FilterOutUntreatedGenes(dataset.treated = cerv.thor.de.uqua.notnorm2.56d.deseq$DOWN, dataset.untreated = cerv.thor.de.uqua.notnorm2.56d.deseq.untr$de.not.na, prefix = "56d down")




FilterOutUntreatedGenesSingleTimePoint <- function(dataset.treated.list, dataset.untreated, prefix, root.dir) {
  filtering.folder <- UpdateFolderPath(root.dir, "Cerv-Thor_Tr_filtering_Untr")
  treated.filtered.list <- list()
  for(lbl in c("SIGN", "UP", "DOWN")) {
    treated.filtered.list[[lbl]] <- FilterOutUntreatedGenesAllThresholds(dataset.treated = dataset.treated.list[[lbl]], dataset.untreated = dataset.untreated$de.not.na, prefix = paste(prefix, lbl), output.dir=filtering.folder)
    
  }
  return(treated.filtered.list)
}


timepoints.results.folder <- UpdateFolderPath(output.results.folder, "Time_by_Time")

output.folder.03 <- UpdateFolderPath(timepoints.results.folder, "03d")
cerv.thor.03.filtered.list <- FilterOutUntreatedGenesSingleTimePoint(dataset.treated.list = cerv.thor.de.uqua.notnorm2.03d.deseq, dataset.untreated = cerv.thor.de.uqua.notnorm2.03d.deseq.untr, prefix = "03d", root.dir=output.folder.03)
output.folder.07 <- UpdateFolderPath(timepoints.results.folder, "07d")
cerv.thor.07.filtered.list <- FilterOutUntreatedGenesSingleTimePoint(dataset.treated.list = cerv.thor.de.uqua.notnorm2.07d.deseq, dataset.untreated = cerv.thor.de.uqua.notnorm2.07d.deseq.untr, prefix = "07d", root.dir=output.folder.07)
output.folder.14 <- UpdateFolderPath(timepoints.results.folder, "14d")
cerv.thor.14.filtered.list <- FilterOutUntreatedGenesSingleTimePoint(dataset.treated.list = cerv.thor.de.uqua.notnorm2.14d.deseq, dataset.untreated = cerv.thor.de.uqua.notnorm2.14d.deseq.untr, prefix = "14d", root.dir=output.folder.14)
output.folder.56 <- UpdateFolderPath(timepoints.results.folder, "56d")
cerv.thor.56.filtered.list <- FilterOutUntreatedGenesSingleTimePoint(dataset.treated.list = cerv.thor.de.uqua.notnorm2.56d.deseq, dataset.untreated = cerv.thor.de.uqua.notnorm2.56d.deseq.untr, prefix = "56d", root.dir=output.folder.56)


SingleTimePointIntersectAndEnrich <- function(cerv.df, thor.df, cerv.thor.filtered.list, main.title, intersection.flag=TRUE, intersection.exclusion.flag=FALSE, enrich.union.intersection=TRUE, conversion.map=NULL, root.dir=NUL) {
  output.dir <- UpdateFolderPath(root.dir, "Intersections")
  for(lbl in c("SIGN", "UP", "DOWN")) {
    cerv <- rownames(cerv.df[[lbl]])
    thor <- rownames(thor.df[[lbl]])
    cerv.thor.list <- cerv.thor.filtered.list[[lbl]]
    
    for(cerv.thor.filtered.ind in names(cerv.thor.list)){
      
      # if(cerv.thor.filtered.ind!= "dimensions.dataframe") {
      if(cerv.thor.filtered.ind == "dataset.tr.minus.untr.padj005") {
        cerv.thor.filtered <- rownames(cerv.thor.list[[cerv.thor.filtered.ind]])
        plot.dir <- UpdateFolderPath(output.dir, lbl, cerv.thor.filtered.ind)
        ct.label <- paste("Cerv-Thor", lbl, gsub(pattern = ".", replacement = " ", cerv.thor.filtered.ind, fixed = TRUE))

        Venn3de(x = cerv, y = thor, z = cerv.thor.filtered, label1 = paste("Cervical", lbl), label2 = paste("Thoracic", lbl), label3 = ct.label, title = paste(main.title, lbl), 
                intersection.flag = intersection.flag, intersection.exclusion.flag = intersection.exclusion.flag, 
                conversion.map=conversion.map, plot.dir = plot.dir, enrich.lists.flag = TRUE)
      }
    }
  }
}

SingleTimePointIntersectAndEnrich(cerv.df = cervical.de.uqua.notnorm2.03d.deseq, thor.df = thoracic.de.uqua.notnorm2.03d, cerv.thor.filtered.list = cerv.thor.03.filtered.list, root.dir=output.folder.03, main.title="03d Intersections", conversion.map = ens.symb.biomart.map)
SingleTimePointIntersectAndEnrich(cerv.df = cervical.de.uqua.notnorm2.07d.deseq, thor.df = thoracic.de.uqua.notnorm2.07d.deseq, cerv.thor.filtered.list = cerv.thor.07.filtered.list, root.dir=output.folder.07, main.title="07d Intersections", conversion.map = ens.symb.biomart.map)
SingleTimePointIntersectAndEnrich(cerv.df = cervical.de.uqua.notnorm2.14d.deseq, thor.df = thoracic.de.uqua.notnorm2.14d.deseq, cerv.thor.filtered.list = cerv.thor.14.filtered.list, root.dir=output.folder.14, main.title="14d Intersections", conversion.map = ens.symb.biomart.map)
SingleTimePointIntersectAndEnrich(cerv.df = cervical.de.uqua.notnorm2.56d.deseq, thor.df = thoracic.de.uqua.notnorm2.56d.deseq, cerv.thor.filtered.list = cerv.thor.56.filtered.list, root.dir=output.folder.56, main.title="56d Intersections", conversion.map = ens.symb.biomart.map)



SingleTimePointIntersectMixAndEnrich <- function(cerv.df, thor.df, cerv.thor.filtered.list, main.title, intersection.flag=TRUE, intersection.exclusion.flag=FALSE, enrich.union.intersection=TRUE, conversion.map=NULL, root.dir=NUL) {
  output.dir <- UpdateFolderPath(root.dir, "Intersections_mix")
  cerv <- rownames(cerv.df[["SIGN"]])
  thor <- rownames(thor.df[["SIGN"]])
  
  for(lbl in c("UP", "DOWN")) {
    
    cerv.thor.list <- cerv.thor.filtered.list[[lbl]]
    
    for(cerv.thor.filtered.ind in names(cerv.thor.list)){
      
      # if(cerv.thor.filtered.ind!= "dimensions.dataframe") {
      if(cerv.thor.filtered.ind == "dataset.tr.minus.untr.padj005") {
        cerv.thor.filtered <- rownames(cerv.thor.list[[cerv.thor.filtered.ind]])
        plot.dir <- UpdateFolderPath(output.dir, lbl, cerv.thor.filtered.ind)
        ct.label <- paste("Cerv-Thor", lbl, gsub(pattern = ".", replacement = " ", cerv.thor.filtered.ind, fixed = TRUE))
        
        Venn3de(x = cerv, y = thor, z = cerv.thor.filtered, label1 = "Cervical", label2 = "Thoracic", label3 = ct.label, title = main.title, 
                intersection.flag = intersection.flag, intersection.exclusion.flag = intersection.exclusion.flag, 
                conversion.map=conversion.map, plot.dir = plot.dir, enrich.lists.flag = TRUE)
      }
    }
  }
}

SingleTimePointIntersectMixAndEnrich(cerv.df = cervical.de.uqua.notnorm2.03d.deseq, thor.df = thoracic.de.uqua.notnorm2.03d, cerv.thor.filtered.list = cerv.thor.03.filtered.list, root.dir=output.folder.03, main.title="03d Intersections", conversion.map = ens.symb.biomart.map)
SingleTimePointIntersectMixAndEnrich(cerv.df = cervical.de.uqua.notnorm2.07d.deseq, thor.df = thoracic.de.uqua.notnorm2.07d.deseq, cerv.thor.filtered.list = cerv.thor.07.filtered.list, root.dir=output.folder.07, main.title="07d Intersections", conversion.map = ens.symb.biomart.map)
SingleTimePointIntersectMixAndEnrich(cerv.df = cervical.de.uqua.notnorm2.14d.deseq, thor.df = thoracic.de.uqua.notnorm2.14d.deseq, cerv.thor.filtered.list = cerv.thor.14.filtered.list, root.dir=output.folder.14, main.title="14d Intersections", conversion.map = ens.symb.biomart.map)
SingleTimePointIntersectMixAndEnrich(cerv.df = cervical.de.uqua.notnorm2.56d.deseq, thor.df = thoracic.de.uqua.notnorm2.56d.deseq, cerv.thor.filtered.list = cerv.thor.56.filtered.list, root.dir=output.folder.56, main.title="56d Intersections", conversion.map = ens.symb.biomart.map)


# 
# dim(cerv.thor.de.uqua.notnorm2.03d.deseq$SIGN)
# dim(cerv.thor.de.uqua.notnorm2.03d.deseq.untr$SIGN)
# cerv.thor.de.uqua.notnorm2.03d.deseq.untr.005 <- cerv.thor.de.uqua.notnorm2.03d.deseq.untr$de.total[cerv.thor.de.uqua.notnorm2.03d.deseq.untr$de.total$pvalue < 0.05,]
# cerv.thor.de.uqua.notnorm2.03d.deseq.untr.005 <- cerv.thor.de.uqua.notnorm2.03d.deseq.untr$de.total[cerv.thor.de.uqua.notnorm2.03d.deseq.untr$de.not.na$pvalue < 0.05,]
# dim(cerv.thor.de.uqua.notnorm2.03d.deseq.untr.005)
# 
# cerv.thor.de.uqua.notnorm2.03d.deseq.sign.filtuntr <- cerv.thor.de.uqua.notnorm2.03d.deseq$SIGN[-which(rownames(cerv.thor.de.uqua.notnorm2.03d.deseq$SIGN) %in% rownames(cerv.thor.de.uqua.notnorm2.03d.deseq.untr.005)),]
# cerv.thor.de.uqua.notnorm2.03d.deseq.sign.filtuntr <- setdiff(x = rownames(cerv.thor.de.uqua.notnorm2.03d.deseq$SIGN), y = rownames(cerv.thor.de.uqua.notnorm2.03d.deseq.untr.005) )
# dim(cerv.thor.de.uqua.notnorm2.03d.deseq.sign.filtuntr)
# 
# cerv.thor.de.uqua.notnorm2.03d.deseq.up.filtuntr <- cerv.thor.de.uqua.notnorm2.03d.deseq$UP[-which(rownames(cerv.thor.de.uqua.notnorm2.03d.deseq$UP) %in% rownames(cerv.thor.de.uqua.notnorm2.03d.deseq.untr.005)),]
# dim(cerv.thor.de.uqua.notnorm2.03d.deseq$UP)
# dim(cerv.thor.de.uqua.notnorm2.03d.deseq.up.filtuntr)
# 
# cerv.thor.de.uqua.notnorm2.03d.deseq.down.filtuntr <- cerv.thor.de.uqua.notnorm2.03d.deseq$DOWN[-which(rownames(cerv.thor.de.uqua.notnorm2.03d.deseq$DOWN) %in% rownames(cerv.thor.de.uqua.notnorm2.03d.deseq.untr.005)),]
# dim(cerv.thor.de.uqua.notnorm2.03d.deseq$DOWN)
# dim(cerv.thor.de.uqua.notnorm2.03d.deseq.down.filtuntr)
# 
# cerv.thor.de.uqua.notnorm2.03d.deseq.filtuntr <- list("SIGN"= cerv.thor.de.uqua.notnorm2.03d.deseq.sign.filtuntr, "UP"=cerv.thor.de.uqua.notnorm2.03d.deseq.up.filtuntr, "DOWN"= cerv.thor.de.uqua.notnorm2.03d.deseq.down.filtuntr)
# 
# VennDE3ListsSuitedCervThorCT(cervical.de.uqua.notnorm2.03d.deseq, thoracic.de.uqua.notnorm2.03d, cerv.thor.de.uqua.notnorm2.03d.deseq.filtuntr, "cerv_thor 03d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/03d/Intersections")
# 
# 
# 
# FilterSingleTimePointByUntreatedResults <- function(list.datasets.to.filter, list.dataset.untreated, pvalue.threshold=0.05) {
#   
#   dataset.untreated.all.sign005 <- list.dataset.untreated$de.not.na[list.dataset.untreated$de.not.na$pvalue < pvalue.threshold,]
#   sign.filtuntr <- list.datasets.to.filter$SIGN[-which(rownames(list.datasets.to.filter$SIGN) %in% rownames(dataset.untreated.all.sign005)),]
#   up.filtuntr <- list.datasets.to.filter$UP[-which(rownames(list.datasets.to.filter$UP) %in% rownames(dataset.untreated.all.sign005)),]
#   down.filtuntr <- list.datasets.to.filter$DOWN[-which(rownames(list.datasets.to.filter$DOWN) %in% rownames(dataset.untreated.all.sign005)),]
#   
#   list.datasets.filtuntr <- list("SIGN"= sign.filtuntr, "UP"=up.filtuntr, "DOWN"= down.filtuntr, "untr.sign"= dataset.untreated.all.sign005)
#   return(list.datasets.filtuntr)
# }
# 
# cerv.thor.de.uqua.notnorm2.03d.deseq.filtuntr <- FilterSingleTimePointByUntreatedResults(list.datasets.to.filter=cerv.thor.de.uqua.notnorm2.03d.deseq, list.dataset.untreated=cerv.thor.de.uqua.notnorm2.03d.deseq.untr)
# cerv.thor.de.uqua.notnorm2.07d.deseq.filtuntr <- FilterSingleTimePointByUntreatedResults(list.datasets.to.filter=cerv.thor.de.uqua.notnorm2.07d.deseq, list.dataset.untreated=cerv.thor.de.uqua.notnorm2.07d.deseq.untr)
# cerv.thor.de.uqua.notnorm2.14d.deseq.filtuntr <- FilterSingleTimePointByUntreatedResults(list.datasets.to.filter=cerv.thor.de.uqua.notnorm2.14d.deseq, list.dataset.untreated=cerv.thor.de.uqua.notnorm2.14d.deseq.untr)
# cerv.thor.de.uqua.notnorm2.56d.deseq.filtuntr <- FilterSingleTimePointByUntreatedResults(list.datasets.to.filter=cerv.thor.de.uqua.notnorm2.56d.deseq, list.dataset.untreated=cerv.thor.de.uqua.notnorm2.56d.deseq.untr)
#  
# VennDE3ListsSuitedCervThorCT(cervical.de.uqua.notnorm2.03d.deseq, thoracic.de.uqua.notnorm2.03d, cerv.thor.de.uqua.notnorm2.03d.deseq.filtuntr, "cerv_thor 03d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/03d/Intersections")
# VennDE3ListsSuitedCervThorCT(cervical.de.uqua.notnorm2.07d.deseq, thoracic.de.uqua.notnorm2.07d.deseq, cerv.thor.de.uqua.notnorm2.07d.deseq.filtuntr, "cerv_thor 07d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/07d/Intersections")
# VennDE3ListsSuitedCervThorCT(cervical.de.uqua.notnorm2.14d.deseq, thoracic.de.uqua.notnorm2.14d.deseq, cerv.thor.de.uqua.notnorm2.14d.deseq.filtuntr, "cerv_thor 14d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/14d/Intersections")
# VennDE3ListsSuitedCervThorCT(cervical.de.uqua.notnorm2.56d.deseq, thoracic.de.uqua.notnorm2.56d.deseq, cerv.thor.de.uqua.notnorm2.56d.deseq.filtuntr, "cerv_thor 56d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/56d/Intersections", intersection.flag=TRUE, intersection.exclusion.flag=TRUE)
# graphics.off()
# VennDE3ListsSuitedCervThorCT(cervical.de.uqua.notnorm2.03d.deseq, thoracic.de.uqua.notnorm2.03d, cerv.thor.de.uqua.notnorm2.03d.deseq.filtuntr, "cerv_thor 03d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/03d/Intersections", intersection.flag=TRUE, intersection.exclusion.flag=TRUE)
# graphics.off()
# VennDE3ListsSuitedCervThorCT(cervical.de.uqua.notnorm2.07d.deseq, thoracic.de.uqua.notnorm2.07d.deseq, cerv.thor.de.uqua.notnorm2.07d.deseq.filtuntr, "cerv_thor 07d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/07d/Intersections", intersection.flag=TRUE, intersection.exclusion.flag=TRUE)
# graphics.off()
# VennDE3ListsSuitedCervThorCT(cervical.de.uqua.notnorm2.14d.deseq, thoracic.de.uqua.notnorm2.14d.deseq, cerv.thor.de.uqua.notnorm2.14d.deseq.filtuntr, "cerv_thor 14d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/14d/Intersections", intersection.flag=TRUE, intersection.exclusion.flag=TRUE)
# graphics.off()
# VennDE3ListsSuitedCervThorCT(cervical.de.uqua.notnorm2.56d.deseq, thoracic.de.uqua.notnorm2.56d.deseq, cerv.thor.de.uqua.notnorm2.56d.deseq.filtuntr, "cerv_thor 56d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/56d/Intersections", intersection.flag=TRUE, intersection.exclusion.flag=TRUE)
# 

VennDE3ListsSuitedCervThorMix <- function(de.df.list1, de.df.list2, de.df.list3, title, venn.plot.dir, intersection.flag=TRUE, intersection.exclusion.flag=FALSE, enrich.lists.flag=TRUE) { 
  
  label1.p = "Cervical"
  label2.p = "Thoracic"
  title.p = paste(title, "Significant DE Genes") 
  
  for (df.label in c("UP", "DOWN")) {
    switch(df.label,
           UP = { 
             title.p = paste(title.p, "and UP cerv-thor") 
             label3.p = paste("Cerv-Thor_filtered", "UP", sep = "_")
           },
           DOWN = { 
             title.p = paste(title.p, " and DOWN cerv-thor") 
             label3.p = paste("Cerv-Thor_filtered", "DOWN", sep = "_")
           }
    )
    
    Venn3de(x = rownames(de.df.list1[["SIGN"]]), y = rownames(de.df.list2[["SIGN"]]), z = rownames(de.df.list3[[df.label]]), 
            label1 = label1.p, label2 = label2.p, label3 = label3.p, title = title.p, 
            intersection.flag = intersection.flag, intersection.exclusion.flag = intersection.exclusion.flag, plot.dir = venn.plot.dir, enrich.lists.flag = TRUE)
    graphics.off()
  }
}

graphics.off()
VennDE3ListsSuitedCervThorMix(cervical.de.uqua.notnorm2.03d.deseq, thoracic.de.uqua.notnorm2.03d, cerv.thor.de.uqua.notnorm2.03d.deseq.filtuntr, "cerv_thor 03d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/03d/Intersections_mix", intersection.flag=TRUE, intersection.exclusion.flag=TRUE)
graphics.off()
VennDE3ListsSuitedCervThorMix(cervical.de.uqua.notnorm2.07d.deseq, thoracic.de.uqua.notnorm2.07d.deseq, cerv.thor.de.uqua.notnorm2.07d.deseq.filtuntr, "cerv_thor 07d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/07d/Intersections_mix", intersection.flag=TRUE, intersection.exclusion.flag=TRUE)
graphics.off()
VennDE3ListsSuitedCervThorMix(cervical.de.uqua.notnorm2.14d.deseq, thoracic.de.uqua.notnorm2.14d.deseq, cerv.thor.de.uqua.notnorm2.14d.deseq.filtuntr, "cerv_thor 14d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/14d/Intersections_mix", intersection.flag=TRUE, intersection.exclusion.flag=TRUE)
graphics.off()
VennDE3ListsSuitedCervThorMix(cervical.de.uqua.notnorm2.56d.deseq, thoracic.de.uqua.notnorm2.56d.deseq, cerv.thor.de.uqua.notnorm2.56d.deseq.filtuntr, "cerv_thor 56d", venn.plot.dir = "/media/dario/dati/time_course/cervical_thoracic_whole/results_intersections_untr_enriched/Time_by_Time/56d/Intersections_mix", intersection.flag=TRUE, intersection.exclusion.flag=TRUE)
