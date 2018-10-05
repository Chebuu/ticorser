## final intersections
# date <- gsub(pattern = " ", replacement = "_", date())
# date <- gsub(pattern = ":", replacement = "_", date)
# output.results.folder <- paste0("cervical_thoracic_whole/", date)


FilterOutUntreatedGenes <- function(dataset.treated, dataset.untreated, column.name = c("pvalue", "padj"), threshold = 0.05, prefix = NULL, output.dir = NULL) {
  
  prefix = UpdatePrefix(prefix, column.name, threshold)
  output.dir = UpdateFolderPath(path = output.dir, prefix)
  print(output.dir)
  
  dataset.untreated.sign <- dataset.untreated[which(dataset.untreated[,column.name] < threshold),]
  filename = UpdateFilename("untreated.sign", prefix)
  # dim(dataset.untreated.sign)
  WriteDataFrameAsTsv(data.frame.to.save = dataset.untreated.sign, file.name.path = file.path(output.dir, filename))
  
  dataset.treated.minus.untreated.rownames <- setdiff(rownames(dataset.treated), rownames(dataset.untreated.sign))
  dataset.treated.minus.untreated <- subset(dataset.treated, rownames(dataset.treated) %in% dataset.treated.minus.untreated.rownames)
  filename = UpdateFilename("treated.minus.untreated", prefix)
  WriteDataFrameAsTsv(data.frame.to.save = dataset.treated.minus.untreated, file.name.path = file.path(output.dir, filename))
  # dim(dataset.treated.minus.untreated)
  
  dataset.treated.intersect.untreated.rownames <- intersect(rownames(dataset.treated), rownames(dataset.untreated.sign))
  dataset.treated.intersect.untreated <- subset(dataset.treated, rownames(dataset.treated) %in% dataset.treated.intersect.untreated.rownames)
  filename = UpdateFilename("treated.intersect.untreated", prefix)
  WriteDataFrameAsTsv(data.frame.to.save = dataset.treated.intersect.untreated, file.name.path = file.path(output.dir, filename))
  
  
  # dim(dataset.treated.intersect.untreated)
  return(list("dataset.untreated.sign"=dataset.untreated.sign, "dataset.treated.minus.untreated"=dataset.treated.minus.untreated, "dataset.treated.intersect.untreated"=dataset.treated.intersect.untreated))
}

FilterOutUntreatedGenesAllThresholds <- function(dataset.treated, dataset.untreated, prefix, output.dir=NULL) {
  
  output.dir <- UpdateFolderPath(output.dir, prefix)
  print(output.dir)
  dataset.padj005.list <- FilterOutUntreatedGenes(dataset.treated, dataset.untreated, column.name = "padj", prefix=prefix, output.dir=output.dir)
  dataset.padj001.list <- FilterOutUntreatedGenes(dataset.treated, dataset.untreated, column.name = "padj", threshold = 0.01, prefix=prefix, output.dir=output.dir)
  dataset.pval005.list <- FilterOutUntreatedGenes(dataset.treated, dataset.untreated, column.name = "pvalue", prefix=prefix, output.dir=output.dir)
  dataset.pval001.list <- FilterOutUntreatedGenes(dataset.treated, dataset.untreated, column.name = "pvalue", threshold = 0.01, prefix=prefix, output.dir=output.dir)
  
  # dataset.untreated.pval005 <- dataset.untreated[dataset.untreated$pvalue < 0.05,]
  # dataset.untreated.pval001 <- dataset.untreated[dataset.untreated$pvalue < 0.01,]
  # dataset.untreated.padj005 <- dataset.untreated[dataset.untreated$padj < 0.05,]
  # dataset.untreated.padj001 <- dataset.untreated[dataset.untreated$padj < 0.01,]
  
  
  # dataset.treated.minus.untreted.pval005 <- setdiff(rownames(dataset.treated), rownames(dataset.untreated.pval005))
  # dataset.treated.minus.untreted.padj005 <- setdiff(rownames(dataset.treated), rownames(dataset.untreated.padj005))
  # dataset.treated.minus.untreted.pval001 <- setdiff(rownames(dataset.treated), rownames(dataset.untreated.pval001))
  # dataset.treated.minus.untreted.padj001 <- setdiff(rownames(dataset.treated), rownames(dataset.untreated.padj001))
  
  # dataset.treated.intersect.untreted.pval005 <- intersect(rownames(dataset.treated), rownames(dataset.untreated.pval005))
  # dataset.treated.intersect.untreted.padj005 <- intersect(rownames(dataset.treated), rownames(dataset.untreated.padj005))
  # dataset.treated.intersect.untreted.pval001 <- intersect(rownames(dataset.treated), rownames(dataset.untreated.pval001))
  # dataset.treated.intersect.untreted.padj001 <- intersect(rownames(dataset.treated), rownames(dataset.untreated.padj001))
  
  
  dimensions.dataframe <- cbind(
                                rep(dim(dataset.treated)[1], times=2),
                                
                                rbind(dim(dataset.pval005.list$dataset.untreated.sign)[1], dim(dataset.padj005.list$dataset.untreated.sign)[1]), #"untr_sign_005"

                                rbind(dim(dataset.pval005.list$dataset.treated.intersect.untreated)[1], dim(dataset.padj005.list$dataset.treated.intersect.untreated)[1]),#tr_int_untr_005
                                rbind(dim(dataset.pval005.list$dataset.treated.minus.untreated)[1], dim(dataset.padj005.list$dataset.treated.minus.untreated)[1]),#tr_min_untr_005
                                
                                rbind(dim(dataset.pval001.list$dataset.untreated.sign)[1], dim(dataset.padj001.list$dataset.untreated.sign)[1]),#untr_sign_001
                                rbind(dim(dataset.pval001.list$dataset.treated.intersect.untreated)[1], dim(dataset.padj001.list$dataset.treated.intersect.untreated)[1]),#tr_min_untr_001
                                rbind(dim(dataset.pval001.list$dataset.treated.minus.untreated)[1], dim(dataset.padj001.list$dataset.treated.minus.untreated)[1])#tr_min_untr_001
                                
  )
  
  colnames(dimensions.dataframe) <- c("dim_tr_sign", 
                                      "untr_sign_005", 
                                      "tr_int_untr_005",
                                      "tr_filt_untr_005", 
                                      "untr_sign_001",
                                      "tr_int_untr_001", 
                                      "tr_filt_untr_001")
  
  rownames(dimensions.dataframe) <- c("pval", "padj")
  
  print(prefix)
  print(dimensions.dataframe)
  WriteDataFrameAsTsv(data.frame.to.save = dimensions.dataframe, file.name.path = file.path(output.dir, paste0(prefix, "_filtered_tr_untr_dimensions_resume")))
  return( list("dataset.tr.minus.untr.pval005" = dataset.pval005.list$dataset.treated.minus.untreated, 
               "dataset.tr.minus.untr.padj005" = dataset.padj005.list$dataset.treated.minus.untreated, 
               # "dataset.treated.intersect.untreted.pval005" = dataset.pval005.list$dataset.treated.intersect.untreted.pval005, 
               # "dataset.treated.intersect.untreted.padj005" = dataset.padj005.list$dataset.treated.intersect.untreted.padj005, 
               "dataset.tr.minus.untr.pval001" = dataset.pval001.list$dataset.treated.minus.untreated, 
               "dataset.tr.minus.untr.padj001" = dataset.padj001.list$dataset.treated.minus.untreated, 
               # "dataset.treated.intersect.untreted.pval001" = dataset.pval001.list$dataset.treated.intersect.untreted.pval001, 
               # "dataset.treated.intersect.untreted.padj001" = dataset.padj001.list$dataset.treated.intersect.untreted.padj001, 
               "dimensions.dataframe" = dimensions.dataframe) )
}

output.filtering.folder <- UpdateFolderPath(output.results.folder, "4Points_Cerv-Thor_Wald_filtering")
# dataset.treated = cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN; dataset.untreated = cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$de.total; prefix="4Points_Cerv-Thor_Wald_filtering"; output.dir = output.results.folder
filtered.datasets.4times.wald.SIGN <- FilterOutUntreatedGenesAllThresholds(dataset.treated = cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN, dataset.untreated = cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$de.not.na, prefix="SIGN", output.dir = output.filtering.folder)
filtered.datasets.4times.wald.UP <- FilterOutUntreatedGenesAllThresholds(dataset.treated = cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$UP, dataset.untreated = cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$de.not.na, prefix="UP", output.dir = output.filtering.folder)
filtered.datasets.4times.wald.DOWN <- FilterOutUntreatedGenesAllThresholds(dataset.treated = cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$DOWN, dataset.untreated = cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$de.not.na, prefix="DOWN", output.dir = output.filtering.folder)

# 
# cerv.thor.de.uqua.notnorm2.des4time.lrt.wald.filtered.pval005.union <- union(rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN), rownames(filtered.datasets.4times.wald.SIGN$dataset.tr.minus.untr.pval005))
# cerv.thor.de.uqua.notnorm2.des4time.lrt.wald.filtered.padj005.union <- union(rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN), rownames(filtered.datasets.4times.wald.SIGN$dataset.tr.minus.untr.padj005))
# cerv.thor.de.uqua.notnorm2.des4time.lrt.wald.filtered.pval001.union <- union(rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN), rownames(filtered.datasets.4times.wald.SIGN$dataset.tr.minus.untr.pval001))
# cerv.thor.de.uqua.notnorm2.des4time.lrt.wald.filtered.padj001.union <- union(rownames(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN), rownames(filtered.datasets.4times.wald.SIGN$dataset.tr.minus.untr.padj001))


UnificateResultsAndEnrich <- function(lrt.results, wald.filtered.dataset, prefix, out.dir) {
  out.dir <- UpdateFolderPath(out.dir, prefix)
  lrt.wald.filtered.union <- union(rownames(lrt.results), rownames(wald.filtered.dataset))
  WriteDataFrameAsTsv(data.frame.to.save = as.data.frame(lrt.wald.filtered.union), file.name.path = file.path(out.dir, "lrt_wald_filtered_union_genes"), row.names = FALSE, col.names = FALSE)
  filename <- UpdatePrefix(prefix, "union")
  enrichSuitedSingleList(de.gene.list = lrt.wald.filtered.union, functional.folder = file.path(out.dir, "functional"), filename = filename)
  return(lrt.wald.filtered.union)
}

UnificateResultsAndEnrichLists <- function(lrt.results, wald.filtered.dataset.list, prefix, root.dir) {

  # new.prefix <- UpdatePrefix(prefix, "lrt_union_wald_filtered_tr_untr_pval_005")
  # lrt.wald.filtered.union.pval005 <- UnificateResultsAndEnrich(lrt.results, wald.filtered.dataset.list$dataset.tr.minus.untr.pval005, prefix=new.prefix, out.dir=file.path(root.dir, prefix))
  # new.prefix <- UpdatePrefix(prefix, "lrt_union_wald_filtered_tr_untr_pval_001")
  # lrt.wald.filtered.union.pval001 <- UnificateResultsAndEnrich(lrt.results, wald.filtered.dataset.list$dataset.tr.minus.untr.pval001, prefix=new.prefix, out.dir=file.path(root.dir, prefix))
  new.prefix <- UpdatePrefix(prefix, "lrt_union_wald_filtered_tr_untr_padj_005")
  lrt.wald.filtered.union.padj005 <- UnificateResultsAndEnrich(lrt.results, wald.filtered.dataset.list$dataset.tr.minus.untr.padj005, prefix=new.prefix, out.dir=file.path(root.dir, prefix))
  # new.prefix <- UpdatePrefix(prefix, "lrt_union_wald_filtered_tr_untr_padj_001")
  # lrt.wald.filtered.union.padj001 <- UnificateResultsAndEnrich(lrt.results, wald.filtered.dataset.list$dataset.tr.minus.untr.padj001, prefix=new.prefix, out.dir=file.path(root.dir, prefix))
  
  return(list(#"lrt.wald.filtered.union.pval005"=lrt.wald.filtered.union.pval005, #"lrt.wald.filtered.union.pval001"=lrt.wald.filtered.union.pval001, 
              "lrt.wald.filtered.union.padj005"=lrt.wald.filtered.union.padj005))#, "lrt.wald.filtered.union.padj001"=lrt.wald.filtered.union.padj001))
}

output.union.folder <- UpdateFolderPath(output.results.folder, "4Points_Cerv-Thor_Wald_union")

filtered.datasets.4times.lrt.wald.filtered.union.sign <- UnificateResultsAndEnrichLists(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN, filtered.datasets.4times.wald.SIGN, prefix = "SIGN", root.dir = output.union.folder)
filtered.datasets.4times.lrt.wald.filtered.union.up <- UnificateResultsAndEnrichLists(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$UP, filtered.datasets.4times.wald.UP, prefix = "UP", root.dir = output.union.folder)
filtered.datasets.4times.lrt.wald.filtered.union.down <- UnificateResultsAndEnrichLists(cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$LRT$DOWN, filtered.datasets.4times.wald.DOWN, prefix = "DOWN", root.dir = output.union.folder)

filtered.datasets.4times.lrt.wald.filtered.union <- list()
filtered.datasets.4times.lrt.wald.filtered.union$SIGN <- filtered.datasets.4times.lrt.wald.filtered.union.sign
filtered.datasets.4times.lrt.wald.filtered.union$UP <- filtered.datasets.4times.lrt.wald.filtered.union.up
filtered.datasets.4times.lrt.wald.filtered.union$DOWN <- filtered.datasets.4times.lrt.wald.filtered.union.down


# output.intersections.folder <- UpdateFolderPath(output.results.folder, "4Points_Itersections")

TimeCourseIntersectAndEnrichUnions <- function(cerv.df, thor.df, cerv.thor.filtered.list, main.title, intersection.flag=TRUE, intersection.exclusion.flag=FALSE, enrich.union.intersection=TRUE, conversion.map=NULL, root.dir=NUL) {
  
  for(lbl in c("SIGN", "UP", "DOWN")) {
    cerv.union <- union(rownames(cerv.df$LRT[[lbl]]), rownames(cerv.df$waldLRT[[lbl]]))
    thor.union <- union(rownames(thor.df$LRT[[lbl]]), rownames(thor.df$waldLRT[[lbl]]))
    cerv.thor.union.list <- cerv.thor.filtered.list[[lbl]]
    
    for(cerv.thor.filtered.ind in names(cerv.thor.union.list)){
      print(cerv.thor.filtered.ind)
      if(cerv.thor.filtered.ind == "lrt.wald.filtered.union.padj005") {
        cerv.thor.filtered.union <- cerv.thor.union.list[[cerv.thor.filtered.ind]]
        plot.dir <- UpdateFolderPath(root.dir, lbl, cerv.thor.filtered.ind)
        ct.label <- paste("Cerv-Thor", gsub(pattern = ".", replacement = " ", cerv.thor.filtered.ind, fixed = TRUE))
        
        Venn3de(x = cerv.union, y = thor.union, z = cerv.thor.filtered.union, label1 = "Cervical_Union", label2 = "Thoracic_Union", label3 = ct.label, title = paste(main.title, lbl), 
                intersection.flag = intersection.flag, intersection.exclusion.flag = intersection.exclusion.flag, 
                conversion.map=conversion.map, plot.dir = plot.dir, enrich.lists.flag = TRUE)
      }
    }
  }
}


output.intersections.folder <- UpdateFolderPath(output.results.folder, "4Points_Itersections")

TimeCourseIntersectAndEnrichUnions(cerv.df = cerv.de.uqua.notnorm2.des4time.lrt.wald, thor.df = thor.de.uqua.notnorm2.des4time.lrt.wald, cerv.thor.filtered.list = filtered.datasets.4times.lrt.wald.filtered.union, root.dir=output.intersections.folder, main.title="4Points Intersections", conversion.map = ens.symb.biomart.map)


# output.filtering.folder <- UpdateFolderPath(output.results.folder, "Subacute_Cerv-Thor_Wald_filtering")
# # dataset.treated = cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN; dataset.untreated = cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$de.total; prefix="4Points_Cerv-Thor_Wald_filtering"; output.dir = output.results.folder
# filtered.datasets.4times.wald.SIGN <- FilterOutUntreatedGenesAllThresholds(dataset.treated = cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN, dataset.untreated = cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$de.not.na, prefix="SIGN", output.dir = output.filtering.folder)
# filtered.datasets.4times.wald.UP <- FilterOutUntreatedGenesAllThresholds(dataset.treated = cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$UP, dataset.untreated = cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$de.not.na, prefix="UP", output.dir = output.filtering.folder)
# filtered.datasets.4times.wald.DOWN <- FilterOutUntreatedGenesAllThresholds(dataset.treated = cerv.thor.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$DOWN, dataset.untreated = cerv.thor.untr.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$de.not.na, prefix="DOWN", output.dir = output.filtering.folder)





















