source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/GlobalVariables.R")
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/PlotFunctions.R")
# source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/BoxplotFunctions.R")
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/TimeCourseFunctions.R")
# source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/MainFunction.R")
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/venn3.R")
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/venn2.R")
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/FunctionalFunctions.R")
source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/HeatmapFunctions.R")
# source("/home/dario/Dropbox/Lavori/IAC/coding/time_course/src/ProjectR6Class.R")
# setwd(dir = "/media/dario/dati/time_course/")
require("org.Rn.eg.db")
# 
# FilterAndNormalizeCounts <- function(counts.dataframe, design.matrix, filter.method="Proportion", normalization.method= "uqua", is.time.course=TRUE, estimated.genes=NULL, ellipse.in.pca=FALSE) {
#   
#   # output.counts.folder <- UpdateFolderPath(output.results.folder, "counts", "FeatureCounts")
#   # output.plots.folder <- UpdateFolderPath(output.results.folder, "plots")
#   
#   # filename <- prefix.output.file
#   
#   message("Filtering counts with ", filter.method, " method\n")
#   filtered.counts <- FilterLowCounts(counts.dataframe = counts.dataframe, design.dataframe = design.matrix, is.normalized = FALSE, method.type = filter.method, cpm.cutoff = 1, cv.percentage = 1)
#   
#   # filename <- UpdateFilename(filename, filter.method)
#   # output.counts.folder <- UpdateFolderPath(output.counts.folder, filter.method)
#   # output.plots.folder <- UpdateFolderPath(output.plots.folder, filter.method)
#   
#   WriteDataFrameAsTsv(filtered.counts, file.name.path = file.path(output.counts.folder, filename))
#   
#   PlotTimesBoxplot(data.frame.to.plot = filtered.counts, design.matrix = design.matrix, output.path = output.plots.folder, prefix.plot = UpdatePrefix(prefix.plot.label, filter.method, "un-normalized"), show.plot.flag = FALSE, save.plot = TRUE, plotly.flag = TRUE)
#   
#   PlotPCAPlotlyFunction(counts.data.frame = filtered.counts, design.matrix = design.matrix, scale = FALSE, plot.folder = output.plots.folder, prefix.plot = UpdatePrefix(prefix.plot.label, filter.method, "un-normalized"), show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, ellipse.flag = ellipse.in.pca)
#   
#   message("Normalizing counts with ", normalization.method, "\n")
#   normalized.filtered.counts <- NormalizeData(filtered.counts, norm.type = normalization.method, estimated.genes = estimated.genes, design.matrix = design.matrix)
#   # print(colnames(normalized.filtered.counts))
#   
#   # filename <- UpdateFilename(filename, normalization.method)
#   # output.counts.folder <- UpdateFolderPath(output.counts.folder, normalization.method)
#   # output.plots.folder <- UpdateFolderPath(output.plots.folder, normalization.method)
#   
#   WriteDataFrameAsTsv(normalized.filtered.counts, file.name.path = file.path(output.counts.folder, filename))
#   
#   PlotTimesBoxplot(data.frame.to.plot = normalized.filtered.counts, design.matrix = design.matrix, output.path = output.plots.folder, prefix.plot = UpdatePrefix(prefix.plot.label, filter.method, normalization.method), show.plot.flag = FALSE, save.plot = TRUE, plotly.flag = TRUE)
#   
#   PlotPCAPlotlyFunction(counts.data.frame = normalized.filtered.counts, design.matrix = design.matrix, scale = FALSE, plot.folder = output.plots.folder, prefix.plot = UpdatePrefix(prefix.plot.label, filter.method, normalization.method), show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, ellipse.flag = ellipse.in.pca)
#   
#   return(list("filtered.counts.df" = filtered.counts, "normalized.counts.df"=normalized.filtered.counts, "file.name"=filename))
# }


FilterAndNormalizeCounts <- function(counts.dataframe, design.matrix, 
                                     output.results.folder, prefix.output.file, 
                                     prefix.plot.label, 
                                     filter.method="Proportion", 
                                     normalization.method="uqua", 
                                     is.time.course=TRUE, estimated.genes=NULL, 
                                     ellipse.in.pca=FALSE) {

  output.counts.folder <- UpdateFolderPath(output.results.folder, "counts", "FeatureCounts")
  output.plots.folder <- UpdateFolderPath(output.results.folder, "plots")

  filename <- prefix.output.file

  message("Filtering counts with ", filter.method, " method\n")
  filtered.counts <- FilterLowCounts(counts.dataframe = counts.dataframe, design.dataframe = design.matrix, is.normalized = FALSE, method.type = filter.method, cpm.cutoff = 1, cv.percentage = 1)

  filename <- UpdateFilename(filename, filter.method)
  output.counts.folder <- UpdateFolderPath(output.counts.folder, filter.method)
  output.plots.folder <- UpdateFolderPath(output.plots.folder, filter.method)

  WriteDataFrameAsTsv(filtered.counts, file.name.path = file.path(output.counts.folder, filename))

  # PlotTimesBoxplot(data.frame.to.plot = filtered.counts, design.matrix = design.matrix, output.path = output.plots.folder, prefix.plot = UpdatePrefix(prefix.plot.label, filter.method, "un-normalized"), show.plot.flag = FALSE, save.plot = TRUE, plotly.flag = TRUE)
  # data.frame.to.plot = filtered.counts; design.matrix = design.matrix; output.path = output.plots.folder; prefix.plot = UpdatePrefix(prefix.plot.label, filter.method, "un-normalized"); show.plot.flag = FALSE; save.plot = TRUE; plotly.flag = TRUE
  PlotPCAPlotlyFunction(counts.data.frame = filtered.counts, design.matrix = design.matrix, scale = FALSE, plot.folder = output.plots.folder, prefix.plot = UpdatePrefix(prefix.plot.label, filter.method, "un-normalized"), show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, ellipse.flag = ellipse.in.pca)

  message("Normalizing counts with ", normalization.method, "\n")
  normalized.filtered.counts <- NormalizeData(filtered.counts, 
                                              norm.type = normalization.method, estimated.genes = estimated.genes, design.matrix = design.matrix)
  # print(colnames(normalized.filtered.counts))

  filename <- UpdateFilename(filename, normalization.method)
  output.counts.folder <- UpdateFolderPath(output.counts.folder, normalization.method)
  output.plots.folder <- UpdateFolderPath(output.plots.folder, normalization.method)

  WriteDataFrameAsTsv(normalized.filtered.counts, file.name.path = file.path(output.counts.folder, filename))

  # PlotTimesBoxplot(data.frame.to.plot = normalized.filtered.counts, design.matrix = design.matrix, output.path = output.plots.folder, prefix.plot = UpdatePrefix(prefix.plot.label, filter.method, normalization.method), show.plot.flag = FALSE, save.plot = TRUE, plotly.flag = TRUE)

  PlotPCAPlotlyFunction(counts.data.frame = normalized.filtered.counts, design.matrix = design.matrix, scale = FALSE, plot.folder = output.plots.folder, prefix.plot = UpdatePrefix(prefix.plot.label, filter.method, normalization.method), show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, ellipse.flag = ellipse.in.pca)

  return(list("filtered.counts.df" = filtered.counts, "normalized.counts.df"=normalized.filtered.counts, "file.name"=filename))
}

########################
PerformDEAnalysis <- function(whole.counts, design.matrix, de.test=c("DeSeqTime_TC", "DeSeqTime_T", "DeSeqTime_NoInteraction", "DeSeq", "NOISeqBio", "NOISeq"), results.folder=NULL, prefix.label = NULL, normalize.data.flag=FALSE, normalization.method="uqua", negative.genes.list=NULL, threshold=0.05, conversion.map=NULL, enrich.results.flag=FALSE) {
  #check input data
  sub.norm.counts <- whole.counts[, which(colnames(whole.counts) %in% rownames(design.matrix))]
  
  
  plot.folder <- UpdateFolderPath(results.folder, "plots")
  de.folder <- UpdateFolderPath(results.folder, "DE_results")
  table.folder <- UpdateFolderPath(results.folder, "tables")
  functional.folder <- UpdateFolderPath(results.folder, "functional")
    
  filename <- prefix.label
  # 
  # if(normalize.data.flag) {
  #   output.counts.folder <- UpdateFolderPath(results.folder, "counts", "FeatureCounts")
  #   message("Normalizing counts with ", normalization.method, "\n")
  #   normalized.counts <- NormalizeData(sub.norm.counts, norm.type = normalization.method, design.matrix = design.matrix, estimated.genes = negative.genes.list)
  #   # print(colnames(normalized.filtered.counts))
  #   filename <- UpdateFilename(filename, normalization.method)
  #   output.counts.folder <- UpdateFolderPath(output.counts.folder, normalization.method)
  #   # plot.folder <- UpdateFolderPath(plot.folder, normalization.method)
  #   WriteDataFrameAsTsv(normalized.counts, file.name.path = file.path(output.counts.folder, filename))
  #   PlotTimesBoxplot(data.frame.to.plot = normalized.counts, design.matrix = design.matrix, output.path = plot.folder, prefix.plot = UpdatePrefix(prefix.label, normalization.method), show.plot.flag = FALSE, save.plot = TRUE, plotly.flag = TRUE)
  #   PlotPCAPlotlyFunction(counts.data.frame = normalized.counts, design.matrix = design.matrix, scale = FALSE, plot.folder = plot.folder, prefix.plot = UpdatePrefix(prefix.label, normalization.method), show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, ellipse.flag = "both")
  #   sub.norm.counts <- round(normalized.counts)
  # }
  # 
  # if(de.test=="DeSeq") {
  #   if(!normalize.data.flag) {
  #     # PlotTimesBoxplot(data.frame.to.plot = sub.norm.counts, design.matrix = design.matrix, output.path = plot.folder, prefix.plot = prefix.label, show.plot.flag = FALSE, save.plot = TRUE, plotly.flag = TRUE)
  #     # PlotPCAPlotlyFunction(counts.data.frame = sub.norm.counts, design.matrix = design.matrix, scale = FALSE, plot.folder = plot.folder, prefix.plot = prefix.label, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, ellipse.flag = "both")
  #     # PlotScatterPlotMatrix(counts.data.frame = sub.norm.counts, design.matrix = design.matrix, plot.folder = plot.folder, prefix.plot = prefix.label, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE)
  #     PlotScatterPlotMatrix(counts.data.frame = sub.norm.counts, design.matrix = design.matrix, plot.folder = plot.folder, prefix.plot = prefix.label, show.plot.flag = FALSE, plotly.flag = FALSE, save.plot = TRUE)
  #   }
  # }
  switch(de.test,
         # 
         # DeSeqTime_TC = {
         #   de.data.frame <- ApplyDeSeq(counts.dataframe = sub.norm.counts, design.dataframe = design.matrix, test.type = "LRT_TC")
         # },
         # DeSeqTime_T = {
         #   de.data.frame <- ApplyDeSeq(counts.dataframe = sub.norm.counts, design.dataframe = design.matrix, test.type = "LRT_T")
         # },
         # DeSeqTime_NoInteraction = {
         #   de.data.frame <- ApplyDeSeq(counts.dataframe = sub.norm.counts, design.dataframe = design.matrix, test.type = "LRT_NoInteraction")
         # },
         DeSeq = {
           de.data.frame <- ApplyDeSeq(counts.dataframe = sub.norm.counts, design.dataframe = design.matrix, test.type = "Wald")
         },
         NOISeqBio = {
           de.data.frame <- ApplyNOISeq(counts.dataframe = sub.norm.counts, design.dataframe = design.matrix, test.type = "NOISeqBio")
         },
         NOISeq = {
           de.data.frame <- ApplyNOISeq(counts.dataframe = sub.norm.counts, design.dataframe = design.matrix, test.type = "NOISeq")
         }
         
  )
  
  de.folder <- UpdateFolderPath(de.folder, de.test)
  prefix.label <- paste(prefix.label, de.test)
  filename <- UpdateFilename(filename, de.test)
  
  if(!is.null(conversion.map)){
    mapped.de.data.frame <- AttachConvertedColumn(de.data.frame, conversion.map)
  } else {
    mapped.de.data.frame <- de.data.frame
  }
  
  #de.results = mapped.de.data.frame; counts.dataframe = sub.norm.counts; design.matrix = design.matrix; show.plot.flag = FALSE; plotly.flag = FALSE; save.plot = TRUE; plot.folder = plot.folder; prefix.plot = prefix.label; threshold = threshold
  PlotMAPlotCounts(de.results = mapped.de.data.frame, counts.dataframe = sub.norm.counts, design.matrix = design.matrix, show.plot.flag = FALSE, plotly.flag = FALSE, save.plot = TRUE, plot.folder = plot.folder, prefix.plot = prefix.label, threshold = threshold)
  PlotMAPlotCounts(de.results = mapped.de.data.frame, counts.dataframe = sub.norm.counts, design.matrix = design.matrix, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = plot.folder, prefix.plot = prefix.label, threshold = threshold)
  PlotVolcanoPlot(de.results = mapped.de.data.frame, counts.dataframe = sub.norm.counts, design.matrix = design.matrix, show.plot.flag = FALSE, plotly.flag = FALSE, save.plot = TRUE, plot.folder = plot.folder, prefix.plot = prefix.label, threshold = threshold)
  PlotVolcanoPlot(de.results = mapped.de.data.frame, counts.dataframe = sub.norm.counts, design.matrix = design.matrix, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = plot.folder, prefix.plot = prefix.label, threshold = threshold)
  # WriteDataFrameAsTsv(de.data.frame = de.data.frame, file.name.path = file.path(de.folder,filename))
  
  # de.genes.list <- SignificantDeGenesPAdj(de.data.frame, threshold)
  de.genes.list <- SignificantDeGenesPAdj(mapped.de.data.frame, threshold)
  
  filename.all <- UpdateFilename(filename, "all")
  WriteDataFrameAsTsv(data.frame.to.save = de.genes.list$de.total, file.name.path = file.path(de.folder, filename.all))
  
  filename.not.na <- UpdateFilename(filename, "all_not_na")
  WriteDataFrameAsTsv(data.frame.to.save = de.genes.list$de.not.na, file.name.path = file.path(de.folder, filename.not.na))
  
  filename.sign <- UpdateFilename(filename, paste("significant", threshold))
  WriteDataFrameAsTsv(data.frame.to.save = de.genes.list$SIGN, file.name.path = file.path(de.folder, filename.sign))
  
  filename.up <- UpdateFilename(filename, "UP")
  WriteDataFrameAsTsv(data.frame.to.save = de.genes.list$UP, file.name.path = file.path(de.folder, filename.up))
  
  filename.down <- UpdateFilename(filename, "DOWN")
  WriteDataFrameAsTsv(data.frame.to.save = de.genes.list$DOWN, file.name.path = file.path(de.folder, filename.down))
  
  message(de.test, " produced differential expressed genes...\nTot: ", dim(de.genes.list$SIGN)[1], "\nUP: ", dim(de.genes.list$UP)[1], "\nDOWN: ", dim(de.genes.list$DOWN)[1])
  de.numbers <- rbind(dim(de.genes.list$SIGN)[1], dim(de.genes.list$UP)[1], dim(de.genes.list$DOWN)[1])
  rownames(de.numbers) <- c("SIGN","UP","DOWN")
  colnames(de.numbers) <- filename
  WriteDataFrameAsTsv(data.frame.to.save = de.numbers, file.name.path = file.path(table.folder, filename))
  if(enrich.results.flag){
    enrichSuited(de.gene.list=de.genes.list, functional.folder=functional.folder, filename=filename)
  }
  
  message("returning de list")
  return(de.genes.list)
}


PerformDEAnalysisLRT <- function(whole.counts, design.matrix, de.test=c("DeSeqTime_TC", "DeSeqTime_T", "DeSeqTime_NoInteraction"), results.folder=NULL, prefix.label = NULL, normalize.data.flag=FALSE, normalization.method="uqua", negative.genes.list=NULL, threshold=0.05, conversion.map=NULL, enrich.results.flag=TRUE) {
  #check input data
  sub.norm.counts <- whole.counts[, which(colnames(whole.counts) %in% rownames(design.matrix))]
  
  
  plot.folder <- UpdateFolderPath(results.folder, "plots")
  de.folder.root <- UpdateFolderPath(results.folder, "DE_results")
  table.folder.root <- UpdateFolderPath(results.folder, "tables")
  functional.folder.root <- UpdateFolderPath(results.folder, "functional")
  
  filename.root <- prefix.label
  prefix.label.root <- prefix.label
  de.data.list <- NULL
  
  switch(de.test,
         
         DeSeqTime_TC = {
           de.data.list <- ApplyDeSeq(counts.dataframe = sub.norm.counts, design.dataframe = design.matrix, test.type = "LRT_TC")
         },
         DeSeqTime_T = {
           de.data.list <- ApplyDeSeq(counts.dataframe = sub.norm.counts, design.dataframe = design.matrix, test.type = "LRT_T")
         },
         DeSeqTime_NoInteraction = {
           de.data.list <- ApplyDeSeq(counts.dataframe = sub.norm.counts, design.dataframe = design.matrix, test.type = "LRT_NoInteraction")
         }
  )

  if(!is.null(de.data.list)) {
    de.genes.data.list <- list()
    i<-1
    for(de.data.frame in de.data.list) {
      de.test <- names(de.data.list)[i]
      table.folder <- UpdateFolderPath(table.folder.root, de.test)
      functional.folder <- UpdateFolderPath(functional.folder.root, de.test)
      de.folder <- UpdateFolderPath(de.folder.root, de.test)
      prefix.label <- paste(prefix.label.root, de.test)
      filename <- UpdateFilename(filename.root, de.test)
      
      if(!is.null(conversion.map)) {
        mapped.de.data.frame <- AttachConvertedColumn(de.data.frame, conversion.map)
      }
      
      #de.results = mapped.de.data.frame; counts.dataframe = sub.norm.counts; design.matrix = design.matrix; show.plot.flag = FALSE; plotly.flag = FALSE; save.plot = TRUE; plot.folder = plot.folder; prefix.plot = prefix.label; threshold = threshold
      PlotMAPlotCounts(de.results = mapped.de.data.frame, counts.dataframe = sub.norm.counts, design.matrix = design.matrix, show.plot.flag = FALSE, plotly.flag = FALSE, save.plot = TRUE, plot.folder = plot.folder, prefix.plot = prefix.label, threshold = threshold)
      PlotMAPlotCounts(de.results = mapped.de.data.frame, counts.dataframe = sub.norm.counts, design.matrix = design.matrix, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = plot.folder, prefix.plot = prefix.label, threshold = threshold)
      PlotVolcanoPlot(de.results = mapped.de.data.frame, counts.dataframe = sub.norm.counts, design.matrix = design.matrix, show.plot.flag = FALSE, plotly.flag = FALSE, save.plot = TRUE, plot.folder = plot.folder, prefix.plot = prefix.label, threshold = threshold)
      PlotVolcanoPlot(de.results = mapped.de.data.frame, counts.dataframe = sub.norm.counts, design.matrix = design.matrix, show.plot.flag = FALSE, plotly.flag = TRUE, save.plot = TRUE, plot.folder = plot.folder, prefix.plot = prefix.label, threshold = threshold)
      # WriteDataFrameAsTsv(de.data.frame = de.data.frame, file.name.path = file.path(de.folder,filename))
      
      if(!is.null(conversion.map)) {
        de.genes.list <- SignificantDeGenesPAdj(mapped.de.data.frame, threshold)
      } else {
        de.genes.list <- SignificantDeGenesPAdj(de.data.frame, threshold)
      }
      
      
      filename.all <- UpdateFilename(filename, "all")
      WriteDataFrameAsTsv(data.frame.to.save = de.genes.list$de.total, file.name.path = file.path(de.folder, filename.all))
      
      filename.not.na <- UpdateFilename(filename, "all_not_na")
      WriteDataFrameAsTsv(data.frame.to.save = de.genes.list$de.not.na, file.name.path = file.path(de.folder, filename.not.na))
      
      filename.sign <- UpdateFilename(filename, paste("significant", threshold))
      WriteDataFrameAsTsv(data.frame.to.save = de.genes.list$SIGN, file.name.path = file.path(de.folder, filename.sign))
      
      filename.up <- UpdateFilename(filename, "UP")
      WriteDataFrameAsTsv(data.frame.to.save = de.genes.list$UP, file.name.path = file.path(de.folder, filename.up))
      
      filename.down <- UpdateFilename(filename, "DOWN")
      WriteDataFrameAsTsv(data.frame.to.save = de.genes.list$DOWN, file.name.path = file.path(de.folder, filename.down))
      
      message(de.test, " produced differential expressed genes...\nTot: ", dim(de.genes.list$SIGN)[1], "\nUP: ", dim(de.genes.list$UP)[1], "\nDOWN: ", dim(de.genes.list$DOWN)[1])
      de.numbers <- rbind(dim(de.genes.list$SIGN)[1], dim(de.genes.list$UP)[1], dim(de.genes.list$DOWN)[1])
      rownames(de.numbers) <- c("SIGN","UP","DOWN")
      colnames(de.numbers) <- filename
      WriteDataFrameAsTsv(data.frame.to.save = de.numbers, file.name.path = file.path(table.folder, filename))
      
      if(enrich.results.flag) {
        enrichSuited(de.gene.list=de.genes.list, functional.folder=functional.folder, filename=filename)
      }
      
      
      de.genes.data.list[[de.test]] <- de.genes.list
      i<-i+1
    }
    
    
    if(enrich.results.flag) {
      de.sign.genes.union <- union(rownames(de.genes.data.list[[names(de.data.list)[1]]][["SIGN"]]), rownames(de.genes.data.list[[names(de.data.list)[2]]][["SIGN"]]))
      functional.folder <- UpdateFolderPath(functional.folder.root, UpdateFilename(names(de.data.list)[1], names(de.data.list)[2], "union"))
      filename <- UpdateFilename(filename.root, names(de.data.list)[1], names(de.data.list)[2], "union")
      enrichSuitedSingleList(de.gene.list = de.sign.genes.union, functional.folder = functional.folder, filename = filename)
    }
    
  }
  message("returning de lists")
  return(de.genes.data.list)
}


enrichSuitedSingleList <- function(de.gene.list, functional.folder, filename, organism="rnorvegicus") {
  
  switch(organism,
         rnorvegicus = {
           org.db = "org.Rn.eg.db"
           org.code = "rno"
         },
         mmusculus = {
           org.db = "org.Mm.eg.db"
           org.code = "mmu"
         }
  )
  
    dir.create(functional.folder, showWarnings = FALSE, recursive = TRUE)
    tryCatch({
      enrichKEGGFunction(gene.list.to.enrich = de.gene.list, functional.folder = functional.folder, filename = filename, organism.db = org.db, organism.code = org.code)
      for(path in c("KEGG", "REAC")) {
        enrichPathwayGProfiler(gene.list.to.enrich = de.gene.list, functional.folder = functional.folder, filename = filename, path.label = path, organism.name = organism)
      }
    }, error=function(e) {
      print(e)
    }
    )
    
    for (go.label in c("BP", "CC", "MF")) {
      tryCatch({
        enrichGOFunction(gene.list.to.enrich = de.gene.list, functional.folder = functional.folder, filename = filename, ontology = go.label, organism.db = org.db, organism.code = org.code)
        enrichGOGProfiler(gene.list.to.enrich = de.gene.list, functional.folder = functional.folder, filename = filename, ontology = go.label, organism.name = organism)
      }, error=function(e) {
        print(e)
      }
      )
    }
  
}

enrichSuited <- function(de.gene.list, functional.folder, filename, organism="rnorvegicus") {
  
  switch(organism,
         rnorvegicus = {
           org.db = "org.Rn.eg.db"
           org.code = "rno"
         },
         mmusculus = {
           org.db = "org.Mm.eg.db"
           org.code = "mmu"
         }
         )
  
  for (df.label in c("SIGN", "UP", "DOWN")) {
    switch(df.label,
           SIGN = { 
             filename.suit <- UpdateFilename(filename, "SIGN")
           },
           UP = { 
             filename.suit <- UpdateFilename(filename, "UP")
           },
           DOWN = { 
             filename.suit <- UpdateFilename(filename, "DOWN")
           }
    )
    
    tryCatch({
      enrichKEGGFunction(gene.list.to.enrich = rownames(de.gene.list[[df.label]]), functional.folder = functional.folder, filename = filename.suit, organism.db = org.db, organism.code = org.code)
      for(path in c("KEGG", "REAC")) {
        enrichPathwayGProfiler(gene.list.to.enrich = rownames(de.gene.list[[df.label]]), functional.folder = functional.folder, filename = filename.suit, path.label = path, organism.name = organism)
      }
      
    }, error=function(e) {
      print(e)
      
    }
    )
    
    for (go.label in c("BP", "CC", "MF")) {
      tryCatch({
        enrichGOFunction(gene.list.to.enrich = rownames(de.gene.list[[df.label]]), functional.folder = functional.folder, filename = filename.suit, ontology = go.label, organism.db = org.db, organism.code = org.code)
        enrichGOGProfiler(gene.list.to.enrich = rownames(de.gene.list[[df.label]]), functional.folder = functional.folder, filename = filename.suit, ontology = go.label, organism.name = organism)
      }, error=function(e) {
        print(e)
      }
      )
    }
    
  }
  
}

VennDE3ListsSuited <- function(de.df.list1, de.df.list2, de.df.list3, title, venn.plot.dir, intersection.flag = TRUE, intersection.exclusion.flag = TRUE){ 
  
  for (df.label in c("SIGN", "UP", "DOWN")) {
    switch(df.label,
           SIGN = { 
             title.p = paste(title, "Significant DE Genes") 
             label1.p = "cervical"
             label2.p = "thoracic"
             label3.p = "cerv-thor"
             },
           UP = { 
             title.p = paste(title, "UP-regulated DE Genes") 
             label1.p = paste("cervical", "UP", sep = "_")
             label2.p = paste("thoracic", "UP", sep = "_")
             label3.p = paste("cerv-thor", "UP", sep = "_")
           },
           DOWN = { 
             title.p = paste(title, "DOWN-regulated DE Genes") 
             label1.p = paste("cervical", "DOWN", sep = "_")
             label2.p = paste("thoracic", "DOWN", sep = "_")
             label3.p = paste("cerv-thor", "DOWN", sep = "_")
           }
    )
    
    Venn3de(x = rownames(de.df.list1[[df.label]]), y = rownames(de.df.list2[[df.label]]), z = rownames(de.df.list3[[df.label]]), 
            label1 = label1.p, label2 = label2.p, label3 = label3.p, title = title.p, 
            intersection.flag = intersection.flag, intersection.exclusion.flag = intersection.exclusion.flag, plot.dir = venn.plot.dir)
    
  }
}

AttachConvertedColumn <- function(de.data.frame, conversion.map) {
  de.data.frame.ord <- de.data.frame[ order(rownames(de.data.frame)),]
  conversion.map.ord <- conversion.map[ order(conversion.map[,1]), ]
  ind.map <- which(conversion.map.ord[,1] %in% rownames(de.data.frame.ord) )
  conversion.map.sub <- conversion.map.ord[ind.map,]
  if(dim(conversion.map.sub)[1]==dim(de.data.frame.ord)[1]) {
    de.data.frame.ord$symbol <- as.character(conversion.map.sub[,2])
  } else {
    stop("gene map and de.dataframe dimensions differs!")
  }
  return(de.data.frame.ord)
}


CreateConvertedDataframe <- function(gene.list, conversion.map) {
  # de.data.frame.ord <- de.data.frame[ order(rownames(de.data.frame)),]
  gene.list.ord <- as.data.frame(gene.list[order(gene.list)])
  conversion.map.ord <- conversion.map[ order(conversion.map[,1]), ]
  ind.map <- which(conversion.map.ord[,1] %in% gene.list.ord[,1] )
  conversion.map.sub <- conversion.map.ord[ind.map,]
  if(dim(conversion.map.sub)[1]==dim(gene.list.ord)[1]) {
    gene.list.ord$symbol <- as.character(conversion.map.sub[,2])
  } else {
    error("gene map and de.dataframe dimensions differs!")
  }
  colnames(gene.list.ord) <- colnames(conversion.map)
  
  return(gene.list.ord)
}






