GeneratePlotStrings <- function(path=NULL, prefix, plot.type) {
  title <- gsub(pattern = "_", replacement = " ", x = UpdatePrefix(prefix, plot.type))

  plot.folder <- gsub(pattern = " ", replacement = "_", x = file.path(path, plot.type))

  plot.file.name <- gsub(pattern = " ", replacement = "_", x = UpdatePrefix(prefix, plot.type))
  if(!is.null(path)) dir.create(plot.folder, showWarnings = FALSE, recursive = TRUE)

  return(list("title"= title, "plot.folder"=plot.folder, "plot.file.name"=plot.file.name))
}


UpdatePrefix <- function(prefix, ...) {
  # new.prefix <- paste(prefix, postix, sep=sep)
  dots <- list(...)
  if(length(dots) != 0) {
    for (str in dots) {
      # str <- gsub(pattern = ".", replacement = "_", str)
      prefix <- paste(prefix, str, sep = " " )
    }
    
  }else {
    stop("provide a string to append to ", new.prefix)
  }
  return(prefix)
}

UpdateFolderPath <- function(path, ...) {
  dots <- list(...)
  if(length(dots) != 0) {
    for (str in dots) {
        str <- gsub(pattern = " ", replacement = "_", str)
        path <- file.path(path, str)
      }
  } else {
    stop("provide a string to append to ", path)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  return(path)
}

UpdateFilename <- function(filename, ..., extension=NULL) {
  dots <- list(...)
  filename <- gsub(pattern = " ", replacement = "_", x = filename)
  if(length(dots) != 0) {
    for (str in dots) filename <- paste(filename, str, sep = "_")
  } else {
    stop("provide a string to append to ", filename)
  }
  if(!is.null(extension)) {
    filename <- paste0(filename, ".", extension)
  }
  return(filename)
}
# 
# ProduceOutputStrings <- function(prefix, out.dir, ..., filename=NULL) {
#   if(class(...)!="character") {
#     dots <- list(...)
#   } else {
#     dots <- 
#   }
#   
#   
#   for(lbl in labels.list) {
#     prefix <- UpdatePrefix(prefix, lbl)
#   }
#   
#   out.dir <- UpdateFolderPath(root.dir, prefix)
#   filename <- UpdateFilename(prefix, "genes")
#   
#   
#   
# }

WriteDataFrameAsTsv <- function(data.frame.to.save, file.name.path, col.names=NA, row.names=TRUE) {
  file.name.path <- gsub(pattern = " ", replacement = "_", x = file.name.path)
  file.name <- paste0(file.name.path, ".tsv")
  write.table(x = data.frame.to.save, file = file.name, quote = FALSE, sep = "\t", col.names = col.names, row.names = row.names)
  message(file.name, " written on disk as TSV file!\n")
}

ReadDataFrameFromTsv <- function(file.name.path, row.names.col=1, header.flag=TRUE, sep="\t", quote.char="") {
  df <- read.table(file = file.name.path, sep = sep, header=header.flag, row.names = row.names.col, quote = quote.char)
  message(file.name.path, " read from disk!\n")
  return(df)
}

# ## global variables
# project.name <<- "cervical_thoracic_whole"
# results.folder <<- UpdateFolderPath(project.name, "results")
# counts.folder <<- UpdateFolderPath(results.folder, "counts")
# plots.folder <<- UpdateFolderPath(results.folder, "plots")
# de.plots.folder <<- UpdateFolderPath(plots.folder, "de_plots")
# boxplot.folder <<- UpdateFolderPath(plots.folder, "boxplots")
# pca.folder <<- UpdateFolderPath(plots.folder, "PCA")
# file.name <<- project.name



