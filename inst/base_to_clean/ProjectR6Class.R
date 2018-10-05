require("R6")
Project <- R6Class( "Project",
                    
                        private = list (
                          
                        
                        ) ,
                    
                        public = list(
                          
                          name = "character",
                          
                          working.path = "character",
                          results.path = "character",
                          counts.path = "character",
                          de.path = "character",
                          plots.path = "character",
                          functional.path = "character",
                          already.exists = "logical",
                          filename = "character",
                          normalization.path = "character",
                          filtering.path = "character",
                          
                          
                          initialize = function(name, main.working.path) {
                            self$name <- self$SetName(name)
                            self$results.path <- UpdateFolderPath(main.working.path, name, "results")
                            self$counts.path <- UpdateFolderPath(self$results.path, "counts")
                            self$plots.path <- UpdateFolderPath(self$results.path, "plots")
                            self$de.path <- UpdateFolderPath(self$results.path, "de_results")
                            self$functional.path <- UpdateFolderPath(self$results.path, "functional_analysis")
                            self$working.path <- self$results.path
                          }
                          ,
                          
                          finalize = function()
                          {
                            
                          }
                          ,
                          
                          GeneratePlotStrings = function(path, prefix, plot.type) {
                            title <- gsub(pattern = "_", replacement = " ", x = UpdatePrefix(prefix, plot.type))
                            
                            plot.folder <- gsub(pattern = " ", replacement = "_", x = file.path(path, plot.type))
                            
                            plot.file.name <- gsub(pattern = " ", replacement = "_", x = UpdatePrefix(prefix, plot.type))
                            
                            dir.create(plot.folder, showWarnings = FALSE, recursive = TRUE)
                            
                            return(list("title"= title, "plot.folder"=plot.folder, "plot.file.name"=plot.file.name))
                          }
                          ,
                          
                          UpdatePrefix = function(prefix, ...) {
                            # new.prefix <- paste(prefix, postix, sep=sep)
                            dots <- list(...)
                            if(length(dots) != 0) {
                              for (str in dots) new.prefix <- paste(prefix, str, sep = " " )
                            }else {
                              stop("provide a string to append to ", new.prefix)
                            }
                            return(new.prefix)
                          }
                          ,
                          
                          UpdateFolderPath = function(path, ...) {
                            dots <- list(...)
                            if(length(dots) != 0) {
                              for (str in dots) path <- file.path(path, str)
                            }else {
                              stop("provide a string to append to ", path)
                            }
                            dir.create(path, recursive = TRUE, showWarnings = FALSE)
                            return(path)  
                          }
                          ,
                          
                          GetProjectPaths = function() {
                            return("name" =self$name
                                    , "working.path" = self$working.path 
                                    , "results.path" = self$results.path 
                                    , "counts.path" = self$counts.path 
                                    , "plots.path" = self$plots.path 
                                    , "de.path" = self$de.path
                                    , "functional.path" = self$functional.path
                                   )
                          }
                          ,
                          
                          UpdateWorkingPath = function(sub.folder.to.generate) {
                            self$working.path <- self$UpdateFolderPath(self$working.path, sub.folder.to.generate)
                            return(self$working.path)
                          }
                          ,
                          
                          GetWorkOnPlots = function() {
                            self$working.path <- self$plots.path
                            return(self$working.path)
                          }
                          ,
                          
                          GetWorkOnCounts = function() {
                            self$working.path <- self$counts.path
                            return(self$working.path)
                          }
                          ,
                          
                          GetWorkOnDE = function() {
                            self$working.path <- self$de.path
                            return(self$working.path)
                          }
                          ,
                          
                          GetWorkOnFunctional = function() {
                            self$working.path <- self$functional.path
                            return(self$working.path)
                          }
                          ,
                          
                          SetName = function(string) {
                            self$name <- string
                          }
                          ,
                          
                          SetNormalizationPath = function(normalization.name) {
                            self$normalization.path <- self$UpdateFolderPath(self$filtering.path, normalization.name)
                          }
                          ,
                          
                          SetFilteringPath = function(filtering.name) {
                            self$filtering.path <- self$UpdateFolderPath(self$counts.path, filtering.name)
                          }
                          
                        )
                        
                          
)


