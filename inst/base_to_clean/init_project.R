IniziateProject <- function(project.name=NULL, prefix=NULL) {
  if(is.null(project.name) || is.null(prefix)) {
    stop("Please provide a project.name and a prefix!")
  }
  
  if(is.character(project.name)) {
    Project.Name <<- project.name
  } else {
    stop("Please provide a project.name as character string!")
  }
    
  if(is.character(prefix)) {
    Prefix <<- prefix
  } else {
    stop("Please provide a prefix as character string!")
  }
  
  main.directory <- file.path(getwd(), "TiCoRSe", Project.Name)
  
  dir.create(path = main.directory, recursive = TRUE)
  
}
