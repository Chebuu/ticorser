initReportFilename=function(filenamepath=NULL, mainTitle=NULL,
                            authors=NULL, documentType="html")
{
    if(is.null(filenamepath)) stop("Report file name is NULL!")
    if(is.null(mainTitle)) stop("Please provide a title for the document!")
    documentType <- match.arg(documentType)
    if(!file.exists(filenamepath))
    {
        if(tools::file_ext(filenamepath) != "Rmd")
        {
            filenamepath <- paste0(filenamepath, ".Rmd")
        }
        file.create(filenamepath)
    }

    private$filenamePath <- filenamepath

    header <- paste0("---
    title: \"", mainTitle, "\"
    author: \"", authors, "\"
    date: \"`r Sys.Date()`\"
    output: rmarkdown::", documentType, "_document
---")

    base::write(header, file = filenamepath,
                append = TRUE, sep = "\n")
    private$markdownSetGlobalOpts()
}

markdownSetGlobalOpts <- function(evalFlag=TRUE, echoFlag=TRUE,
                                  warningFlag=FALSE, showMessages=FALSE,
                                  includeFlag=FALSE)
{
    options <- paste0("```{r global_options, include=FALSE}\n",
                    "knitr::opts_chunk$set(",
                    "eval=", evalFlag,
                    ", echo=", echoFlag,
                    ", warning=", warningFlag,
                    ", message=", showMessages,
                    ", include=", includeFlag,
                    ")\n```")
    base::write(options, file = private$filenamePath,
                append = TRUE, sep = "\n")

}
