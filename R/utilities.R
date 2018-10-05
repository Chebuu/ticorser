
#' WriteDataFrameAsTsv
#' just save a dataframe on disk in TSV format
#' @param data.frame.to.save the data frame to store on disk
#' @param file.name.path the file path where to store the dataset
#' @param col.names see ?write.table description (default is NA)
#' @param row.names see ?write.table description (default is TRUE)
#'
#' @return none
#' @export
#'
#' @examples
#' df <- data.frame(1, 1:10, rownames=paste0(rep("rn", 10), 1:10),
#'     colnames=c("c1","c2"))
#' WriteDataFrameAsTsv(df)
WriteDataFrameAsTsv <- function(data.frame.to.save=NULL,
                            file.name.path=file.path(tempdir(),tempfile()),
                            col.names=NA, row.names=TRUE)
{
    if(is.null(data.frame.to.save)) stop("Please provide a data frame!")
    file.name.path <- gsub(pattern=" ", replacement="_", x=file.name.path)
    file.name <- paste0(file.name.path, ".tsv")
    write.table(x=data.frame.to.save, file=file.name, quote=FALSE, sep="\t",
                col.names=col.names, row.names=row.names)
    message(file.name, " written on disk as TSV file!\n")
}

#' ReadDataFrameFromTsv
#'
#' @param file.name.path a TSV format filename with path
#' @param row.names.col the column number where the rownames are located (def.1)
#' @param header.flag is an HEADER present? (default is TRUE)
#' @param sep the column separator (default is "\t")
#' @param quote.char the character used for quoting characters (defualt is none)
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' TBW
ReadDataFrameFromTsv <- function(file.name.path, row.names.col=1,
                                header.flag=TRUE, sep="\t", quote.char="")
{
    if(!file.exists(file.name.path)) stop("Please provide an existing filename")
    df <- read.table(file=file.name.path, sep=sep, header=header.flag,
                     row.names=row.names.col, quote=quote.char)
    message(file.name.path, " read from disk!\n")
    return(df)
}


#' convertGenesViaMouseDb
#'
#' @description converts genes from SYMBOL or ENTREZID to SYMBOL or ENTREZID for
#' mouse genome as defined in org.Mm.eg.db
#'
#' @param gene.list a list of genes
#' @param fromType one of SYMBOL or ENTREZID
#' @param toType one of SYMBOL or ENTREZID
#'
#' @return the gene.map with only converted values.
#'
#' @export
#' @import org.Mm.eg.db
#' @importFrom AnnotationDbi select
#' @examples
convertGenesViaMouseDb <- function(gene.list, fromType=c("SYMBOL", "ENTREZID"),
                                   toType=c("SYMBOL", "ENTREZID"))
{
    require("org.Mm.eg.db")
    annotated.map <- AnnotationDbi::select(org.Mm.eg.db,
                                           keys=gene.list,
                                           keytype=fromType,
                                           columns=c("SYMBOL", "ENTREZID"))
    col.check <- colnames(annotated.map)[-which(colnames(annotated.map) ==
                                                    fromType)]
    indsna <- is.na(annotated.map[col.check])
    if(sum(indsna) > 0) annotated.map <- annotated.map[-which(indsna),]
    return(annotated.map)

}


#' convertGenesViaBiomart
#' @description Converts a list of genes using biomaRt package. (See biomaRt for
#' further usage description)
#' @param specie one of "hg38", "mm10", "rnor6".
#' @param attrs one or more attributes to be returned
#' as defined in biomaRt::listAttrs
#' @param filter one or more filters as defined in biomaRt::listFilters
#' @param filter.values a list of values for the defined filters
#'
#' @return the gene.map of attributes
#'
#' @importFrom biomaRt useMart getBM
#' @export
#'
#' @examples
#' TBW
convertGenesViaBiomart <- function(specie=c("hg38", "mm10", "rnor6"),
                                   attrs=NULL, filter=NULL, filter.values=NULL)
{
    specie <- match.arg(specie)
    stopifnot((length(attrs) != 0))
    stopifnot( is.null(filter) || (!is.null(filter.values)) )
    stopifnot( (!is.null(filter)) || is.null(filter.values) )

    if(length(which(attrs %in% filter))==0)
    {
        stop("please use a filter matching the attributes!")
    }

    switch (specie,
            "mm10"={ ds <- "mmusculus_gene_ensembl"},
            "hg38"={ ds <- "hsapiens_gene_ensembl"},
            "rnor6"={ ds <- "rnorvegicus_gene_ensembl"}
    )

    mart <- biomaRt::useMart("ensembl", dataset=ds)

    # listAttributes(mart)[grep(external", listAttributes(mart)[,1]),1]
    # listFilters(mart)[grep("external", listFilters(mart)[,1]),1]
    # attrs <- c("ensembl_gene_id", "external_gene_name", "entrezgene")
    gene.map <- biomaRt::getBM(attributes=attrs, mart=mart,
                               filters=filter, values=filter.values)
    if(!is.null(filter))
    {
        idx.dp <- which(duplicated(gene.map[[filter]]))
        if(length(idx.dp) > 0 )
        {
            gene.map <- gene.map[-idx.dp,]
        }
        gene.map <- gene.map[order(gene.map[[filter]]),]
    }
    idx.entr <- grep(pattern="entrezgene", x=colnames(gene.map))
    if(length(idx.entr) > 0 )
    {
        gene.map[,idx.entr] <- as.character(gene.map[,idx.entr])
    }

    return(gene.map)
}


#' attachGeneColumnToDf
#' attaches an additional column to a dataframe using its rownames as
#' identifiers
#' @param mainDf the dataframe where the
#' @param genesMap
#' @param rowNamesIdentifier
#' @param mapFromIdentifier
#' @param mapToIdentifier
#'
#' @return
#' @export
#'
#' @examples
attachGeneColumnToDf <- function(mainDf, genesMap,
                     rowNamesIdentifier=c("entrezgene", "ensembl", "symbol"),
                     mapFromIdentifier=NULL, mapToIdentifier=NULL)
{
    match.arg(rowNamesIdentifier)
    stopifnot(!is.null(mapFromIdentifier))
    stopifnot(!is.null(mapToIdentifier))

    mainDf <- mainDf[order(rownames(mainDf)),]
    rownames <- rownames(mainDf)

    mainDf$check <- NA
    mainDf$gene <- NA

    genesMap <- genesMap[order(genesMap[[mapFromIdentifier]]),]

    idx.re <- which(rownames %in% genesMap[[mapFromIdentifier]])
    idx.er <- which(genesMap[[mapFromIdentifier]] %in% rownames)

    mainDf$check[idx.re] <- genesMap[[mapFromIdentifier]][idx.er]
    mainDf$gene[idx.re] <- genesMap[[mapToIdentifier]][idx.er]
    ## removing NA (not mapped) lines
    # noNaMainDf <- mainDf
    # idx.na <- which(is.na(mainDf$gene))
    # if(length(idx.na) > 0) noNaMainDf <- mainDf[-idx.na,]
    # idx.e <- which(noNaMainDf$gene == "")
    # if(length(idx.e) > 0) noNaMainDf <- noNaMainDf[-idx.e,]
    # return(noNaMainDf)
    return(mainDf)

}


#' subsetDfByCol
#'
#' @description Subsets a dataframe by a list of colnames or rownames
#'
#' @param df a dataframe.
#' @param list a list of identifiers.
#' @param colname the df column where to check the list in the df.
#' If NULL, the rownames will be used.
#'
#' @return the subsetted df
#' @export
#'
#' @examples
subsetDfByCol <- function(df, list, colname=NULL)
{
    if(is.null(colname))
    {
        idx <- which(rownames(df) %in% list)
    } else {
        idx <- which(df[[colname]] %in% list)
    }
    df <- df[idx, , drop=FALSE]
}


#' isCol
#' @description checks if a column is present in a dataframe/matrix.
#' @param df.to.check the dataframe/matrix
#' @param colname the colname to check
#'
#' @return boolean value
#' @export
#'
#' @examples
#'
isCol <- function(df.to.check, colname)
{
    idx <- which(colnames(df.to.check) %in% colname)
    if(length(idx) > 0)
    {
        return(TRUE)
    }else{
        return(FALSE)
    }
}


#' generatePlotStrings
#' @description generates a list of title, plot.folder and plot.file.name
#' strings. Typically used in plot functions.
#'
#' @param path the starting path
#' @param prefix the prefix for the file
#' @param plot.type the type of the plot
#'
#' @return a list of strings title, plot.folder and plot.file.name.
#' @keywords internal
#'
#' @examples
generatePlotStrings <- function(path=NULL, prefix, plot.type)
{
    title <- gsub(pattern="_", replacement=" ",
                x=UpdatePrefix(prefix, plot.type))

    plot.folder <- gsub(pattern=" ", replacement="_",
                        x=file.path(path, plot.type))

    plot.file.name <- gsub(pattern=" ", replacement="_",
                        x=UpdatePrefix(prefix, plot.type))

    if(!is.null(path)) dir.create(plot.folder, showWarnings=FALSE,
                                recursive=TRUE)

    return(list("title"= title,
                "plot.folder"=plot.folder,
                "plot.file.name"=plot.file.name))
}


#' updatePrefix
#' @description given an input string, it appends one or more strings to it,
#' separated by spaces. (tipically used by generatePlotStrings function)
#' @param prefix the string of the prefix to update
#' @param ... a list of strings to append to the prefix
#'
#' @return the updated prefix
#' @export
#'
#' @examples
#'
updatePrefix <- function(prefix, ...)
{
    dots <- unlist(list(...))
    if(length(dots) != 0)
    {
        for (str in dots)
        {
            # str <- gsub(pattern=".", replacement="_", str)
            prefix <- paste(prefix, str, sep=" ")
        }

    } else {
        stop("provide a string to append to ", new.prefix)
    }
    return(prefix)
}

#' updateFolderPath
#' @description appends one or more strings to the path and creates all
#' the directories in the paths recursively.
#' Additionally, replaces all the whitespaces with underscores.
#' @param path a string representing a path
#' @param ... a list of one or more elements to append to the path
#'
#' @return the string of the updated path
#' @export
#'
#' @examples
#' pth <- "./old/path"
#' updateFolderPath(path=pth, c("new", "directories", "in_the", "path"))
#'
updateFolderPath <- function(path, ...)
{
    dots <- unlist(list(...))
    if(length(dots) != 0)
    {
        for (str in dots)
        {
            str <- gsub(pattern=" ", replacement="_", str)
            path <- file.path(path, str)
        }
    } else {
        stop("provide a string to append to ", path)
    }
    dir.create(path, recursive=TRUE, showWarnings=FALSE)
    message("Recursively created ", path, " on disk")
    return(path)
}

#' updateFilename
#' @description appends one or more strings to a filename string
#' @param filename the string representing the starting filename
#' @param ... a list of strings to append to filename
#' @param extension an extension to add to the filename without . (optional)
#'
#' @return a string with the updated filename
#' @export
#'
#' @examples
#' fn <- file1
#' newfn <- updateFilename(filename=fn, c("with", "more", "informations"),
#'                         extension="pdf")
#' print(newfn)
#'
updateFilename <- function(filename, ..., extension=NULL)
{
    dots <- unlist(list(...))
    print(dots)
    filename <- gsub(pattern=" ", replacement="_", x=filename)
    if(length(dots) != 0)
    {
        for (str in dots) filename <- paste(filename, str, sep="_")
    } else {
        stop("provide a string to append to ", filename)
    }
    if(!is.null(extension))
    {
        filename <- paste0(filename, ".", extension)
    }
    return(filename)
}
