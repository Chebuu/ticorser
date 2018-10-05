
#' filterLowCounts
#'
#' @param counts.dataframe
#' @param is.normalized
#' @param design.dataframe
#' @param cond.col.name
#' @param method.type
#' @param cv.percentage
#' @param cpm.cutoff
#' @param seq.depth
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
filterLowCounts <- function(counts.dataframe, is.normalized = c(TRUE, FALSE),
                            design.dataframe, cond.col.name=NULL,
                            method.type = c("CPM", "Wilcoxon", "Proportion", "quantile"),
                            cv.percentage=100, cpm.cutoff=1, seq.depth=NULL,
                            quartile.threshold=0.99,
                            verbose=TRUE)
{

    method.type <- match.arg(method.type)
    stopifnot(is.logical(is.normalized))
    stopifnot(!is.null(cond.col.name))
    stopifnot((cond.col.name %in% colnames(design.dataframe)))


    design.dataframe <- design.dataframe[order(rownames(design.dataframe)), ,
                                         drop=FALSE]
    counts.dataframe <- counts.dataframe[, order(colnames(counts.dataframe)),
                                         drop=FALSE]

    idx.cols <- which(colnames(counts.dataframe) %in% rownames(design.dataframe))
    sub.counts.dataframe <- counts.dataframe[, idx.cols]

    conditions <- design.dataframe[, cond.col.name]

    if(method.type == "Proportion")
    {
        if(is.normalized) {
            ## calcolare seq.depth
            if(is.null(seq.depth))
                stop("Proportion test cannot be performed on normalized counts",
                     " without sequencing depth!\nYou need column totals ",
                     " before normalizing the data.")
        }
        else
        {
            seq.depth <- NULL
        }
    }

    if(verbose) message("features amount before normalization: ",
                        dim(sub.counts.dataframe)[1])

    switch( method.type,
            CPM = {
                method.number <- 1
            },
            Wilcoxon = {
                method.number <- 2
            },
            Proportion = {
                method.number <- 3
            },
            quantile = {
                filtered.data.frame <- filterOutGenesOverQuartile(
                    samples.data.frame=sub.counts.dataframe,
                    quartile.threshold=quartile.threshold
                )
                if(verbose) message("features amount before normalization: ",
                                    dim(filtered.data.frame)[1])
                return(filtered.data.frame)
            }
    )



    filtered.dataframe <- NOISeq::filtered.data(dataset = sub.counts.dataframe,
                                                factor = conditions,
                                                norm = is.normalized,
                                                depth = seq.depth,
                                                method = method.number,
                                                cv.cutoff = cv.percentage,
                                                cpm = cpm.cutoff,
                                                p.adj = "BH")
    if(verbose) message("features amount before normalization: ",
                        dim(filtered.data.frame)[1])
    return(filtered.dataframe)
}


filterOutGenesOverQuartile <- function(samples.data.frame,
                                       quartile.threshold=0.99)
{########### TO REFINE 5-OCT-18
    temporary.data.frame <- samples.data.frame
    quartiles <- apply(X=samples.data.frame, MARGIN=2,
                       quantile,
                       quartile.threshold)

    samples <- colnames(samples.data.frame)
    for(sample in samples)
    {
        temporary.data.frame[
            which(samples.data.frame[,sample] >
                      quartiles[sample]),sample] <- NA
    }
    arr.ind <- which(is.na(temporary.data.frame), arr.ind=TRUE)
    arr.ind.ro <- arr.ind[order(arr.ind[,"row"]),]
    row.unq <- unique(arr.ind.ro[,"row"])

    numcols <- dim(temporary.data.frame)[2]
    another.df <- temporary.data.frame
    for(row.ind in row.unq)
    {
        row.ind.ind <- which(arr.ind.ro[,"row"] %in% row.ind)
        if(length(row.ind.ind) == numcols)
        {
            message("deleting ", row.ind, " row")
            another.df <- another.df[-row.ind, ]
        }
    }
    return(another.df)
}
