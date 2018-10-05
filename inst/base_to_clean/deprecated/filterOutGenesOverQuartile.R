filterOutGenesOverQuartile <- function(samples.data.frame, quartile.threshold=0.99) {
  temporary.data.frame <- samples.data.frame
  quartiles <- apply(X = samples.data.frame, MARGIN = 2, quantile, quartile.threshold)
  samples <- colnames(samples.data.frame)
  for(sample in samples) {
    temporary.data.frame[which(samples.data.frame[,sample] > quartiles[sample]),sample] <- NA
  }
  arr.ind <- which(is.na(temporary.data.frame), arr.ind = TRUE)
  arr.ind.ro <- arr.ind[order(arr.ind[,"row"]),]
  row.unq <- unique(arr.ind.ro[,"row"])
  
  numcols <- dim(temporary.data.frame)[2]
  another.df <- temporary.data.frame
  for(row.ind in row.unq) {
    row.ind.ind <- which(arr.ind.ro[,"row"] %in% row.ind)
    if(length(row.ind.ind) == numcols) {
      message("deleting ", row.ind, " row")
      another.df <- another.df[-row.ind, ]
    }
  }
  return(another.df)
}

filtered.data.try <- filterOutGenesOverQuartile(normalized.counts.prop.uqua)
