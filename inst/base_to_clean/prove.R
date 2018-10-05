source("https://bioconductor.org/biocLite.R")
biocLite("fission")
library("fission")
data("fission")
head(fission)
library("DESeq2")
?DESeqDataSet
head(fission@colData)
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)
str(ddsTC)
ddsTC <- DESeq(ddsTC, test = "LRT", reduced = ~ strain + minute)
#resTC <- results(ddsTC) 
resTC$symbol <- mcols(ddsTC)$symbol

setwd("Dropbox/Lavori/IAC/time_course/")

############# design
times <- c(rep(c("4d", "7d", "14d", "56d"), each=5), rep(c("4d", "7d", "14d", "56d"), each=3))
conditions <- c(rep("Treated", each=5*4), rep("control", each=3*4))
replicates <- c( rep(LETTERS[1:5], times=4), rep(LETTERS[1:3], times=4))
names <- paste(conditions, replicates, times, sep="_")

data.design <- data.frame(Times=times, Conditions=conditions, row.names = names)
head(data.design)
write.table(x=data.design, file = "design_file.txt", sep = "\t")

###########

###### dataset
counts.t1 <- read.table(file = "counts_FeatureCounts_t1.txt", header = TRUE, sep = "\t")
head(counts.t1)
counts.t1 <- cbind(counts.t1, counts.t1[,1:2])

#colnames(counts.t1) <- c("t1_1", "t1_2", "t1_3", "t1_4", "t1_5", "c1_1", "c1_2", "c1_3")

counts <- cbind(counts.t1, counts.t1, counts.t1, counts.t1)
colnames(counts) <- names
head(counts)

counts.dataframe <- counts
design.dataframe <- data.design
  
DeSeqTimeCourse <- function(counts.dataframe, design.dataframe) {
  try( 
    if( (!("Times" %in% colnames(design.dataframe)) ||  !("Conditions" %in% colnames(design.dataframe))))
      stop("Design data frame must contain Times and Conditions colnames") 
    )
  
  require("DESeq2")
  deseq.dataset <- DESeq2::DESeqDataSetFromMatrix(countData = counts.dataframe, colData = design.dataframe, design = ~ Times + Conditions + Times:Conditions)
  de.results <- DESeq2::DESeq(deseq.dataset, test = "LRT", reduced = ~ Times + Conditions)
  readable.results <- DESeq2::results(de.results)
  return(readable.results)
}

# dataset <- DESeqDataSetFromMatrix(countData = counts, colData = data.design, design = ~ Times + Condition + Times:Condition)
# res <- DESeq(dataset, test = "LRT", reduced = ~ Times + Condition)
# res.f <- results(res)

data <- plotCounts(dataset, which.min(res.f$padj), intgroup=c("Times","Condition"), returnData=TRUE)
pp <- ggplot(data, aes(x=Times, y=count, color=Condition, group=Condition)) + geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()

plot_ly(pp)
# counts.t2 <- read.table(file = "counts_FeatureCounts_t2.txt", header = TRUE, sep = "\t")
# head(counts.t2)

FilterLowCounts <- function(counts.dataframe, design.dataframe, is.normalized = c(TRUE, FALSE), method.type = c("CPM", "Wilcoxon", "Proportion"), cv.percentage, cpm.cutoff, seq.depth=NULL) {
  try( 
    if( !("Conditions" %in% colnames(design.dataframe)))
      stop("Design dataframe must contain Conditions on colnames"),
    if( sum(rownames(design.dataframe) %in% colnames(counts.dataframe)) != length(rownames(design.dataframe)) )
      stop("Rownames of Design dataframe must be equal to Colnames of Counts dataframe ")
  )
  # design.dataframe <- design.dataframe[rand,]
  # design.dataframe[which( rownames(design.dataframe) %in% colnames(counts.dataframe)), "Conditions"]
  
  design.dataframe <- design.dataframe[order(rownames(design.dataframe)), ]
  counts.dataframe <- counts.dataframe[, order(colnames(counts.dataframe))]
  
  conditions <- design.dataframe$Conditions
  
  if(method.type == "Proportion") {
    if(is.normalized) {
      ## calcolare seq.depth
      if(is.null(seq.depth)) stop("Proportion test cannot be performed on normalized counts without sequencing depth!")
    } else {
      seq.depth <- NULL
    }
  }
  
  switch( method.type,
          CPM = {
            method.number <- 1
          },
          Wilcoxon = {
            method.number <- 2
          },
          Proportion = {
            method.number <- 3
          }
  )
  
  require("NOISeq")
  filtered.dataframe <- NOISeq::filtered.data(dataset = counts.dataframe, factor = conditions, norm = is.normalized, depth = seq.depth, method = method.number, cv.cutoff = cv.percentage, cpm = cpm.cutoff, p.adj = "BH")
  
  return(filtered.dataframe)
}


biocLite("maSigPro")
library(maSigPro)
data(NBdata)
data("NBdesign")


########################################################## PCA
PlotPCA <- function(samples.dataframe, design.dataframe, plot.label="PCA", save.plot=FALSE, legendpos = "bottomright") {
  
  samples.dataframe.reo <- samples.dataframe[, order(colnames(samples.dataframe))]
  design.dataframe.reo <- design.dataframe[order(rownames(design.dataframe)),]

  PCA=prcomp(t(log(samples.dataframe.reo+1)))
  
  times <- unique(design.dataframe.reo$Times)
  rand.cols <- rainbow(length(times))
  
  colors <- c()
  pch.c <- c()
  i <- 1
  j <- 1
  for (time in times) {
    # print(time)
    ind.row.t <- which(design.dataframe.reo$Times %in% time)
    conditions <- unique(design.dataframe.reo$Conditions[ind.row.t])
    colors <- c(colors, rep(rand.cols[i], length(ind.row.t)))
    # print(colors)
    for(condition in conditions) {
      ind.row.c <- which(design.dataframe.reo$Conditions[ind.row.t] %in% condition)
      pch.c <- c(pch.c, rep(j, length(ind.row.c)))
      j <- j+1
      # print(pch.c)
    }
    i <- i+1
  }
  
  par(xpd = T, mar = par()$mar + c(0,0,0,7))
  plot(PCA$x, pch=pch.c, col=colors, cex=1.5, main = plot.label, lwd=1.5)
  legend(max(PCA$x)+10, 29, cex=0.8 , legend=colnames(samples.dataframe.reo), pch=pch.c, col=colors)## trovare il massimo sulla y per la legend
  ## implementare il salvataggio del plot
  par(mar=c(5, 4, 4, 2) + 0.1)
  
}


PlotLoadingsPCA <- funtion(samples.dataframe, plot.label="PCA Loadings", save.plot=FALSE) {
  
  screeplot(princomp(samples.dataframe), main=plot.label)
}


### filter high counts
FilterHighCounts <- function(samples.dataframe, numerical.threshold=NULL, quantile.threshold=.99) {
  
  q <- c()
  for(j in 1:dim(samples.dataframe)[2]){
    q <- c(q, quantile(samples.dataframe[,j], probs=quantile.threshold))
  }
  print(q)
  q <- mean(q)
  print(q)
  
  mat.ind <- which(samples.dataframe >= q, arr.ind = TRUE)
  mat.ind[,1]
  head(mat.ind)
  sd <- samples.dataframe[-which(rownames(samples.dataframe) %in% rownames(mat.ind)), ]
  
}


PlotPCAPlotlyFunction <- function(counts.data.frame, design.matrix, colour.design.column.str, shape.design.column.str, scale=FALSE, output.path=NULL, prefix.plot=NULL) {
  ## check colnames design matrix
  require("ggfortify")
  require("plotly")
  
  PCA <- prcomp(t(counts.data.frame), scale.=scale, center=TRUE)
  ggpca <- autoplot(PCA, 
                    data=design.matrix, colour=colour.design.column.str, 
                    label=FALSE, #label.size=3,
                    shape=FALSE)#, 
  # shape=as.factor(thoracic.design.matrix$Conditions))#, 
  #loadings=TRUE, loadings.label = TRUE)
  condition <- design.matrix[,shape.design.column.str]
  time <- design.matrix[,colour.design.column.str]
  name <- rownames(design.matrix)
  ggpca.sh <- ggpca + geom_point(aes(shape=condition, colour=time, name=name), size=2)
  if(!(is.null(output.path) && is.null(prefix.plot))){
    require("htmlwidgets")
    filenamepath <- file.path(output.path, paste0(prefix.plot, "_PCA.html"))
    htmlwidgets::saveWidget(ggplotly(ggpca.sh), filanemapath)
  }
  ggplotly(ggpca.sh)
}

# library("ggfortify")
# scale=FALSE
# samples.dataframe.reo <- ruved.thoracic.set@assayData$normalizedCounts
# PCA <- prcomp(t(samples.dataframe.reo), scale.=scale, center=TRUE)
# autoplot(PCA, data=thoracic.design.matrix, colour="Times", label=FALSE, label.size=3, shape=TRUE, loadings=FALSE, shape="Conditions")
# ggpca <- autoplot(PCA, 
#                 data=thoracic.design.matrix, colour="Times", 
#                 label=FALSE, #label.size=3,
#                 shape=FALSE)#, 
#                 # shape=as.factor(thoracic.design.matrix$Conditions))#, 
#                 #loadings=TRUE, loadings.label = TRUE)
# condition <- thoracic.design.matrix$Conditions
# time <- thoracic.design.matrix$Times
# name <- rownames(thoracic.design.matrix)
# 
# ggpca.sh <- ggpca + geom_point(aes(shape=condition, colour=time, name=name), size=2)
# # install.packages("plotly")
# require("plotly")
# 
# # ggpcaly <- plotly_build(ggpca.sh)
# # head(ggpcaly)
# save(ggpca.sh, file = "ggpca.RData")
# ggplotly(ggpca.sh)
# 







