## functions heatmap


SubsetCountsByDesignMatrix <- function(counts.df, design.matrix) {
  counts.df <- counts.df[, which(colnames(counts.df) %in% rownames(design.matrix))]
  return(counts.df)
}

CalculateLogFoldChangeForEachTimeNew <- function(subset.counts, design.matrix, log.base=2) {
  times <- unique(design.matrix$Times)
  conditions <- unique(design.matrix$Conditions)
  
  c.fc <- c()
  i <- 1
  for(time in times) {
    #print(time)
    sampl.names <- rownames(design.matrix)[which(design.matrix$Times == time ), drop=FALSE]
    sub.cols <- subset.counts[,which(colnames(subset.counts) %in% sampl.names), drop=FALSE]
    m.cond <- list()
    for(cond in conditions){
      cond.sample.names <- rownames(design.matrix)[which(design.matrix$Conditions == cond ), drop=FALSE]
      sub.sub.cols <- sub.cols[, which(colnames(sub.cols) %in% cond.sample.names), drop=FALSE]
      m.cond[[cond]] <-  apply(sub.sub.cols, 1, mean) 
    }
    fc <- log( (m.cond[[1]]/m.cond[[2]]), base=log.base )
    
    c.fc <- cbind(c.fc, fc)
    colnames(c.fc)[i] <- time
    i<- i+1
  }
  return(c.fc)
} 


OrderByRownames <- function(named.df) {
  return(named.df[order(rownames(named.df)),])
}



PlotQpcrRnaHeatmap <- function(qpcr, rna, threshold=Inf, heatmap.title="heatmap", scale.heat="none", color.palette= c("red", "black", "green"), row.labels.percent=0.9) {
  qpcr[qpcr>threshold] <- threshold
  qpcr[qpcr<(-threshold)] <- (-threshold)
  
  rna[rna>threshold]<-threshold
  rna[rna<(-threshold)]<- (-threshold)
  
  heat.bind<-as.matrix(cbind(qpcr,rna))
  PlotHeatmap(fold.changes = heat.bind, row.labels.percent = row.labels.percent, scale = scale.heat, title = heatmap.title, color.palette = color.palette)
  
}

ComputeCorrelation <- function(qpcr, rna, threshold=Inf){
  qpcr[qpcr>threshold] <- threshold
  qpcr[qpcr<(-threshold)] <- (-threshold)
  
  rna[rna>threshold]<-threshold
  rna[rna<(-threshold)]<- (-threshold)
  
  corr <- list()
  for(i in 1:dim(qpcr)[1]) {
    # print(as.numeric(binded.values[i,]))
    # print(as.numeric(binded.values[i,]))
    corr[[i]]<-cor(as.numeric(qpcr[i,]), as.numeric(rna[i,]))
    # print(corr[[i]])
    # print("---")
  }
  names(corr) <- rownames(qpcr)
  
  corr.unl <- unlist(corr)
  cat("---\nNAs:\n ")
  print( corr.unl[is.na(corr.unl)] )
  cat("---\ncorrelatations:\n")
  print( corr.unl )
  cat("---\nSummary:\n")
  print(summary(corr.unl))
  # boxplot(corr.unl)
  uncorr.genes <- corr.unl[which(corr.unl < 0)]
  return(uncorr.genes)
}

PlotHeatmap <- function(fold.changes, scale="row", dendrogram="none", 
                title="heatmap", row.labels.percent=0.5, save.png=FALSE, 
                file.name=NULL, color.palette = c("red", "black", "green")) 
{
  if(save.png ) {
    if(!is.null(file.name)){
      # pdf(file.name,    # create PNG for the heat map        
      #     width = 5*300,        # 5 x 300 pixels
      #     height = 5*300)       
      png(paste0(file.name,".png"),    # create PNG for the heat map
          width = 5*300,        # 5 x 300 pixels
          height = 5*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8)
    } else {
      error("Please set a file.name for the png file!")
    }
  }
  require("gplots")
  pal <- colorRampPalette(color.palette)(n = 1000)
  # pal <- colorRampPalette(c("yellow","red"))(n = 100)
  lmat = rbind(c(0,3,0), c(2,1,0), c(0,4,0))
  lwid = c(0.5, 4, 0.5)
  lhei = c(1, 4, 1)
  
  heatmap.2(
    fold.changes
    , col = pal
    , scale=scale, 
    dendrogram = dendrogram, ## manage dendrogram
    density.info = "none", ## manage density in legend
    trace="none", ## manage trace along columns
    # margins=c(10, 10), ## manage margins
    Colv=NA, ## manage clustering on columns
    main=title
     , lmat=lmat
     , lhei=lhei
     , lwid = lwid
    , cexRow = row.labels.percent
    , colsep = c(4,8)
    # ,key.par = list(cex=0.5)
    # ,keysize = 0.3
    # , breaks = 5000
    
    # , key=FALSE
  )
  
  if(save.png) {
    dev.off()
  }
}

ComputeZscore <- function(log.values, byRowCol = c(1,2)) {
  log.values.z <-  apply(log.values, 1, function(x) 
  { 
    (x-mean(x))/sd(x)
  })
  return(log.values.z)
}

ComputeLogFCZscores <- function(counts.data.frame, design.matrix, log.base=2, ZbyRowCol = 1) {
  nn.counts <- SubsetCountsByDesignMatrix(counts.data.frame, design.matrix)
  
  l.fc <- CalculateLogFoldChangeForEachTimeNew(subset.counts = nn.counts, design.matrix = design.matrix, log.base = log.base)
  
  l.fc.z <- t(ComputeZscore(l.fc, byRowCol = ZbyRowCol))
  return(list("counts"=nn.counts, "lfc"=l.fc, "lfcz"=l.fc.z))
}



CalculateLogFoldChangeOfFoldChangesForEachTime <- function(subset.counts, design.matrix, log.base=2) {
  times <- as.character(unique(design.matrix$Times))
  
  tissues <- as.character(unique(design.matrix$Tissue))
  
  tissues.fcs <- list()
  for(tissue in tissues) {
    design.tissue.matrix <- design.matrix[which(design.matrix$Tissue %in% tissue ), , drop=FALSE]
    tissue.samples <- rownames(design.tissue.matrix)
    tissue.counts <- subset.counts[,which(colnames(subset.counts) %in% tissue.samples), drop=FALSE]
    t.df.fc <- c()
    conditions <- as.character(unique(design.tissue.matrix$Conditions))
    
    times.fc <- list()
    
    for(time in times) {

      sampl.names <- rownames(design.tissue.matrix)[which(design.tissue.matrix$Times == time ), drop=FALSE]
      sub.cols <- tissue.counts[,which(colnames(tissue.counts) %in% sampl.names), drop=FALSE]
      
      m.cond <- list()
      
      for(cond in conditions) {
        cond.sample.names <- rownames(design.tissue.matrix)[which(design.tissue.matrix$Conditions == cond ), drop=FALSE]
        sub.sub.cols <- sub.cols[, which(colnames(sub.cols) %in% cond.sample.names), drop=FALSE]
        m.cond[[cond]] <-  apply(sub.sub.cols, 1, mean) 
      }
      time.fc <- m.cond[[1]]/m.cond[[2]]
      t.df.fc <- cbind(t.df.fc, time.fc)
      
    }
    colnames(t.df.fc) <- times
    tissues.fcs[[tissue]] <- t.df.fc
  }
  
  df.log.fc.fc <- c()
  for(time in times) {
    time.log.fc.fc <- log( (tissues.fcs[[1]][,time] / tissues.fcs[[2]][,time]) , base = log.base)
    df.log.fc.fc <- cbind(df.log.fc.fc, time.log.fc.fc)
  }
  colnames(df.log.fc.fc) <- times
  
  return(df.log.fc.fc)

}




