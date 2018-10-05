# ## boxplot
# 
# filtered.counts <- ReadDataFrameFromTsv("/media/dario/dati/time_course/Thoracic/results/counts/FeatureCounts/proportion/uqua/thoracic__filtered_counts__Proportion__normalized__uqua.tsv")
# design.matrix <- ReadDataFrameFromTsv("/media/dario/dati/time_course/Thoracic/design_file/thoracic_design_file.txt")

PrepareDataFrameForGGBoxplot <- function(data.frame.to.plot, design.matrix, to.log=TRUE) {
  ## control if design matrix contains times and conditions
  
  if(to.log){
    new.df <- stack(log(data.frame.to.plot+1))
  } else{
    new.df <- stack(data.frame.to.plot)
  }
  
  times.p.s <- c(rep(NA, dim(new.df)[1]))
  conditions.p.s <- c(rep(NA, dim(new.df)[1]))
  new.df <- cbind(new.df, times.p.s, conditions.p.s)
  colnames(new.df) <- c("values", "samples", "times", "conditions")
  head(new.df)
  
  times <- unique(design.matrix$Times)
  for(time in times) {
    samples.t <- rownames(design.matrix)[which(design.matrix$Times %in% time)]
    new.df$times[which(new.df$samples %in% samples.t)] <- as.character(time)
  }
  new.df$times <- as.factor(new.df$times)
  rm(times)
  
  conditions <- unique(design.matrix$Conditions)
  for(condition in conditions) {
    samples.t <- rownames(design.matrix)[which(design.matrix$Conditions %in% condition)]
    new.df$conditions[which(new.df$samples %in% samples.t)] <- as.character(condition)
  }
  new.df$conditions <- as.integer(new.df$conditions)
  rm(conditions)
  return(new.df)
}

PlotTimesBoxplot <- function(data.frame.to.plot, design.matrix, output.path=NULL, prefix.plot=NULL, show.plot.flag=TRUE, save.plot=FALSE, plotly.flag=FALSE) {
  ## control if design matrix contains times and conditions
  require("plotly")
  strings <- GeneratePlotStrings(path = output.path, prefix = prefix.plot, plot.type = "Boxplot")
  
  new.df <- PrepareDataFrameForGGBoxplot(data.frame.to.plot = data.frame.to.plot, design.matrix = design.matrix)
  ggp <- ggplot(new.df, aes(x=samples, y=values, fill=times )) + geom_boxplot(position=position_dodge(2)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(strings$title)
  
  if(save.plot) {
    if(is.null(output.path)) {
      stop("Please set a folder where to plot the boxplot!")
    }
    if(!is.null(strings$filename)){
      SaveGGplot(ggplot.to.save = ggp, plot.folder = strings$plot.folder, plot.file.name = strings$plot.file.name, plotly.flag = plotly.flag)
    }
    
  } 
  
  if(show.plot.flag) {
    if(plotly.flag) {
      ggplotly(ggp)
    } else {
      plot(ggp)
    }
  }
  
    # if(!(is.null(output.path) || is.null(prefix.plot))){
  #   require("htmlwidgets")
  #   filenamepath <- file.path(output.path, paste0(prefix.plot, "_boxplot.html"))
  #   saveWidget(ggplotly(ggp), filenamepath)
  # }
  
  
  ggplotly(ggp)
}


# PrepareDataFrameForGGPCA <- function(data.frame.to.plot, design.matrix, to.log=TRUE) {
#   ## control if design matrix contains times and conditions
#   
#   
#   new.df <- as.data.frame(stack(data.frame.to.plot))
#   
#   new.df <- new.df[,-3]
#   
#   times.p.s <- c(rep(NA, dim(new.df)[1]))
#   conditions.p.s <- c(rep(NA, dim(new.df)[1]))
#   new.df <- cbind(new.df, times.p.s, conditions.p.s)
#   colnames(new.df) <- c("samples", "PCA", "values", "times", "conditions")
#   head(new.df)
#   
#   times <- unique(design.matrix$Times)
#   for(time in times) {
#     samples.t <- rownames(design.matrix)[which(design.matrix$Times %in% time)]
#     new.df$times[which(new.df$samples %in% samples.t)] <- as.character(time)
#   }
#   new.df$times <- as.factor(new.df$times)
#   rm(times)
#   
#   conditions <- unique(design.matrix$Conditions)
#   for(condition in conditions) {
#     samples.t <- rownames(design.matrix)[which(design.matrix$Conditions %in% condition)]
#     new.df$conditions[which(new.df$samples %in% samples.t)] <- as.character(condition)
#   }
#   new.df$conditions <- as.factor(new.df$conditions)
#   rm(conditions)
#   return(new.df)
# }

