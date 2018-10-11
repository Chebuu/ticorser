##################### GENERIC FUNCTIONS
# UpdatePrefix <- function(prefix, postix, sep=" ") {
#     new.prefix <- paste(prefix, postix, sep=sep)
#     return(new.prefix)
# }

SaveGGplot <- function(ggplot.to.save, plot.folder, 
                       plot.file.name, plotly.flag=FALSE) {
    #' Save a plot generated with GGPlot on disk.
    #' 
    #' @param ggplot.to.save a ggplot or ggplotly object.
    #' @param plot.folder Character. the folder where to save the plot.
    #' @param plot.file.name Character. The file.name for the plot file.
    #' @param plotly.flag logical. Is it a ggplotly object?
    #' 
    #' @return none
    #' 
    #' @export 
    #' 
    #' @importFrom grDevices pdf dev
    #' @import htmlwidgets
    
    require("htmlwidgets")
    plotname <- file.path(plot.folder, paste0(plot.file.name, ".pdf"))
    i<-1
    while(file.exists(plotname)) {
        plotname <- file.path(plot.folder, paste0(plot.file.name, i,".pdf"))
        i<-i+1
    }
    
    if(plotly.flag) {
        htmlwidgets::saveWidget(ggplotly(ggplot.to.save), 
                                paste0(getwd(),"/",plotname, ".html"))
                                #paste0(getwd(), "/", plotname, ".html"))
        message("Plot saved as ", paste0(plotname, ".html"))
    } else {
        pdf(plotname)
        # plot(ggplot.to.save)
        print(ggplot.to.save)
        dev.off()
        message("Plot saved as ", plotname)
    }
}


########################################################## DE FUNCTIONS
# ProcessDEResultsForPlot <- function(de.results, threshold) {
#     de.results.new <- de.results
#     if("baseMean" %in% colnames(de.results.new)) { ## working on DESeq results
#         de.results.new <- de.results.new[, c("baseMean", "log2FoldChange", "padj") ]
#         # de.results.new$log2Counts <- log2(de.results.new$baseMean+1)
#         de.results.new$minuslog10PAdj <- (-1)*log10(de.results.new$padj)
#         de.results.new$significance <- "not-significative"
#         de.results.new$significance[which(de.results.new$padj<threshold)] <- "significative"
#         de.results.new$gene <- rownames(de.results.new)
#         de.results.new <- de.results.new[order(de.results.new$padj, decreasing=FALSE),]
#     } else if("prob" %in% colnames(de.results.new)) { ## working on NOISeq results
#         
#         de.results.new <- de.results.new[, "prob", drop=FALSE]
#         de.results.new$log2FoldChange <- de.results.new$log2FC
#         de.results.new$log2Counts <- (1/2) * log2( (de.results[,1] * de.results[,2]) )
#         de.results.new$minuslog10PAdj <- (-1)*log10(de.results.new$padj)
#         de.results.new$significance <- "not-significative"
#         de.results.new$significance[which(de.results.new$prob > threshold)] <- "significative"
#         
#     }
#     return(de.results.new)
# }
# 
# ProcessDEResultsForMAPlot <- function(de.results, counts.dataframe=NULL, design.matrix=NULL, threshold) {
#     de.results.new <- de.results
#     
#     if("baseMean" %in% colnames(de.results.new)) { ## working on DESeq results
#         
#         conds <- design.matrix$Conditions
#         de.results.new <- de.results.new[, c("baseMean", "log2FoldChange", "padj") ]
#         de.results.new$meanA <- apply(counts.dataframe[, which(colnames(counts.dataframe) %in% rownames(subset(design.matrix, Conditions %in% conds[1])) ) ], 1, mean)
#         de.results.new$meanB <- apply(counts.dataframe[, which(colnames(counts.dataframe) %in% rownames(subset(design.matrix, Conditions %in% conds[2])) ) ], 1, mean)
#         de.results.new$log2Counts <- (1/2)*log2((de.results.new$meanA * de.results.new$meanB))
#         de.results.new$significance <- "not-significative"
#         de.results.new$significance[which(de.results.new$padj < threshold)] <- "significative"
#         de.results.new <- de.results.new[order(de.results.new$padj, decreasing=FALSE),]
#         
#     } else if("prob" %in% colnames(de.results.new)) { ## working on NOISeq results
#         
#         log10(results_noi[,2]/results_noi[,1]), - log10(1 - results_noi$prob + 0.000001)
#         
#         de.results.new <- de.results.new[, "prob", drop=FALSE]
#         de.results.new$log2FoldChange <- de.results.new$log2FC
#         de.results.new$log2Counts <- (1/2) * log2( (de.results[,1] * de.results[,2]) )
#         de.results.new$significance <- "not-significative"
#         de.results.new$significance[which(de.results.new$prob > threshold)] <- "significative"
#         
#     }
#     
#     de.results.new$gene <- rownames(de.results.new)
#     
#     return(de.results.new)
# }

################ MA PLOT FUNCTIONS


GenerateGGMA <- function(processed.de.results, strings, plotly.flag=FALSE) {
    #' Generate an MA-Plot using ggplot
    #' 
    #' @param processed.de.results dataframe. D-E results processed with ProcessDEResultsForPlot
    #' @param strings character list. Returned by GeneratePlotStrings
    #' @param plotly.flag logical. Is it a ggplotly object?
    #' 
    #' @return none
    #' 
    #' @export 
    #' 
    #' @importFrom ggplot2 geom_point aes scale_color_manual labs
    #' @import ggplot2
    
    if(plotly.flag) {
        xlabl <- paste("log<sub>2</sub>(counts)")
        ylabl <- paste("-log<sub>2</sub>(FC)")
    }else {
        xlabl <- bquote(~log[2]~"(counts)")
        ylabl <- bquote(~log[2]~"(FC)")
    }
    
    switch(processed.de.results$method[1],
                 DESeq={
                     ggp <- ggplot2::ggplot(processed.de.results) + 
                            geom_point(aes(x=log2Counts, y=log2FoldChange, 
                                           color=significance, ensembl=gene, 
                                           symbol=symbol, 
                                           padj=format(padj, nsmall=10)), 
                                           size=0.7)    
                            + labs( list(title=strings$title, x=xlabl, y=ylabl))
                            + scale_color_manual(values=c("blue2", "red2")) 

                     if(!plotly.flag) {
                         ggp <- ggp + 
                                geom_point(data=subset(processed.de.results, 
                                             significance=="significative"), 
                                             aes(x=log2Counts, y=log2FoldChange, 
                                             color=significance, ensembl=gene, 
                                             symbol=symbol, 
                                             padj=format(padj, nsmall=10)),
                                             size=0.7 ) 
                     }
                 },
                 NOISeqBio={
                     ggp <- ggplot2::ggplot(processed.de.results)
                            + geom_point(aes(x=log2Counts, y=log2FoldChange, 
                                         color=significance, ensembl=gene, 
                                         symbol=symbol, prob=prob), size=0.7)
                            + labs( list(title=strings$title, x=xlabl, y=ylabl))
                            + scale_color_manual(values=c("blue2", "red2")) 
                     if(!plotly.flag) {
                         ggp <- ggp 
                                + geom_point(data=subset(processed.de.results, 
                                             significance=="significative"), 
                                             aes(x=log2Counts, y=log2FoldChange,
                                                 color=significance, 
                                                 ensembl=gene, 
                                                 symbol=symbol, 
                                                 prob=prob), 
                                             size=0.7 ) 
                     }     
                 },
                 NOISeq={
                     ggp <- ggplot2::ggplot(processed.de.results) 
                            + geom_point(aes(x=log2Counts, y=log2FoldChange, 
                                             color=significance, ensembl=gene, 
                                             symbol=symbol, prob=prob), 
                                             size=0.7)
                            + labs( list(title=strings$title, x=xlabl, y=ylabl))
                            + scale_color_manual(values=c("blue2", "red2")) 
                     if(!plotly.flag) {
                         ggp <- ggp 
                                + geom_point(data=subset(processed.de.results, 
                                             significance=="significative"), 
                                             aes(x=log2Counts, y=log2FoldChange,
                                                 color=significance, 
                                                 ensembl=gene, symbol=symbol, 
                                                 prob=prob), size=0.7 ) 
                     }     
                 }
                 
    )
    ggp <- ggp 
           + geom_hline(yintercept=0) 
           + geom_hline(yintercept=1, colour="darkgreen", linetype="dashed") 
           + geom_hline(yintercept=-1, colour="darkgreen", linetype="dashed")
    
    return(ggp)
}


PlotMAPlotCounts <- function(de.results, counts.dataframe, design.matrix, show.plot.flag=TRUE, plotly.flag=FALSE, save.plot=FALSE, plot.folder=NULL, prefix.plot="MA_Plot", threshold) {
    require("plotly")
    
    # title <- paste0(prefix.plot, " MA-Plot")
    
    strings <- GeneratePlotStrings(path=plot.folder, prefix=prefix.plot, plot.type="MAPlot")
    #de.results=de.results; threshold=threshold; counts.dataframe=counts.dataframe; design.matrix=design.matrix
    processed.de.results <- ProcessDEResultsForPlot(de.results, threshold=threshold, counts.dataframe=counts.dataframe, design.matrix=design.matrix)
    # if(plotly.flag) {
    #     ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2Counts, y=log2FoldChange, color=significance, name=gene, padj=padj), size=0.7)    + ggtitle(strings$title) + scale_color_manual(values=c("blue2", "red2")) 
    # } else {
    #     ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2Counts, y=log2FoldChange, color=significance, name=gene, padj=padj), size=0.7) + scale_color_manual(values=c("blue2", "red2"))    + geom_point(data= subset(processed.de.results, significance=="significative"), aes(x=log2Counts, y=log2FoldChange, color=significance, name=gene, padj=padj), size=0.7 )    + ggtitle(strings$title) 
    # }
    ggp <- GenerateGGMA(processed.de.results=processed.de.results, strings=strings, plotly.flag=plotly.flag)
    
    if(save.plot) {
        if(is.null(plot.folder)) {
            stop("Please set a folder where to plot the MA-Plot!")
        }
        if(!is.null(strings$plot.file.name)){
            SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=strings$plot.file.name, plotly.flag=plotly.flag)
        }
        
    } 
    
    if(show.plot.flag) {
        if(plotly.flag) {
            ggplotly(ggp)
        } else {
            plot(ggp)
        }
    }
    
    ##return(ggp)
}

# 
# PlotMAPlot <- function(de.results, show.plot.flag=TRUE, plotly.flag=FALSE, save.plot=FALSE, plot.folder=NULL, prefix.plot="MA_Plot", threshold) {
#     require("plotly")
#     
#     # title <- paste0(prefix.plot, " MA-Plot")
#     
#     strings <- GeneratePlotStrings(path=plot.folder, prefix=prefix.plot, plot.type="MAPlot")
#     
#     processed.de.results <- ProcessDEResultsForPlot(de.results, threshold)
#     if(plotly.flag) {
#         ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2Counts, y=log2FoldChange, color=significance, name=gene, padj=padj), size=0.7) + labs(list(title=strings$title, x="log'[2](Counts)", y="log'[2](FC)")) + scale_color_manual(values=c("blue2", "red2")) 
#     } else {
#         ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2Counts, y=log2FoldChange, color=significance, name=gene, padj=padj), size=0.7) + scale_color_manual(values=c("blue2", "red2"))    + geom_point(data= subset(processed.de.results, significance=="significative"), aes(x=log2Counts, y=log2FoldChange, color=significance, name=gene, padj=padj), size=0.7 )    + labs(list(title=strings$title, x="log'[2](Counts)", y="log'[2](FC)")) 
#     }
#     
#     
#     if(save.plot) {
#         if(is.null(plot.folder)) {
#             stop("Please set a folder where to plot the MA-Plot!")
#         }
#         if(!is.null(strings$plot.file.name)){
#             SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=strings$plot.file.name, plotly.flag=plotly.flag)
#         }
#         
#     } 
#     
#     if(show.plot.flag) {
#         if(plotly.flag) {
#             ggplotly(ggp)
#         } else {
#             plot(ggp)
#         }
#     }
#     
#     ##return(ggp)
# }

#########################
ProcessDEResultsForPlot <- function(de.results, threshold, counts.dataframe=NULL, design.matrix=NULL) {
    de.results.new <- de.results
    de.results.new <- de.results.new[order(rownames(de.results.new)),]
    counts.dataframe.ord <- counts.dataframe[order(rownames(counts.dataframe)),]
    
    if("baseMean" %in% colnames(de.results.new)) { ## working on DESeq results
        
        de.results.new <- de.results.new[, c("padj", "log2FoldChange", "symbol"), drop=FALSE ]
        
        conds <- as.character(unique(design.matrix$Conditions))
        sub.des.A <- subset(design.matrix, design.matrix$Conditions %in% conds[1])
        sub.counts.A <- counts.dataframe.ord[, which(colnames(counts.dataframe.ord) %in% rownames(sub.des.A) ), drop=FALSE ]
        sub.des.B <- subset(design.matrix, design.matrix$Conditions %in% conds[2])
        sub.counts.B <- counts.dataframe.ord[, which(colnames(counts.dataframe.ord) %in% rownames(sub.des.B) ), drop=FALSE ]
        de.results.new$meanA <- apply(sub.counts.A, 1, mean)
        de.results.new$meanB <- apply(sub.counts.B, 1, mean)
        #########
        de.results.new$log2Counts <- (1/2) * log2((de.results.new$meanA * de.results.new$meanB))
        de.results.new$log10FoldChange <- log10( (de.results.new$meanA / de.results.new$meanB ) )
        de.results.new$minuslog10PAdj <- (-1) * log10(de.results.new$padj)
        de.results.new$significance <- "not-significative"
        de.results.new$significance[which(de.results.new$padj < threshold)] <- "significative"
        de.results.new <- de.results.new[order(de.results.new$padj, decreasing=FALSE),]
        de.results.new$method <- rep(x="DESeq", times=dim(de.results.new)[1])
        
    } else if("theta" %in% colnames(de.results.new)) { ## working on NOISeqbio results
        
        de.results.new <- de.results.new[, c(1, 2, 6)]
        de.results.new$prob <- de.results$prob
        de.results.new$log2FoldChange <- de.results$log2FC
        de.results.new$log10FoldChange <- log10( (de.results[,1]/de.results[,2]) )
        de.results.new$log2Counts <- (1/2) * log2( (de.results[,1] * de.results[,2]) )
        de.results.new$significance <- "not-significative"
        de.results.new$significance[which(de.results.new$prob > threshold)] <- "significative"
        de.results.new <- de.results.new[order(de.results.new$prob, decreasing=TRUE),]
        de.results.new$minuslog101minuspp <- (-1) * log10( (1 - de.results.new$prob + 0.000001))
        de.results.new$method <- rep(x="NOISeqBio", times=dim(de.results.new)[1])
        
    } else if("M" %in% colnames(de.results.new)) { ## working on NOISeq results
        
        de.results.new <- de.results.new[, c(1, 2, 7)]
        de.results.new$prob <- de.results$prob
        de.results.new$log2FoldChange <- de.results$M
        de.results.new$log10FoldChange <- log10( (de.results[,1]/de.results[,2]) )
        de.results.new$log2Counts <- (1/2) * log2( (de.results[,1] * de.results[,2]) )
        de.results.new$significance <- "not-significative"
        de.results.new$significance[which(de.results.new$prob > threshold)] <- "significative"
        de.results.new <- de.results.new[order(de.results.new$prob, decreasing=TRUE),]
        de.results.new$minuslog101minuspp <- (-1) * log10( (1 - de.results.new$prob + 0.000001))
        de.results.new$method <- rep(x="NOISeq", times=dim(de.results.new)[1])
        
    }
    
    de.results.new$gene <- rownames(de.results.new)
    
    return(de.results.new)
}



############ VOLCANO FUNCTIONS 


GenerateGGVolcano <- function(processed.de.results, strings, plotly.flag) {
    require(ggplot2)
    switch(processed.de.results$method[1],
                 DESeq={
                     if(plotly.flag) {
                         xlabl <- paste0("log<sub>2</sub>(FC)")
                         ylabl <- paste("-log<sub>10</sub>(padj)")
                     }else {
                         xlabl <- bquote(~log[2]~"(FC)")
                         ylabl <- bquote(~-log[10]~"(padj)")
                     }
                     
                     ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2FoldChange, y=minuslog10PAdj, color=significance, ensembl=gene, symbol=symbol, padj=format(padj, nsmall=10) ), size=0.7) + labs(list(title=strings$title, x=xlabl, y=ylabl)) + scale_color_manual(values=c("blue2", "red2"))
                     if(!plotly.flag) {
                         ggp <- ggp + geom_point(data= subset(processed.de.results, significance=="significative"), aes(x=log2FoldChange, y=minuslog10PAdj, color=significance, ensembl=gene, symbol=symbol, padj=format(padj, nsmall=10)), size=0.7 )
                     }     
                 },
                 NOISeqBio={
                     if(plotly.flag) {
                         xlabl <- paste0("log<sub>2</sub>(FC)")
                         ylabl <- paste0("-log<sub>10</sub>(1-prob)")
                     }else {
                         xlabl <- bquote(~log[2]~"(FC)")
                         ylabl <- bquote(~-log[10]~"(1-Prob)")
                     }
                     
                     ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2FoldChange, y=minuslog101minuspp, color=significance, ensembl=gene, symbol=symbol, prob=format(prob, nsmall=10)), size=0.7) + labs(list(title=strings$title, x=xlabl, y=ylabl)) + scale_color_manual(values=c("blue2", "red2")) 
                     if(!plotly.flag) {
                         ggp <- ggp + geom_point(data= subset(processed.de.results, significance=="significative"), aes(x=log2FoldChange, y=minuslog101minuspp, color=significance, ensembl=gene, symbol=symbol, prob=format(prob, nsmall=10)), size=0.7 )
                     }     
                 },
                 NOISeq={
                     if(plotly.flag) {
                         xlabl <- paste0("log<sub>2</sub>(FC)")
                         ylabl <- paste0("-log<sub>10</sub>(1-prob)")
                     }else {
                         xlabl <- bquote(~log[2]~"(FC)")
                         ylabl <- bquote(~-log[10]~"(1-Prob)")
                     }
                     
                     ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2FoldChange, y=minuslog101minuspp, color=significance, ensembl=gene, symbol=symbol, prob=format(prob, nsmall=10)), size=0.7) + labs(list(title=strings$title, x=xlabl, y=ylabl)) + scale_color_manual(values=c("blue2", "red2")) 
                     if(!plotly.flag) {
                         ggp <- ggp + geom_point(data= subset(processed.de.results, significance=="significative"), aes(x=log2FoldChange, y=minuslog101minuspp, color=significance, ensembl=gene, symbol=symbol, prob=format(prob, nsmall=10)), size=0.7 )
                     }     
                 }
    )
    
    ggp <- ggp + geom_vline(xintercept=0) + geom_vline(xintercept=1, colour="darkgreen", linetype="dashed") + geom_vline(xintercept=-1, colour="darkgreen", linetype="dashed")
    
    return(ggp)
}

PlotVolcanoPlot <- function(de.results, counts.dataframe=NULL, design.matrix=NULL, show.plot.flag=TRUE, plotly.flag=FALSE, save.plot=FALSE, plot.folder=NULL, prefix.plot=NULL, threshold) {
    require("plotly")
    # title <- paste0(prefix.plot, " Volcano Plot")
    strings <- GeneratePlotStrings(path=plot.folder, prefix=prefix.plot, plot.type="VolcanoPlot")
    
    processed.de.results <- ProcessDEResultsForPlot(de.results=de.results, threshold=threshold, counts.dataframe=counts.dataframe, design.matrix=design.matrix)
    
    ggp <- GenerateGGVolcano(processed.de.results, strings, plotly.flag)
    
    
    if(save.plot) {
        if(is.null(plot.folder)) {
            stop("Please set a folder where to plot the Volcano-Plot!")
        }
        if(!is.null(strings$plot.file.name)){
            SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=strings$plot.file.name, plotly.flag=plotly.flag)
        } else {
            stop("Please set a name for the Volcano-Plot!")
        }
        
        
    } 
    
    if(show.plot.flag) {
        if(plotly.flag) {
            ggplotly(ggp)
        } else {
            plot(ggp)
        }
    }
}


###### PCA FUNCTIONS
PlotPCAPlotlyFunction <- function(counts.data.frame, design.matrix, scale=FALSE, plot.folder=NULL, prefix.plot=NULL, show.plot.flag=TRUE, plotly.flag=FALSE, save.plot=FALSE, ellipse.flag=FALSE) {
    
    # counts.data.frame=cervical.wilcoxon.unnormalized; design.matrix=cervical.design.matrix; colour.design.column.str="Times"; shape.design.column.str="Conditions"; 
    # output.path=cervical.output.plots.folder; prefix.plot="cervical_unnormalized"; title="Cervical Un-Normalized PCA"
    
    ## check colnames design matrix
    # require("ggfortify")
    require("plotly")
    strings <- GeneratePlotStrings(path=plot.folder, prefix=prefix.plot, plot.type="PCA")
    sub.counts.dataframe <- counts.data.frame[ , which(colnames(counts.data.frame) %in% rownames(design.matrix)), drop=FALSE]
    PCA <- prcomp(t(sub.counts.dataframe), scale.=scale, center=TRUE)
    
    sub.pca <- as.data.frame(PCA$x[,c("PC1", "PC2")])
    sub.pca$samples <- rownames(sub.pca)
    sub.pca <- sub.pca[order(sub.pca$samples), ]
    design.matrix.o <- design.matrix[order(rownames(design.matrix)), , drop=FALSE]
    
    if(dim(sub.pca)[1] == dim(design.matrix.o)[1]) {
        new.sub.pca <- cbind(sub.pca, design.matrix.o)
    } else {
        sub.design.matrix <- subset(design.matrix.o, rownames(design.matrix.o) %in% rownames(sub.pca))
        new.sub.pca <- cbind(sub.pca, sub.design.matrix)
    }
    
    if(ellipse.flag=="both") {
        PlotPCAPlotlyFunction(counts.data.frame, design.matrix, scale, plot.folder, prefix.plot, show.plot.flag, plotly.flag, save.plot, ellipse.flag=FALSE)
        PlotPCAPlotlyFunction(counts.data.frame, design.matrix, scale, plot.folder, prefix.plot, show.plot.flag, plotly.flag, save.plot, ellipse.flag=TRUE)
        return()
    } else if(ellipse.flag) { 
        if(length(which(colnames(new.sub.pca) %in% "Times")) >0 ) {
            ggp <- ggplot2::ggplot(new.sub.pca) + geom_point(aes(x=PC1, y=PC2, color=Times, shape=Conditions, name=samples), size=2) + stat_ellipse(aes(x=PC1, y=PC2, color=Conditions))+ ggtitle(strings$title)
        } else {
            ggp <- ggplot2::ggplot(new.sub.pca) + geom_point(aes(x=PC1, y=PC2, color=Conditions, shape=Conditions, name=samples), size=2) + stat_ellipse(aes(x=PC1, y=PC2, color=Conditions))+ ggtitle(strings$title)
        }
        
        strings$plot.file.name <- paste0(strings$plot.file.name, "_ellipse")
    } else if(!ellipse.flag) {
        if(length(which(colnames(new.sub.pca) %in% "Times")) >0 ) {
        ggp <- ggplot2::ggplot(new.sub.pca) + geom_point(aes(x=PC1, y=PC2, color=Times, shape=Conditions, name=samples), size=2) + ggtitle(strings$title)
        } else {
            ggp <- ggplot2::ggplot(new.sub.pca) + geom_point(aes(x=PC1, y=PC2, color=Conditions, shape=Conditions, name=samples), size=2) + stat_ellipse(aes(x=PC1, y=PC2, color=Conditions))+ ggtitle(strings$title)
        }
    }
    
    
    if(save.plot) {
        if(is.null(plot.folder)) {
            stop("Please set a folder where to plot the PCA!")
        }
        if(!is.null(strings$plot.file.name)) {
            SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=strings$plot.file.name, plotly.flag=plotly.flag)
        }
        
    } 
    
    if(show.plot.flag) {
        if(plotly.flag) {
            ggplotly(ggp)
        } else {
            plot(ggp)
        }
    }
    
}


# 
# PlotPCAFunctionBu <- function(samples.dataframe, design.dataframe, plot.label="PCA", save.plot=FALSE) {
#     
#     samples.dataframe.reo <- samples.dataframe[, order(colnames(samples.dataframe))]
#     design.dataframe.reo <- design.dataframe[order(rownames(design.dataframe)),]
#     
#     PCA=prcomp(t(samples.dataframe.reo), scale.=TRUE, center=TRUE)
#     
#     times <- unique(design.dataframe.reo$Times)
#     rand.cols <- rainbow(length(times))
#     
#     colors <- c()
#     pch.c <- c()
#     i <- 1
#     j <- 1
#     for (time in times) {
#         # print(time)
#         ind.row.t <- which(design.dataframe.reo$Times %in% time)
#         conditions <- unique(design.dataframe.reo$Conditions[ind.row.t])
#         colors <- c(colors, rep(rand.cols[i], length(ind.row.t)))
#         
#         # print(colors)
#         for(condition in conditions) {
#             ind.row.c <- which(design.dataframe.reo$Conditions[ind.row.t] %in% condition)
#             for(ind in ind.row.c) {
#                 pch.c <- c(pch.c, j)
#                 j <- j+1
#             }
#             
#             # if(j == 2) {
#             #     j <- 1
#             # } else {
#             #     j <- 2
#             # }
#             
#             # print(pch.c)
#         }
#         j<-1
#         i <- i+1
#     }
#     
#     par(xpd=T, mar=par()$mar + c(0,0,0,7))
#     plot(PCA$x, pch=pch.c, col=colors, cex=1.5, main=plot.label, lwd=1.5)
#     legend(max(PCA$x[,1])+10, max(PCA$x[,2])+10, cex=0.8 , legend=colnames(samples.dataframe.reo), pch=pch.c, col=colors)## trovare il massimo sulla y per la legend
#     ## implementare il salvataggio del plot
#     par(mar=c(5, 4, 4, 2) + 0.1)
#     
# }
# 
# PlotPCAFunction <- function(samples.dataframe, design.dataframe, scale=TRUE, plot.label="PCA", save.plot=FALSE) {
#     
#     # samples.dataframe.reo <- samples.dataframe[, order(colnames(samples.dataframe))]
#     # design.dataframe.reo <- design.dataframe[order(rownames(design.dataframe)),]
#     
#     PCA=prcomp(t(samples.dataframe), scale.=scale, center=TRUE)
#     
#     times <- unique(design.dataframe$Times)
#     rand.cols <- rainbow(length(times))
#     
#     colors <- c()
#     pch.c <- c()
#     i <- 1
#     
#     for (time in times) {
#         # print(time)
#         ind.row.t <- which(design.dataframe$Times %in% time)
#         colors <- c(colors, rep(rand.cols[i], length(ind.row.t)))
#         
#         # # print(colors)
#         # for(condition in conditions) {
#         #     ind.row.c <- which(design.dataframe.reo$Conditions[ind.row.t] %in% condition)
#         #     pch.c <- c(pch.c, rep(j, length(ind.row.c)))
#         #     if(j == 2) {
#         #         j <- 1
#         #     } else {
#         #         j <- 2
#         #     }
#         #     
#         #     # print(pch.c)
#         # }
#         i <- i+1
#     }
#     # print(colors)
#     conditions <- unique(design.dataframe$Conditions)
#     j <- 1
#     for(condition in conditions) {
#         ind.row.c <- which(design.dataframe$Conditions %in% condition)
#         pch.c <- c(pch.c, rep(j, length(ind.row.c)))
#         if(j == 2) {
#             j <- 1
#         } else {
#             j <- 2
#         }
#         
#         # print(pch.c)
#     }
#     
#     par(xpd=T, mar=par()$mar + c(0,0,0,7))
#     plot(PCA$x, pch=pch.c, col=colors, cex=1.5, main=plot.label, lwd=1.5)
#     legend(max(PCA$x[,1])+50, max(PCA$x[,2])+50, cex=0.8 , legend=colnames(samples.dataframe), pch=pch.c, col=colors)## trovare il massimo sulla y per la legend
#     ## implementare il salvataggio del plot
#     par(mar=c(5, 4, 4, 2) + 0.1)
#     
# }

PlotVariancePCA <- function(samples.dataframe, plot.label="PCA Variability", save.plot=FALSE) {
    screeplot(princomp(samples.dataframe), main=plot.label)
}


############## BOXPLOT FUNCTIONS
PrepareDataFrameForGGBoxplot <- function(data.frame.to.plot, design.matrix, to.log=TRUE) {
    ## control if design matrix contains times and conditions
    
    if(to.log){
        new.df <- as.data.frame(stack(log(data.frame.to.plot+1)))
    } else{
        new.df <- as.data.frame(stack(data.frame.to.plot))
    }
    new.df <- new.df[,c(1:2,4)]
    times.p.s <- c(rep(NA, dim(new.df)[1]))
    conditions.p.s <- c(rep(NA, dim(new.df)[1]))
    new.df <- cbind(new.df, times.p.s, conditions.p.s)
    colnames(new.df) <- c("names", "samples", "values", "times", "conditions")
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
    # new.df$conditions <- as.integer(new.df$conditions)
    rm(conditions)
    return(new.df)
}

    

PlotTimesBoxplot <- function(data.frame.to.plot, design.matrix, output.path=NULL, prefix.plot=NULL, show.plot.flag=TRUE, save.plot=FALSE, plotly.flag=FALSE) {
    ## control if design matrix contains times and conditions
    require("plotly")
    strings <- GeneratePlotStrings(path=output.path, prefix=prefix.plot, plot.type="Boxplot")
    
    new.df <- PrepareDataFrameForGGBoxplot(data.frame.to.plot=data.frame.to.plot, design.matrix=design.matrix)
    
    ggp <- ggplot(new.df, aes(x=samples, y=names, fill=times )) + geom_boxplot(position=position_dodge(2)) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle(strings$title)
    
    
    if(save.plot) {
        if(is.null(strings$plot.folder)) {
            stop("Please set a folder where to plot the boxplot!")
        }
        if(!is.null(strings$plot.file.name)){
            SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=strings$plot.file.name, plotly.flag=plotly.flag)
        }
        
    }
    
    if(show.plot.flag) {
        if(plotly.flag) {
            ggplotly(ggp)
        } else {
            plot(ggp)
        }
    }
    
}

### scatter plot matrix
PlotScatterPlotMatrix <- function(counts.data.frame, design.matrix, plot.folder=NULL, prefix.plot=NULL, show.plot.flag=TRUE, plotly.flag=FALSE, save.plot=FALSE) {

    require("GGally")
    require("plotly")
    strings <- GeneratePlotStrings(path=plot.folder, prefix=prefix.plot, plot.type="ScatterPlotMatrix")
    
    df <- counts.data.frame[, which(colnames(counts.data.frame) %in% rownames(design.matrix)) ]
    
    ggp <- ggpairs(df)
    
    if(save.plot) {
        if(is.null(strings$plot.folder)) {
            stop("Please set a folder where to plot the boxplot!")
        }
        if(!is.null(strings$plot.file.name)){
            SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=strings$plot.file.name, plotly.flag=plotly.flag)
        }
    } 
    
    if(show.plot.flag) {
        if(plotly.flag) {
            ggplotly(ggp)
        } else {
            plot(ggp)
        }
    }
    
    
}



###### plot counts along times
ProcessCountDataFrameForPlotCountsAcrossTimes <- function(gene.name.normalized.counts, design.matrix, gene.name) {
    gene.name.norm.counts.t <- as.data.frame(t(gene.name.normalized.counts))
    
    gene.name.norm.counts.t$GeneName <- gene.name
    gene.name.norm.counts.t$Times <- as.character(design.matrix[which(rownames(design.matrix) %in% rownames(gene.name.norm.counts.t)), "Times"])
    gene.name.norm.counts.t$Conditions <- as.character(design.matrix[which(rownames(design.matrix) %in% rownames(gene.name.norm.counts.t)), "Conditions"])
    gene.name.norm.counts.t$Counts <- gene.name.norm.counts.t[,1]
    return(gene.name.norm.counts.t)
}

PlotCountsAlongTimes <- function(normalized.counts, design.matrix, gene.name, gene.name.column.name="gene.names", show.plot.flag=TRUE, plotly.flag=FALSE, save.plot=FALSE, plot.folder=NULL, prefix.plot="CountsTimePlot"){
    
    strings <- GeneratePlotStrings(path=plot.folder, prefix=prefix.plot, plot.type="CountsTimePlot")
    
    sub.normalized.counts <- normalized.counts[,which(colnames(normalized.counts) %in% rownames(design.matrix))]
    sub.normalized.counts[,gene.name.column.name] <-    normalized.counts[,gene.name.column.name]
    
    gene.name.norm.counts <- sub.normalized.counts[which(tolower(sub.normalized.counts[,gene.name.column.name]) %in% tolower(gene.name)), ]
    
    gene.name.norm.counts <- gene.name.norm.counts[, which(colnames(gene.name.norm.counts) %in% rownames(design.matrix))]
    
    if(dim(gene.name.norm.counts)[1] == 1) {
        processed.counts.df <- ProcessCountDataFrameForPlotCountsAcrossTimes(gene.name.norm.counts, design.matrix, gene.name)
    } else if(dim(gene.name.norm.counts)[1] == 0) {
        stop(gene.name, " gene Not Found!")
    } else if(dim(gene.name.norm.counts)[1] > 1) {
        stop(gene.name, " founded in more than one row!")
    }
    require("plotly")
    ggp <- ggplot(processed.counts.df, mapping=aes(x=Times, y=Counts, color=Conditions, group=Conditions)) + geom_point() + stat_smooth(se=FALSE, method="loess") + scale_y_log10() + ggtitle(paste(strings$title, gene.name, "gene", sep=" "))
    
    if(save.plot) {
        if(is.null(strings$plot.folder)) {
            stop("Please set a folder where to plot the boxplot!")
        }
        if(!is.null(strings$plot.file.name)){
            SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=paste(strings$plot.file.name, gene.name, sep="_"), plotly.flag=plotly.flag)
        }
        
    } 
    
    if(show.plot.flag) {
        if(plotly.flag) {
            ggplotly(ggp)
        } else {
            plot(ggp)
        }
    }
    
}


PlotSummary <- function(summary.df, dot.size=0.3) {
    require("ggplot2")
    p <- ggplot(summary.df, aes(x=times, y=genes, fill=values)) + geom_dotplot(binaxis='y', stackdir='center',    dotsize=dot.size) + scale_fill_manual(values=c("DOWN"="red", "NO_CHANGE"="white", "UP"="green", "LRT"="turquoise", "WALD"="plum1", "LRT+WALD"="orange"))
    plot(p)
}



