createSignLabelsNumbers <- function(proc.df, column.name)
{
    values <- unique(proc.df[[column.name]])
    tot.de <- sum(proc.df[[column.name]] == values[1])
    tot.not.de <- dim(proc.df)[1] - tot.de 
    
    idxsign <- which(proc.df[[column.name]]==values[1])
    proc.df[[column.name]][idxsign] <- paste0(values[1], 
                                            " [", tot.de, "]")
    proc.df[[column.name]][-idxsign] <- paste0(values[2], 
                                            " [", tot.not.de, "]")
    return(proc.df)
}


luciaVolcanoPlot <- function(res.o, positive.controls.df, prefix)
{
    xlabl <- paste0("log<sub>2</sub>(FC)")
    ylabl <- paste("-log<sub>10</sub>(PValue)")
    title <- paste(prefix, "Volcano Plot")
    
    new.de <- ProcessDEResultsForPlot(de.results=res.o, threshold=0.01, 
                                      design.matrix=desMat)
    pos.contr <- positive.controls.df
    with.pos.de <- new.de
    with.pos.de$padj <- format(round(with.pos.de$padj, 5), nsmall=5)
    without.pos.de <- createSignLabelsNumbers(with.pos.de, "significance")
    pp <- ggplot(data=without.pos.de) +
        geom_point(aes(x=log2FoldChange, y=minuslog10pval, color=significance, 
                        text=paste0("padj=", padj, "\nname=", gene)), 
                        size=0.7) 
    if(!is.null(positive.controls.df) )
    {
        with.pos.de$hit <- NA
        with.pos.de$hit[which(tolower(with.pos.de$gene) %in% 
                                tolower(pos.contr[,1]))] <- "pctr"
        
        sub.de <- with.pos.de[which(with.pos.de$hit=="pctr"),] 
        sub.de$lit <- "est"
        pos.lit <- pos.contr[pos.contr[,2]=="lit",, drop=FALSE]
        idx.lit <- which(tolower(sub.de$gene) %in% tolower(pos.lit[,1]))
        sub.de$lit[idx.lit] <- "lit"
        
        sign <- which(sub.de$significance=="padj < 0.01")
        sub.de$hit[sign] <- "pctr < 0.01"
        notsign <- which(sub.de$significance=="padj >= 0.01")
        sub.de$hit[notsign] <- "pctr >= 0.01"
        
        sub.de <- createSignLabelsNumbers(sub.de, "hit")
        pp <- pp + 
            geom_point(data=sub.de, aes(x=log2FoldChange, y=minuslog10pval, 
                                color=hit, shape=lit, 
                                text=paste0("padj=", padj, "\nname=", gene)), 
                                size=0.9)
    }
   
    
    
    pp <- pp + labs(list(title=title, x=xlabl, y=ylabl))
    return(pp)
}










