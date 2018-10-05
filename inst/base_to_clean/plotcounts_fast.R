

ens.gene <- "ENSRNOG00000019689"

gene.name <- as.character(norm.converted[which(rownames(norm.converted)==ens.gene),"symbol"])

prefix <- paste0(gene.name, "_gene")

PlotCountsAlongTimes(normalized.counts = norm.converted, design.matrix = cervical.design.matrix, gene.name = gene.name, 
                     gene.name.column.name = "symbol", show.plot.flag = TRUE, plotly.flag = TRUE, prefix.plot = prefix, plot.folder = "./", save.plot = TRUE)


PlotCountsAlongTimes(normalized.counts = normalized.counts.prop.uqua.names, design.matrix = cervical.design.matrix, gene.name = "stat3", 
                     gene.name.column.name = "symbol", show.plot.flag = TRUE, plotly.flag = TRUE)#, prefix.plot = prefix, plot.folder = "./", save.plot = FALSE)

res.cdtr.condtr <- results(de.results, name="Conditions_treated_vs_untreated", test="Wald")
subres.cdtr.condtr<-subset(res.cdtr.condtr, padj<0.05)
subres.cdtr.condtr <- subres.cdtr.condtr[order(subres.cdtr.condtr$padj),]


res.07.condtr <- results(de.results, name="Times07d.Conditionstreated", test="Wald")

dim(subres.07.condtr)

subres.07.condtr<-subset(res.07.condtr, padj<0.05)
subres.07.condtr <- subres.07.condtr[order(subres.07.condtr$padj),]



res.14.condtr <- results(de.results, name="Times14d.Conditionstreated", test="Wald")

dim(subres.14.condtr)

subres.14.condtr<-subset(res.14.condtr, padj<0.05)
subres.14.condtr <- subres.14.condtr[order(subres.14.condtr$padj),]


res.56.condtr <- results(de.results, name="Times56d.Conditionstreated", test="Wald")
subres.56.condtr<-subset(res.56.condtr, padj<0.05)
subres.56.condtr <- subres.56.condtr[order(subres.56.condtr$padj),]
dim(subres.56.condtr)


x <- rownames(subres.cdtr.condtr)
y<- rownames(subres.14.condtr)
z<- rownames(subres.56.condtr)

Venn3de(x,y,z, label1 = "wald", label2="07d", label3 = "03d", title = "venn", intersection.flag = TRUE, intersection.exclusion.flag = TRUE, plot.dir = "./")


x<-rownames(cerv.de.uqua.notnorm2.des4time.mid.formula$SIGN)
y<- rownames(subres.cdtr.condtr)
z<- rownames(cerv.de.uqua.notnorm2.des4time.first.formula$SIGN)


x<- rownames(subres.cdtr.condtr)
y<- rownames(cervical.de.uqua.notnorm2.07d$SIGN)
z<- rownames(cervical.de.uqua.notnorm2.03d.deseq$SIGN)

dim(subres.cdtr.condtr)
dim(cervical.de.uqua.notnorm2.07d$SIGN)
dim(cervical.de.uqua.notnorm2.14d$SIGN)


y<-rownames(cerv.de.uqua.notnorm2.des4time.final.formula2$SIGN)


PlotCountsAlongTimes(normalized.counts = norm.converted, design.matrix = cerv.thor.design.matrix, gene.name = "Serpina3n", 
                     gene.name.column.name = "symbol", show.plot.flag = TRUE, plotly.flag = TRUE, prefix.plot = "serpina_down_cerv_thor", plot.folder = "./", save.plot = TRUE)






#################################
d <- make.design.matrix(NBdesign)
NBp <- p.vector(NBdata, d, counts=TRUE)
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all")
get$summary

d
cc.dm

ccc.dm <- cervical.design.matrix
ccc.dm$Treated <- 0
ccc.dm$Untreated <- 0
ccc.dm$Replicates <- 0
conditions <- unique(ccc.dm$Conditions)

ccc.dm$Treated[which(ccc.dm$Conditions=="treated")] <- 1 
ccc.dm$Untreated[which(ccc.dm$Conditions=="untreated")] <- 1 

ccc.dm$Replicates <- as.numeric(as.factor(paste(ccc.dm$Times, ccc.dm$Conditions, sep = "_")))

ccc.dm$Times <- as.numeric(gsub("d", "", ccc.dm$Times))
ccc.dm.m <- as.matrix(ccc.dm[, c(1,4,3,5)])

ccc.dsm.m<-make.design.matrix(ccc.dm.m)

cervical.time.norm.data <- normalized.counts.prop.uqua[, which(colnames(normalized.counts.prop.uqua) %in% rownames(cervical.design.matrix))]

cerv.p.masig <- p.vector(cervical.time.norm.data, ccc.dsm.m, counts=TRUE, theta = 10)

cerv.p.masig.t <- T.fit(cerv.p.masig)

get<-get.siggenes(cerv.p.masig.t, vars="all")
length(get$summary)

dim(cerv.p.masig.t$influ.info)




# 
# see.genes(get$sig.genes, k = 1)
# 
# masig.res <- cbind(cerv.p.masig$p.vector, cerv.p.masig$p.adjusted)
# 
# colnames(masig.res) <- c("p.value", "p.adjusted")
# 
# masig.res <- as.data.frame(masig.res)
# dim(masig.res)
# 
# sign.masig.res <- subset(masig.res, p.adjusted < 0.05)
# 
# dim(sign.masig.res)
# 













