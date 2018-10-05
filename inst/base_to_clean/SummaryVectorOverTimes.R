
cerv.03.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/03d/Cervical_deseq/DE_results/DeSeq/Cervical_03d_prop_uqua_DeSeq_all.tsv")
cerv.07.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/07d/Cervical_deseq/DE_results/DeSeq/Cervical_07d_prop_uqua_DeSeq_all.tsv")
cerv.14.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/14d/Cervical_deseq/DE_results/DeSeq/Cervical_14d_prop_uqua_DeSeq_all.tsv")
cerv.56.deseq <- ReadDataFrameFromTsv("/media/dario/dati/time_course/cervical_thoracic_whole/results_final_20/Time_by_Time/56d/Cervical_deseq/DE_results/DeSeq/Cervical_56d_prop_uqua_DeSeq_all.tsv")

head(cerv.03.deseq)

times <- c("03d", "07d", "14d", "56d")

cervical.deseq.list <- list(cerv.03.deseq, cerv.07.deseq, cerv.14.deseq, cerv.56.deseq)
names(cervical.deseq.list) <- times

genes <- rownames(cerv.03.deseq)

summary.df <- data.frame(row.names = genes)
i=1
for (deseq.df in cervical.deseq.list) {
  time <- times[i]
  dde.genes <- SignificantDeGenesPAdj(deseq.df)
  summary.df[,time] <- 0
  summary.df[which(rownames(summary.df)%in%rownames(dde.genes$UP)), time] <- 1
  summary.df[which(rownames(summary.df)%in%rownames(dde.genes$DOWN)), time] <- -1
  i=i+1
}

cerv.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN
summary.df[, "LRT"] <- 0
summary.df[which(rownames(summary.df) %in% rownames(cerv.de.uqua.notnorm2.des4time.lrt.wald$LRT$SIGN)), "LRT"] <- 1

summary.df[, "Wald"] <- 0
summary.df[which(rownames(summary.df) %in% rownames(cerv.de.uqua.notnorm2.des4time.lrt.wald$waldLRT$SIGN)), "Wald"] <- 1

summary.df$both <- apply(summary.df, 1, function(x) { sum(x["Wald"], x["LRT"]) } )
summary.df$both[summary.df$both < 2] <- 0


summary.df$test <- apply(summary.df, 1, function(x) { sum(x["both"], x["LRT"]) } )
ind.lrt <- which(summary.df$test==1)
summary.df$test[ind.lrt] <- 0

summary.df$test <- apply(summary.df, 1, function(x) { sum(x["test"], x["Wald"]) } )
ind.wald <- which(summary.df$test==1)

# which(summary.df$tmp==4)
summary.df[summary.df==0] <- "NO_CHANGE"
summary.df$test[ind.lrt] <- "LRT"
summary.df$test[ind.wald] <- "WALD"
summary.df$test[summary.df$test == 4] <- "LRT+WALD"
summary.df$`03d`[summary.df$`03d`==1] <- "UP"
summary.df$`03d`[summary.df$`03d`==-1] <- "DOWN"

summary.df$`07d`[summary.df$`07d`==1] <- "UP"
summary.df$`07d`[summary.df$`07d`==-1] <- "DOWN"
summary.df$`14d`[summary.df$`14d`==1] <- "UP"
summary.df$`14d`[summary.df$`14d`==-1] <- "DOWN"
summary.df$`56d`[summary.df$`56d`==1] <- "UP"
summary.df$`56d`[summary.df$`56d`==-1] <- "DOWN"

# summary.df$LRT[summary.df$LRT == 1] <- "LRT"
# summary.df$Wald[summary.df$Wald == 1] <- "Wald"

summ.df.genes <- AttachConvertedColumn(de.data.frame = summary.df, conversion.map = ens.symb.biomart.map)

summary.converted<-AttachConvertedColumn(de.data.frame = summary.df, conversion.map = ens.symb.biomart.map)
WriteDataFrameAsTsv(file.name.path = "cervical_thoracic_whole/summary_results/summary_df_31_03_2017",data.frame.to.save = summary.df)
WriteDataFrameAsTsv(file.name.path = "cervical_thoracic_whole/summary_results/summary_df_genesymbols_31_03_2017",data.frame.to.save = summary.converted)
# heatmap.2(x = as.matrix(summary.df))




tmp.test <- summ.df.genes[, c("03d", "07d", "14d", "56d", "test")]
ddf <- as.data.frame(rep(summ.df.genes$symbol, times=dim(tmp.test)[2]))
colnames(ddf) <- "genes"
ddf$times <- sort(rep(colnames(tmp.test), each=dim(tmp.test)[1]))
ddf$values <- as.factor(c(tmp.test$`03d`, tmp.test$`07d`, tmp.test$`14d`, tmp.test$`56d`, tmp.test$test))
cervical.summary.ddf <- ddf
save("cervical.summary.ddf", file="heatmap_files/report_files/cervical_summary_df.RData")
# plotSummary(
# p <- ggplot(ddf, aes(x=times, y=genes, fill=values)) + geom_dotplot(binaxis='y', stackdir='center',  dotsize=0.2) + scale_fill_manual(values = c("DOWN"="red", "0"="white", "UP"="green", "LRT"="turquoise", "Wald"="plum1", "both"="orange"))

# pdf("prova.pdf", width = 20, height = 20)
# print(p)
# dev.off()

# 
# sd.test <- summary.df[1:100,]
# 
# new.df <- as.data.frame(rep(rownames(sd.test), times=dim(sd.test)[2]))
# colnames(new.df) <- "genes"
# 
# sd.test <- sd.test[, c("03d", "07d", "14d", "56d", "tmp")]

# new.df$times <- sort(rep(colnames(sd.test), 100))

# new.df$values <- as.factor(c(sd.test$`03d`, sd.test$`07d`, sd.test$`14d`, sd.test$`56d`, sd.test$tmp))


# dotchart(t(as.matrix(summary.df[1:100,])),  cex=.7, main="Gas Milage for Car Models", xlab="Miles Per Gallon")

plotSummary <- function(summary.df) {
  p <- ggplot(summary.df, aes(x=times, y=genes, fill=values)) + geom_dotplot(binaxis='y', stackdir='center',  dotsize=0.25) + scale_fill_manual(values = c("DOWN"="red", "0"="white", "UP"="green", "LRT"="turquoise", "Wald"="plum1", "both"="orange"))
  plot(p)
}


# p<-ggplot(new.df, aes(x=times, y=genes, fill=values)) + geom_dotplot(binaxis='y', stackdir='center',  dotsize=0.25) + scale_fill_manual(name = "Genes in Times", values = c("DOWN"="red", "0"="white", "UP"="green")) + scale_colour_manual(name = "Tests", values = c("LRT"="turquoise", "Wald"="plum1", "both"="orange"))

p#+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# ggplotly(p)












