
### detect gene symbols in map and in pathways 

ens.specific.genes <- ens.symb.biomart.map[ which(tolower(ens.symb.biomart.map$Associated.Gene.Name) %in% c("serpina3n", "gfap", "lcn2")),]


synaptic.vesicle.cycle.kegg <- c("ENSRNOG00000000060", "ENSRNOG00000000840", "ENSRNOG00000003905", "ENSRNOG00000004560", "ENSRNOG00000006426", "ENSRNOG00000006542",
                                 "ENSRNOG00000011000", "ENSRNOG00000015420", "ENSRNOG00000015865", "ENSRNOG00000017220", "ENSRNOG00000019193", "ENSRNOG00000019317", 
                                 "ENSRNOG00000019433", "ENSRNOG00000021013", "ENSRNOG00000026490", "ENSRNOG00000030862", "ENSRNOG00000033835", "ENSRNOG00000036814")

sub.norm.counts.synaptic.vesicle.cycle <- normalized.counts.prop.uqua[which(rownames(normalized.counts.prop.uqua) %in% synaptic.vesicle.cycle.kegg), ]
dim(sub.norm.counts.synaptic.vesicle.cycle)
length(synaptic.vesicle.cycle.kegg)

syn.ves.cyc.l.fc.fc <- CalculateLogFoldChangeOfFoldChangesForEachTime(subset.counts = sub.norm.counts.synaptic.vesicle.cycle, design.matrix = all.des.mat, log.base = 2)

syn.ves.cyc.l.fc.fc.symb <- AttachConvertedColumn(de.data.frame = as.data.frame(syn.ves.cyc.l.fc.fc), conversion.map = ens.symb.biomart.map)

rownames(syn.ves.cyc.l.fc.fc.symb) <- syn.ves.cyc.l.fc.fc.symb$symbol
syn.ves.cyc.l.fc.fc.symb <- syn.ves.cyc.l.fc.fc.symb[, -5]

# which(is.infinite(syn.ves.cyc.l.fc.fc))
PlotHeatmap(as.matrix(syn.ves.cyc.l.fc.fc.symb), scale = "row", row.labels.percent = 1.5, dendrogram = "row", title = "4Points LUGs KEGG Synaptic Vescicle Cycle")

synaptic.vesicle.cycle.kegg.03d <- c( "ENSRNOG00000000060", "ENSRNOG00000000840", "ENSRNOG00000003905", "ENSRNOG00000004560", "ENSRNOG00000011000", "ENSRNOG00000015420", "ENSRNOG00000016147", "ENSRNOG00000019193", "ENSRNOG00000019317", "ENSRNOG00000019433", "ENSRNOG00000026490", "ENSRNOG00000030862", "ENSRNOG00000033835", "ENSRNOG00000036814")
length(synaptic.vesicle.cycle.kegg.03d)

sub.norm.counts.synaptic.vesicle.cycle.03d <- normalized.counts.prop.uqua[which(rownames(normalized.counts.prop.uqua) %in% synaptic.vesicle.cycle.kegg.03d), ]
dim(sub.norm.counts.synaptic.vesicle.cycle.03d)

syn.ves.cyc.l.fc.fc.03d <- CalculateLogFoldChangeOfFoldChangesForEachTime(subset.counts = sub.norm.counts.synaptic.vesicle.cycle.03d, design.matrix = all.des.mat, log.base = 2)

syn.ves.cyc.l.fc.fc.symb.03d <- AttachConvertedColumn(de.data.frame = as.data.frame(syn.ves.cyc.l.fc.fc.03d), conversion.map = ens.symb.biomart.map)

rownames(syn.ves.cyc.l.fc.fc.symb.03d) <- syn.ves.cyc.l.fc.fc.symb.03d$symbol
syn.ves.cyc.l.fc.fc.symb.03d <- syn.ves.cyc.l.fc.fc.symb.03d[, -5]

# which(is.infinite(syn.ves.cyc.l.fc.fc))
PlotHeatmap(as.matrix(syn.ves.cyc.l.fc.fc.symb.03d), scale = "row", row.labels.percent = 1.5, dendrogram = "row", title = "03d LUGs KEGG Synaptic Vescicle Cycle")


which(!rownames(syn.ves.cyc.l.fc.fc.symb.03d) %in% rownames(syn.ves.cyc.l.fc.fc.symb))


enrichment.filename <- "/media/dario/dati/time_course/cervical_thoracic_whole/Thu_Jun_15_11_18_26_2017/4Points/4Points_Itersections/SIGN/lrt.wald.filtered.union.padj005/Venn3/venn3_union_of_LUGs_Cervical_Union_Thoracic_Union_Cerv-Thor_lrt_wald_filtered_union_padj005/functional/gProfiler/KEGG/venn3_union_of_LUGs_Cervical_Union_Thoracic_Union_Cerv-Thor_lrt_wald_filtered_union_padj005_genes_KEGG__PATHWAYS_significant.tsv"
enrichment.file <- ReadDataFrameFromTsv(enrichment.filename, row.names.col = 1, header.flag = TRUE)

ens.genes <- strsplit(x = as.character(enrichment.file$intersection), split = "; ")
names(ens.genes) <- enrichment.file$term.name


# plotInterestingHeatmap(counts = normalized.counts.prop.uqua, interesting.gene.list = as.character(unlist(interesting.term.genes)), design.matrix=all.des.mat, gene.map=ens.symb.biomart.map, heat.title=paste0("4Points LUGs ", interesting.term.name))
# interesting.term.name = "GABAergic synapse"; ens.genes = ens.genes; counts = normalized.counts.prop.uqua; design.matrix = all.des.mat; gene.map = ens.symb.biomart.map
identifyTermAndPlotHeatmap(interesting.term.name = "GABAergic synapse", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title)

identifyTermAndPlotHeatmap <- function(interesting.term.name, ens.genes, counts, design.matrix, gene.map, log.base=2, prefix.title) {
  # interesting.term.name = "GABAergic synapse"
  
  interesting.term.genes <- ens.genes[which(names(ens.genes) == interesting.term.name)]
  
  plotInterestingHeatmap(counts = counts, interesting.gene.list = as.character(unlist(interesting.term.genes)), design.matrix=design.matrix, gene.map=gene.map, heat.title=paste(prefix.title, interesting.term.name, sep = " "))
}

plotInterestingHeatmap <- function(counts, interesting.gene.list, design.matrix, gene.map, heat.title, log.base=2) {
  
  interesting.counts <- normalized.counts.prop.uqua[which(rownames(counts) %in% interesting.gene.list), ]
  # dim(sub.norm.counts.synaptic.vesicle.cycle.03d)
  
  interesting.l.fc.fc <- CalculateLogFoldChangeOfFoldChangesForEachTime(subset.counts = interesting.counts, design.matrix = design.matrix, log.base = log.base)
  
  interesting.l.fc.fc.symb <- AttachConvertedColumn(de.data.frame = as.data.frame(interesting.l.fc.fc), conversion.map = gene.map)
  
  if( dim(which(is.infinite(interesting.l.fc.fc), arr.ind = TRUE))[1] > 0 ) {
    message("infinite present in lfcfc")
    print(interesting.l.fc.fc.symb)
    message("removing infinite row")
    print(which(is.infinite(interesting.l.fc.fc), arr.ind = TRUE))
    # interesting.l.fc.fc.symb <- interesting.l.fc.fc.symb[-which(is.infinite(interesting.l.fc.fc), arr.ind = TRUE)[1,1], ] 
    print(interesting.l.fc.fc.symb)
    }
  
  rownames(interesting.l.fc.fc.symb) <- interesting.l.fc.fc.symb$symbol
  interesting.l.fc.fc.symb <- interesting.l.fc.fc.symb[, -5]
  
  # which(is.infinite(syn.ves.cyc.l.fc.fc))
  filename <- file.path(getwd(), gsub(pattern = " ", replacement = "_", x = heat.title))
  PlotHeatmap(as.matrix(interesting.l.fc.fc.symb), scale = "row", row.labels.percent = 0.45, dendrogram = "row", title = heat.title, save.png = TRUE, file.name = filename)
}


# identifyTermAndPlotHeatmap(interesting.term.name = "Cardiac muscle contraction", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map)

# identifyTermAndPlotHeatmap(interesting.term.name = "Cardiac muscle contraction", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map)


# enrichment.path <- "/media/dario/dati/time_course/cervical_thoracic_whole/Thu_Jun_15_11_18_26_2017/Time_by_Time/03d/Intersections_mix/UP/dataset.tr.minus.untr.padj005/Venn3/venn3_union_of_LUGs_Cervical_Thoracic_Cerv-Thor_UP_dataset_tr_minus_untr_padj005/functional/gProfiler/GO/"
enrichment.filename <- "/media/dario/dati/time_course/cervical_thoracic_whole/Thu_Jun_15_11_18_26_2017/Time_by_Time/03d/Intersections/SIGN/dataset.tr.minus.untr.padj005/Venn3/venn3_union_of_LUGs_Cervical_SIGN_Thoracic_SIGN_Cerv-Thor_SIGN_dataset_tr_minus_untr_padj005/functional/gProfiler/KEGG/venn3_union_of_LUGs_Cervical_SIGN_Thoracic_SIGN_Cerv-Thor_SIGN_dataset_tr_minus_untr_padj005_genes_KEGG__PATHWAYS_significant.tsv"
enrichment.file.3d <- ReadDataFrameFromTsv(enrichment.filename, row.names.col = 1, header.flag = TRUE)

ens.genes <- strsplit(x = as.character(enrichment.file.3d$intersection), split = "; ")
names(ens.genes) <- enrichment.file.3d$term.name

identifyTermAndPlotHeatmap(interesting.term.name = "GABAergic synapse", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="03d LUGs")
identifyTermAndPlotHeatmap(interesting.term.name = "Cardiac muscle contraction", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="03d LUGs")
identifyTermAndPlotHeatmap(interesting.term.name = "Adrenergic signaling in cardiomyocytes", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="03d LUGs")
identifyTermAndPlotHeatmap(interesting.term.name = "Synaptic vesicle cycle", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="03d LUGs")

setwd("/media/dario/dati/time_course/cervical_thoracic_whole/Thu_Jun_15_11_18_26_2017/heatmaps/7d")
enrichment.filename.7d <- "/media/dario/dati/time_course/cervical_thoracic_whole/Thu_Jun_15_11_18_26_2017/Time_by_Time/07d/Intersections/SIGN/dataset.tr.minus.untr.padj005/Venn3/venn3_union_of_LUGs_Cervical_SIGN_Thoracic_SIGN_Cerv-Thor_SIGN_dataset_tr_minus_untr_padj005/functional/gProfiler/KEGG/venn3_union_of_LUGs_Cervical_SIGN_Thoracic_SIGN_Cerv-Thor_SIGN_dataset_tr_minus_untr_padj005_genes_KEGG__PATHWAYS_significant.tsv"
enrichment.file.7d <- ReadDataFrameFromTsv(enrichment.filename.7d, row.names.col = 1, header.flag = TRUE)

ens.genes <- strsplit(x = as.character(enrichment.file.7d$intersection), split = "; ")
names(ens.genes) <- enrichment.file.7d$term.name
identifyTermAndPlotHeatmap(interesting.term.name = "GABAergic synapse", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="07d LUGs")
identifyTermAndPlotHeatmap(interesting.term.name = "Adrenergic signaling in cardiomyocytes", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="07d LUGs")
identifyTermAndPlotHeatmap(interesting.term.name = "Synaptic vesicle cycle", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="07d LUGs")
identifyTermAndPlotHeatmap(interesting.term.name = "ECM-receptor interaction", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="07d LUGs")


plotInterestingHeatmap <- function(counts, interesting.gene.list, design.matrix, gene.map, heat.title, log.base=2) {
  
  interesting.counts <- normalized.counts.prop.uqua[which(rownames(counts) %in% interesting.gene.list), ]
  # dim(sub.norm.counts.synaptic.vesicle.cycle.03d)
  
  interesting.l.fc.fc <- CalculateLogFoldChangeOfFoldChangesForEachTime(subset.counts = interesting.counts, design.matrix = design.matrix, log.base = log.base)
  
  
  interesting.l.fc.fc.symb <- AttachConvertedColumn(de.data.frame = as.data.frame(interesting.l.fc.fc), conversion.map = gene.map)
  
  if( dim(which(is.infinite(interesting.l.fc.fc), arr.ind = TRUE))[1] > 0 ) {
    message("infinite present in lfcfc")
    print(interesting.l.fc.fc.symb)
    message("removing infinite row")
    ar.ind <- which(is.infinite(interesting.l.fc.fc), arr.ind = TRUE)
    # print(ar.ind)
    # print(ar.ind[,1])
    for(row.ind in ar.ind[,1]) {
      interesting.l.fc.fc <- interesting.l.fc.fc[-row.ind, ] 
    }
    interesting.l.fc.fc.symb <- AttachConvertedColumn(de.data.frame = as.data.frame(interesting.l.fc.fc), conversion.map = gene.map)
    print(interesting.l.fc.fc.symb)
    
  }
  
  rownames(interesting.l.fc.fc.symb) <- interesting.l.fc.fc.symb$symbol
  interesting.l.fc.fc.symb <- interesting.l.fc.fc.symb[, -5]
  
  # which(is.infinite(syn.ves.cyc.l.fc.fc))
  filename <- file.path(getwd(), gsub(pattern = " ", replacement = "_", x = heat.title))
  PlotHeatmap(as.matrix(interesting.l.fc.fc.symb), scale = "row", row.labels.percent = 1.5, dendrogram = "row", title = heat.title, save.png = TRUE, file.name = filename)
}


setwd("/media/dario/dati/time_course/cervical_thoracic_whole/Thu_Jun_15_11_18_26_2017/heatmaps/14d")
enrichment.filename.14d <- "/media/dario/dati/time_course/cervical_thoracic_whole/Thu_Jun_15_11_18_26_2017/Time_by_Time/14d/Intersections/SIGN/dataset.tr.minus.untr.padj005/Venn3/venn3_union_of_LUGs_Cervical_SIGN_Thoracic_SIGN_Cerv-Thor_SIGN_dataset_tr_minus_untr_padj005/functional/gProfiler/KEGG/venn3_union_of_LUGs_Cervical_SIGN_Thoracic_SIGN_Cerv-Thor_SIGN_dataset_tr_minus_untr_padj005_genes_KEGG__PATHWAYS_significant.tsv"
enrichment.file.14d <- ReadDataFrameFromTsv(enrichment.filename.14d, row.names.col = 1, header.flag = TRUE)

ens.genes <- strsplit(x = as.character(enrichment.file.14d$intersection), split = "; ")
names(ens.genes) <- enrichment.file.14d$term.name
identifyTermAndPlotHeatmap(interesting.term.name = "Phagosome", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="14d LUGs")

identifyTermAndPlotHeatmap(interesting.term.name = "Adrenergic signaling in cardiomyocytes", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="14d LUGs")

identifyTermAndPlotHeatmap(interesting.term.name = "Synaptic vesicle cycle", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="14d LUGs")

identifyTermAndPlotHeatmap(interesting.term.name = "Complement and coagulation cascades", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="14d LUGs")


setwd("/media/dario/dati/time_course/cervical_thoracic_whole/Thu_Jun_15_11_18_26_2017/heatmaps/56d")
enrichment.filename.56d <- "/media/dario/dati/time_course/cervical_thoracic_whole/Thu_Jun_15_11_18_26_2017/Time_by_Time/56d/Intersections/SIGN/dataset.tr.minus.untr.padj005/Venn3/venn3_union_of_LUGs_Cervical_SIGN_Thoracic_SIGN_Cerv-Thor_SIGN_dataset_tr_minus_untr_padj005/functional/gProfiler/KEGG/venn3_union_of_LUGs_Cervical_SIGN_Thoracic_SIGN_Cerv-Thor_SIGN_dataset_tr_minus_untr_padj005_genes_KEGG__PATHWAYS_significant.tsv"
enrichment.file.56d <- ReadDataFrameFromTsv(enrichment.filename.56d, row.names.col = 1, header.flag = TRUE)

ens.genes <- strsplit(x = as.character(enrichment.file.56d$intersection), split = "; ")
names(ens.genes) <- enrichment.file.56d$term.name
identifyTermAndPlotHeatmap(interesting.term.name = "ECM-receptor interaction", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="56d LUGs")
identifyTermAndPlotHeatmap(interesting.term.name = "Protein digestion and absorption", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="56d LUGs")

ens.symb.biomart.map[which(tolower(ens.symb.biomart.map$Associated.Gene.Name) %in% tolower(pos.reg.necr$V1)),]
dim(pos.reg.necr)

apoptotic.mapped.ens <- ens.symb.biomart.map[which(tolower(ens.symb.biomart.map$Associated.Gene.Name) %in% tolower(unique(apoptotic$V1))),"Ensembl.Gene.ID"]
identifyTermAndPlotHeatmap(interesting.term.name = "ECM-receptor interaction", ens.genes = ens.genes, counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, gene.map = ens.symb.biomart.map, prefix.title="genes")
