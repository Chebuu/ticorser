library(pathview)
data("gse16873.d")
# 
# sim.cpd.data=sim.mol.data(mol.type="cpd", nmol=3000)
# 
# set.seed(10)
# sim.cpd.data2 = matrix(sample(sim.cpd.data, 18000, replace = T), ncol = 6)
all.genes <- CalculateLogFoldChangeOfFoldChangesForEachTime(subset.counts = normalized.counts.prop.uqua, design.matrix = all.des.mat, log.base = 2)
all.genes.symb <- AttachConvertedColumn(de.data.frame = as.data.frame(all.genes), conversion.map = ens.symb.biomart.map)


synaptic.vesicle.cycle.kegg <- c("ENSRNOG00000000060", "ENSRNOG00000000840", "ENSRNOG00000003905", "ENSRNOG00000004560", "ENSRNOG00000006426", "ENSRNOG00000006542",
                                 "ENSRNOG00000011000", "ENSRNOG00000015420", "ENSRNOG00000015865", "ENSRNOG00000017220", "ENSRNOG00000019193", "ENSRNOG00000019317", 
                                 "ENSRNOG00000019433", "ENSRNOG00000021013", "ENSRNOG00000026490", "ENSRNOG00000030862", "ENSRNOG00000033835", "ENSRNOG00000036814")
interested.genes.kegg.lfcfc <- all.genes.symb[which(rownames(all.genes.symb) %in% synaptic.vesicle.cycle.kegg),]

rownames(interested.genes.kegg.lfcfc) <- interested.genes.kegg.lfcfc$symbol

interested.genes.kegg.lfcfc <- interested.genes.kegg.lfcfc[,-5]
# colnames(interested.genes.kegg.lfcfc)
pathview(gene.data = interested.genes.kegg.lfcfc, pathway.id = "04721", species = "rno", gene.idtype = "SYMBOL", 
         kegg.dir = "./", out.suffix="synaptic_vescicle_cycle_tc", #keys.align="y", 
         # is.signal=FALSE,
         low = c("gene"="red", cpd="red"),
         mid = c("gene"="grey", cpd="grey"),
         high = c("gene"="green", cpd="green"),
         # cols.ts.gene= as.vector(c("Red", "Grey", "Green")),
         kegg.native=TRUE, muti.state=TRUE, 
         same.layer=TRUE)


PlotKeggMapTimeCoursePathview <- function(interested.genes.kegg.lfc, kegg.id, specie, gene.idtype, output.dir, prefix, plot.col.key.flag=TRUE, is.signal.flag=TRUE, low.color.list = c("gene"="red", cpd="red"), mid.color.list= c("gene"="grey", cpd="grey"), high.color.list = c("gene"="green", cpd="green")) {

  
  pathview(gene.data = interested.genes.kegg.lfcfc, pathway.id = kegg.id, species = specie, 
           gene.idtype = gene.idtype, kegg.dir = output.dir, out.suffix=prefix, 
           plot.col.key=plot.col.key.flag,
           low = low.color.list,
           mid = mid.color.list,
           high = high.color.list,
           is.signal=is.signal.flag,
           # keys.align="y", 
           kegg.native=TRUE, muti.state=TRUE, same.layer=TRUE)  
}
