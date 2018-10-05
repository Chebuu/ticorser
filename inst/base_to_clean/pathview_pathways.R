oxphosgenelist <- c("ENSRNOG00000000557", "ENSRNOG00000000840", "ENSRNOG00000001170", "ENSRNOG00000006542", "ENSRNOG00000006939", "ENSRNOG00000007235", "ENSRNOG00000012550", 
                    "ENSRNOG00000014625", "ENSRNOG00000015320", "ENSRNOG00000016660", "ENSRNOG00000016952", "ENSRNOG00000017220", "ENSRNOG00000017446", "ENSRNOG00000020602", "ENSRNOG00000021177", 
                    "ENSRNOG00000024539", "ENSRNOG00000026616", "ENSRNOG00000030862", "ENSRNOG00000036814", "ENSRNOG00000040005", "ENSRNOG00000048174", "ENSRNOG00000049912")

interested.genes.kegg.lfcfc <- all.genes.symb[which(rownames(all.genes.symb) %in% oxphosgenelist),]

rownames(interested.genes.kegg.lfcfc) <- interested.genes.kegg.lfcfc$symbol

interested.genes.kegg.lfcfc <- interested.genes.kegg.lfcfc[,-5]

PlotKeggMapTimeCoursePathview(interested.genes.kegg.lfc = interested.genes.kegg.lfcfc, kegg.id = "00190", specie = "rno", gene.idtype = "SYMBOL", output.dir = "./", prefix = " ")


gabaergic.cerv.thor.tc.lugs <- c("ENSRNOG00000000007; ENSRNOG00000002229; ENSRNOG00000002349; ENSRNOG00000002360; ENSRNOG00000002636; ENSRNOG00000003241; ENSRNOG00000003491; ENSRNOG00000003905; ENSRNOG00000004560; ENSRNOG00000005291; ENSRNOG00000005697; ENSRNOG00000006108; ENSRNOG00000008431; ENSRNOG00000012876; ENSRNOG00000017417; ENSRNOG00000018111; ENSRNOG00000018200; ENSRNOG00000019482")
gabaergic.cerv.thor.tc.lugs <- strsplit(x = gabaergic.cerv.thor.tc.lugs, split = "; ", fixed = TRUE)

interested.genes.kegg.lfcfc <- all.genes.symb[which(rownames(all.genes.symb) %in% gabaergic.cerv.thor.tc.lugs[[1]]),]
dim(interested.genes.kegg.lfcfc)
rownames(interested.genes.kegg.lfcfc) <- interested.genes.kegg.lfcfc$symbol
interested.genes.kegg.lfcfc <- interested.genes.kegg.lfcfc[,-5]
require("pathview")
PlotKeggMapTimeCoursePathview(interested.genes.kegg.lfc = interested.genes.kegg.lfcfc, kegg.id = "04727", specie = "rno", gene.idtype = "SYMBOL", output.dir = "./", prefix = " ")

oxphos.03d.lugs <- c("ENSRNOG00000001170; ENSRNOG00000006939; ENSRNOG00000007235; ENSRNOG00000012550; ENSRNOG00000014625; ENSRNOG00000016660; ENSRNOG00000016952; ENSRNOG00000017446; ENSRNOG00000020602; ENSRNOG00000021177; ENSRNOG00000024568; ENSRNOG00000026616; ENSRNOG00000028717; ENSRNOG00000040005; ENSRNOG00000048174; ENSRNOG00000048320; ENSRNOG00000049912")
oxphos.03d.lugs <- strsplit(x = oxphos.03d.lugs, split = "; ", fixed = TRUE)
interested.genes.kegg.lfcfc <- all.genes.symb[which(rownames(all.genes.symb) %in% oxphos.03d.lugs[[1]]),]
dim(interested.genes.kegg.lfcfc)
rownames(interested.genes.kegg.lfcfc) <- interested.genes.kegg.lfcfc$symbol
interested.genes.kegg.lfcfc <- interested.genes.kegg.lfcfc[,-5]
require("pathview")
PlotKeggMapTimeCoursePathview(interested.genes.kegg.lfc = interested.genes.kegg.lfcfc, kegg.id = "00190", specie = "rno", gene.idtype = "SYMBOL", output.dir = "./", prefix = " ")


parkinson.03d.lugs <- c("ENSRNOG00000001170; ENSRNOG00000006939; ENSRNOG00000007235; ENSRNOG00000012550; ENSRNOG00000014625; ENSRNOG00000016660; ENSRNOG00000016952; ENSRNOG00000017446; ENSRNOG00000020602; ENSRNOG00000021177; ENSRNOG00000024568; ENSRNOG00000026616; ENSRNOG00000028717; ENSRNOG00000040005; ENSRNOG00000048174; ENSRNOG00000048320; ENSRNOG00000049912")
parkinson.03d.lugs <- strsplit(x = parkinson.03d.lugs, split = "; ", fixed = TRUE)
interested.genes.kegg.lfcfc <- all.genes.symb[which(rownames(all.genes.symb) %in% parkinson.03d.lugs[[1]]),]
dim(interested.genes.kegg.lfcfc)
rownames(interested.genes.kegg.lfcfc) <- interested.genes.kegg.lfcfc$symbol
interested.genes.kegg.lfcfc <- interested.genes.kegg.lfcfc[,-5]
require("pathview")
PlotKeggMapTimeCoursePathview(interested.genes.kegg.lfc = interested.genes.kegg.lfcfc, kegg.id = "05012", specie = "rno", gene.idtype = "SYMBOL", output.dir = "./", prefix = " ")



oxphos.tc.lugs <- c("ENSRNOG00000001170; ENSRNOG00000006939; ENSRNOG00000007235; ENSRNOG00000014625; ENSRNOG00000016660; ENSRNOG00000016952; ENSRNOG00000017446; ENSRNOG00000020602; ENSRNOG00000021177; ENSRNOG00000024539; ENSRNOG00000026616; ENSRNOG00000040005; ENSRNOG00000048174; ENSRNOG00000049912")
oxphos.tc.lugs <- strsplit(x = oxphos.tc.lugs, split = "; ", fixed = TRUE)
interested.genes.kegg.lfcfc <- all.genes.symb[which(rownames(all.genes.symb) %in% oxphos.tc.lugs[[1]]),]
length(oxphos.tc.lugs[[1]])
dim(interested.genes.kegg.lfcfc)
rownames(interested.genes.kegg.lfcfc) <- interested.genes.kegg.lfcfc$symbol
interested.genes.kegg.lfcfc <- interested.genes.kegg.lfcfc[,-5]
require("pathview")
PlotKeggMapTimeCoursePathview(interested.genes.kegg.lfc = interested.genes.kegg.lfcfc, kegg.id = "00190", specie = "rno", gene.idtype = "SYMBOL", output.dir = "./", prefix = " ")


parkinson.tc.lugs <- c("ENSRNOG00000001170; ENSRNOG00000006939; ENSRNOG00000007235; ENSRNOG00000012181; ENSRNOG00000014625; ENSRNOG00000016660; ENSRNOG00000016952; ENSRNOG00000017446; ENSRNOG00000020602; ENSRNOG00000021177; ENSRNOG00000024539; ENSRNOG00000026616; ENSRNOG00000038202; ENSRNOG00000040005; ENSRNOG00000048174; ENSRNOG00000049912")
parkinson.tc.lugs <- strsplit(x = parkinson.tc.lugs, split = "; ", fixed = TRUE)
interested.genes.kegg.lfcfc <- all.genes.symb[which(rownames(all.genes.symb) %in% parkinson.tc.lugs[[1]]),]
length(parkinson.tc.lugs[[1]])
dim(interested.genes.kegg.lfcfc)
rownames(interested.genes.kegg.lfcfc) <- interested.genes.kegg.lfcfc$symbol
interested.genes.kegg.lfcfc <- interested.genes.kegg.lfcfc[,-5]
require("pathview")
PlotKeggMapTimeCoursePathview(interested.genes.kegg.lfc = interested.genes.kegg.lfcfc, kegg.id = "05012", specie = "rno", gene.idtype = "SYMBOL", output.dir = "./", prefix = " ")



alzheimer.03d.lugs <- c("ENSRNOG00000001170; ENSRNOG00000006939; ENSRNOG00000007235; ENSRNOG00000012550; ENSRNOG00000014625; ENSRNOG00000016660; ENSRNOG00000016952; ENSRNOG00000017446; ENSRNOG00000020602; ENSRNOG00000021177; ENSRNOG00000024568; ENSRNOG00000026616; ENSRNOG00000028717; ENSRNOG00000038202; ENSRNOG00000040005; ENSRNOG00000048174; ENSRNOG00000048320; ENSRNOG00000049912")
alzheimer.03d.lugs <- strsplit(x = alzheimer.03d.lugs, split = "; ", fixed = TRUE)
interested.genes.kegg.lfcfc <- all.genes.symb[which(rownames(all.genes.symb) %in% alzheimer.03d.lugs[[1]]),]
dim(interested.genes.kegg.lfcfc)
rownames(interested.genes.kegg.lfcfc) <- interested.genes.kegg.lfcfc$symbol
interested.genes.kegg.lfcfc <- interested.genes.kegg.lfcfc[,-5]
require("pathview")
PlotKeggMapTimeCoursePathview(interested.genes.kegg.lfc = interested.genes.kegg.lfcfc, kegg.id = "05010", specie = "rno", gene.idtype = "SYMBOL", output.dir = "./", prefix = " ")


alzheimer.tc.lugs <- c("ENSRNOG00000001170; ENSRNOG00000006939; ENSRNOG00000007235; ENSRNOG00000012181; ENSRNOG00000014625; ENSRNOG00000016660; ENSRNOG00000016952; ENSRNOG00000017446; ENSRNOG00000020602; ENSRNOG00000021177; ENSRNOG00000024539; ENSRNOG00000026616; ENSRNOG00000038202; ENSRNOG00000040005; ENSRNOG00000048174; ENSRNOG00000049912")
alzheimer.tc.lugs <- strsplit(x = alzheimer.tc.lugs, split = "; ", fixed = TRUE)
interested.genes.kegg.lfcfc <- all.genes.symb[which(rownames(all.genes.symb) %in% alzheimer.tc.lugs[[1]]),]
dim(interested.genes.kegg.lfcfc)
rownames(interested.genes.kegg.lfcfc) <- interested.genes.kegg.lfcfc$symbol
interested.genes.kegg.lfcfc <- interested.genes.kegg.lfcfc[,-5]
require("pathview")
PlotKeggMapTimeCoursePathview(interested.genes.kegg.lfc = interested.genes.kegg.lfcfc, kegg.id = "05010", specie = "rno", gene.idtype = "SYMBOL", output.dir = "./", prefix = " ")
