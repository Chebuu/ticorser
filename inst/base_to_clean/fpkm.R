
cerv.thor.cnts.t <- read.table("/media/dario/dati/time_course_4/cervical_thoracic/cervical_thoracic_results/results_t/counts/Cervical_Thoracic_Counts_Renamed_FeatureCounts.tsv")
cerv.thor.cnts.u <- read.table("/media/dario/dati/time_course_4/cervical_thoracic/cervical_thoracic_results/results_u/counts/FeatureCounts/Cervical_Thoracic_Untr_Counts_Renamed_FeatureCounts.tsv")
cerv.thor.des <- read.table("/media/dario/dati/time_course_4/cervical_thoracic/design_file/cervical_thoracic_all_design_file.txt")

cerv.thor.cnts <- cbind(cerv.thor.cnts.t, cerv.thor.cnts.u)
cerv.thor.cnts <- cerv.thor.cnts[,order(colnames(cerv.thor.cnts))]

library("DESeq2")

# se <- SummarizedExperiment(list(counts=cerv.thor.cnts), colData=cerv.thor.des)
# 
# dse <- DESeqDataSetFromMatrix(countData=as.matrix(cerv.thor.cnts), colData=(cerv.thor.des), ~1)



annotation.file <- "/media/dario/dati/time_course_4/downloaded_references/ensembl_site/Rattus_norvegicus.Rnor_5.0.79_gene_biotype_protein_coding_header.gtf"
all.bam.folder <- "/media/dario/dati/time_course_4/cervical_thoracic/all_bams/"

bam.files <- list.files(all.bam.folder, "bam$", full.names=TRUE)

fc_SE <- Rsubread::featureCounts(files = bam.files, 
                                 annot.ext = annotation.file, 
                                 isGTFAnnotationFile=TRUE, 
                                 GTF.featureType = 'gene', ## def
                                 GTF.attrType = "gene_id", ## gene_id def
                                 useMetaFeatures = TRUE, ## def
                                 allowMultiOverlap = FALSE, ## def
                                 nthreads = 8, 
                                 strandSpecific = 0, ## 0 def
                                 countMultiMappingReads = FALSE, ## def
                                 isPairedEnd = FALSE) ## FALSE def

save(fc_SE, file="fc_SE_for_fpkm.RData")
saveRDS(fc_SE, file="fc_SE_for_fpkm.Rds")

grl <- lapply(c(1:dim(fc_SE$annotation)[1]), function(i) 
{
    chr <- unique(strsplit(fc_SE$annotation$Chr[i], ";")[[1]])
    st <- strsplit(fc_SE$annotation$Start[i], ";")[[1]][1]
    len <- fc_SE$annotation$Length[i]
    # cat(chr, "-", st, "-", len, "\n")
    GRanges(seqnames=chr, ranges=IRanges(start=as.numeric(st), width=len))
})

grl.anno <- GRangesList(grl)


cnts <- fc_SE$counts[, order(colnames(fc_SE$counts))]

des <- cerv.thor.des[order(rownames(cerv.thor.des)),]

colnames(cnts) <- rownames(des)



se <- SummarizedExperiment(list(counts=cnts), colData=des, rowRanges=grl.anno)

dds <- DESeqDataSet(se, ~1)

fpkm <- fpkm(dds)

class(fpkm)
head(fpkm)[,1:5]
head(cnts)[,1:5]
write.table(x=fpkm, file=)


grl.anno <- GRanges(seqnames=fc_SE_gene$annotation$Chr, ranges=IRanges(start=as.numeric(fc_SE_gene$annotation$Start), 
                                                                       width=fc_SE_gene$annotation$Length))
cnts_gene <- fc_SE$counts[, order(colnames(fc_SE$counts))]
colnames(cnts_gene) <- rownames(des)
se_gene <- SummarizedExperiment(list(counts=cnts_gene), colData=des, rowRanges=grl.anno)
dds_gene <- DESeqDataSet(se_gene, ~1)

fpkm_gene <- fpkm(dds_gene)
head(fpkm_gene)[,1:5]
head(cnts_gene)[,1:5]
write.table(x=fpkm_gene, file="/media/dario/dati/iSync/works/IAC/coding/time_course/fpkm_on_gene.tsv",quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
