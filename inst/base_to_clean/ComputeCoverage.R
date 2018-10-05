## compute coverages

genome.file <- "./cervical_thoracic/all_bams/cov/rn5.chrom.sizes.txt"

bam.file.path <- "./cervical_thoracic/all_bams"
bam.files.names <- list.files(path = file.path(bam.file.path), pattern="bam$")

ConvertColAndMTChr <- function(filename) {
  
  bed.file<-read.table(filename, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  
  bed.file$V1[which(bed.file$V1 %in% "MT")] <- "M"
  
  
  bed.file$V1 <- paste0("chr", bed.file$V1)
  filename <- paste0(filename, "_changed_chr.bed")
  write.table(x=bed.file, file=filename, sep = "\t", col.names = FALSE, row.names = FALSE, qmethod = "double", quote = FALSE)
  return(filename)
}

for( file in bam.files.names ) {
  out.file <- paste0(file, "depth.txt")
  outfile.path.depth <- file.path(bam.file.path, "depths", out.file)
  depth.command <- paste("samtools view", file.path(bam.file.path, file), "| awk '{print $1}' | uniq | wc -l > ", outfile.path.depth)
  system(depth.command)
}



url <- "http://140.164.12.166/~FTP/UHN_time_course/coverage/"
coverage.filename <- "/media/dario/dati/time_course/cervical_thoracic/covs.txt"

for(file in bam.files.names) {
  out.file1 <- paste0(file, "_no_junct.bed")
  out.file.path1 <- file.path(bam.file.path, "cov", out.file1)
  command <- paste("bamToBed -i", file.path(bam.file.path, file), "-cigar | awk '$7 !~ /N/' > ", out.file.path1)
  # system(command)
  in.file.path1 <- ConvertColAndMTChr(out.file.path1)
  out.file2 <- paste0(out.file1, "_coverage.bedGraph")
  out.file.path2 <- file.path(bam.file.path, "cov", out.file2)
  command <- paste("genomeCoverageBed -bg -i", in.file.path1, "-g", genome.file, ">", out.file.path2)
  # system(command)
  out.file3 <- paste0(out.file1, "_coverage.bw")
  out.file.path3 <- file.path(bam.file.path, "cov", "bw", out.file3)
  command <- paste("bedGraphToBigWig", out.file.path2, genome.file, out.file.path3)
  system(command)
  ##add track header
  file.url <- paste0(url, out.file3)
  sample.name <- paste(strsplit(x = file, split = "_", fixed = TRUE)[[1]][1:3], collapse = " ")
  track.header <- paste0("track type=bigWig name=\"", sample.name, "\" description=\"", sample.name, "\" visibility=full color=200,100,0 altColor=0,100,200 priority=20 bigDataUrl=", file.url,"\n")
  # write(x=track.header, file=coverage.filename, append = TRUE)
  print(track.header)
}




for( file in bam.files.names ) {
  
  out.file1 <- paste0(file, "_no_junct.bed")
  out.file.path1 <- file.path(bam.file.path, "cov", out.file1)
  command <- paste("bamToBed -i", file.path(bam.file.path, file), "-cigar | awk '$7 !~ /N/' > ", out.file.path1)
  system(command)
  in.file.path1 <- ConvertColAndMTChr(out.file.path1)
  out.file2 <- paste0(out.file1, "_coverage.bedGraph")
  out.file.path2 <- file.path(bam.file.path, "cov", out.file2)
  command <- paste("genomeCoverageBed -bg -i", in.file.path1, "-g", genome.file, ">", out.file.path2)
  system(command)
  depth.command <- paste("samtools view", file.path(bam.file.path, file), "| awk '{print $1}' | uniq | wc -l")
  file.depth <- system(depth.command, intern=TRUE)
  
  bg.file <- read.table(file = out.file.path2, header = FALSE, sep = "\t", quote = "")
  
  bg.file$V4 <- bg.file$V4 * (1000000/as.numeric(file.depth))
  out.file3 <- paste0(out.file1, "_coverage_normalized.bedGraph")
  out.file.path3 <- file.path(bam.file.path, "cov", out.file3)
  write.table(x = bg.file, file = out.file.path3, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  out.file4 <- paste0(out.file1, "_normalized_coverage.bw")
  out.file.path4 <- file.path(bam.file.path, "cov", "bw", out.file4)
  command <- paste("bedGraphToBigWig", out.file.path3, genome.file, out.file.path4)
  system(command)
  
}
# 
# file.name <- "~/Scaricati/140.164.12.60.txt"
# 
# file.rows <- read.table(file.name, stringsAsFactors = FALSE)
# 
# strsplit(x = file.rows[,1], split = "/", fixed = TRUE)
# 
# bw.r.i<-grep(pattern = ".bw", x = file.rows[,1], fixed=TRUE)
# nomi.bw <- file.rows[bw.r.i,1]
# 
# fattori.bw.folder[duplicated(fattori.bw.folder)]
# splitted.names <- strsplit(x = nomi.bw, split = "/", fixed = TRUE)
# fattori.bw.folder <- unlist(lapply(splitted.names, function(x) {return(x[length(x)])}))
# length(fattori.bw.folder)
# 
# bigwig.file <- "noemi_bigwig.txt"
# for(i in 1:length(fattori.bw.folder)) {
#   track.header <- paste0("track type=bigWig name=\"", fattori.bw.folder[i], "\" description=\"", fattori.bw.folder[i], "\" visibility=full color=200,100,0 bigDataUrl=", nomi.bw[i],"\n")
#   write(x=track.header, file=bigwig.file, append = TRUE)
# }
# track.header <- paste0("track type=bigWig name=\"", sample.name, "\" description=\"", sample.name, "\" visibility=full color=200,100,0 altColor=0,100,200 priority=20 bigDataUrl=", file.url,"\n")

prev.wd <- getwd()
setwd("/media/dario/dati/time_course/cervical_thoracic/all_bams/depths/")
depths.files<-list.files("/media/dario/dati/time_course/cervical_thoracic/all_bams/depths/")

covs.data.frame <- data.frame()
for(file in depths.files) {
  clean.filename <- strsplit(x = file, split = ".", fixed = TRUE)[[1]][1]
  covs.data.frame <- rbind(covs.data.frame, cbind(clean.filename, readLines(file)))
}

write.table(x = covs.data.frame, file = "samples_unique_reads.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = c("sample_id", "unique_reads"))
setwd(prev.wd)

bam.file.path <- "/media/dario/dati/time_course/cervical_thoracic/all_bams"
bam.files.names <- list.files(path = file.path(bam.file.path), pattern="bam$")

reads.data.frame <- data.frame()
for( file in bam.files.names ) {
  out.file <- paste0(file, "reads.txt")
  outfile.path.reads <- file.path(bam.file.path, "total_reads", out.file)
  depth.command <- paste("samtools view", file.path(bam.file.path, file), "| awk '{print $1}' | wc -l > ", outfile.path.reads)
  system(depth.command)
  
  clean.filename <- strsplit(x = file, split = ".", fixed = TRUE)[[1]][1]
  reads.data.frame <- rbind(reads.data.frame, cbind(clean.filename, readLines(outfile.path.reads)))
}

write.table(x = reads.data.frame, file = "/media/dario/dati/time_course/cervical_thoracic/all_bams/samples_total_reads.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = c("sample_id", "total_reads"))

reads <- cbind(reads.data.frame, covs.data.frame[,2])

colnames(reads) <- c("sample_id", "total_reads", "unique_reads")
write.table(x = reads, file = "/media/dario/dati/time_course/cervical_thoracic/all_bams/samples_reads.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# prova.depth<-read.table("cervical_thoracic/all_bams/prova.depth.txt", header = TRUE, sep = "\t", quote = "")
# 
# sum(prova.depth$X1.1)
# 
# prova.depth<-read.table("cervical_thoracic/all_bams/aa.out", header = TRUE, sep = "\t", quote = "")
# 
# tail(prova.depth)
# unique(prova.depth$X1)
# # 
# # listed.names <- strsplit(x = bam.files.names, split = "_", fixed = TRUE)
# for(name in listed.names) {
#   print(paste(name[1:3], collapse = " "))
# }










# samtools view accepted_hits.bam | grep -w 'NH:i:1' | awk '{if ($6!~/N/) print}' | 
#   cat header.sam - | samtools view -Sb - | 
#   genomeCoverageBed -ibam stdin -g /home/solid/ColorspacePackages/bedtools-2.17.0/genomes/human.hg19.genome -bg > Reads_unique_no_junct.wiggle
