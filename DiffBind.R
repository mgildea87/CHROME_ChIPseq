library(org.Hs.eg.db)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPQC)
library(DiffBind)
library(ChIPseeker)
library(plyr)
library(BiocParallel)
register(SerialParam())
library(parallel)

#Load fraction of reads in peaks file
files <- list.files('alignment/frag_len', full.names = T, recursive = F)
FRP <- read.table('FRP.txt', header = TRUE)
FRP <- data.frame(sample = FRP$Sample, FRP = FRP$FRP)
FRP$Condition <- c(rep("Control", 2), rep("KD", 3))
FRP$sample <- factor(FRP$sample, levels = FRP$sample)

#Generate sample sheet for ChIPQC and DIFFbind
samples <- data.frame(SampleID = FRP$sample, Tissue = rep(NA, 5), Factor = rep(NA, 5), Condition = rep(NA, 5), Replicate = c("1","2","1", "2", "3"))
bam_files_antibody <- list.files('alignment', pattern = 'Antibody_sorted.bam$', full.names = T, recursive = F)
samples$bamReads <- bam_files_antibody
samples$ControlID <- paste(FRP$Condition,"c",1:3, sep = "")
bam_files_con <- list.files('alignment', pattern ='Control_sorted.bam$', full.names = T, recursive = F)
samples$bamControl <- bam_files_con
peak_files <- list.files('peaks', pattern = "narrowPeak", recursive = F, full.names = T)
samples$Peaks <- peak_files
samples$PeakCaller <- rep("narrow", 5)
samples$Tissue <- c(rep("Control", 2), rep("KD", 3))
samples$Factor <- c(rep("Untreated", 5))
samples$Condition <- paste(samples$Tissue, samples$Factor, sep = "_")

#Generate diffbind object
diff_bind <- dba(sampleSheet = samples)

#Black list filter
diff_bind <- dba.blacklist(diff_bind, blacklist = T, greylist = F, cores = 40)

#Generate consensus peak set
diff_bind_consensus <- dba.peakset(diff_bind, consensus = DBA_CONDITION, minOverlap = .90)
diff_bind_consensus <- dba(diff_bind_consensus, mask = diff_bind_consensus$masks$Consensus, minOverlap = 1)
diff_bind_consensus_peaks <- dba.peakset(diff_bind_consensus, bRetrieve = TRUE)

#Count reads in peaks
diff_bind_consensus_counts <- dba.count(diff_bind, peaks = diff_bind_consensus_peaks, summits = 400, bParallel=FALSE)

#Normalize
diff_bind_consensus_counts <- dba.normalize(diff_bind_consensus_counts, normalize=DBA_NORM_LIB, library=DBA_LIBSIZE_FULL)

#Diff peak strength
dm <- dba.contrast(diff_bind_consensus_counts, group1=diff_bind_consensus_counts$masks$Control_Untreated, group2=diff_bind_consensus_counts$masks$KD_Untreated, name1 = "Control", name2 = "KD")
dm <- dba.analyze(dm, method = DBA_DESEQ2, bGreylist = F)
DESeq_report <- dba.report(dm, bNormalized = T, th = 1)

#Annotate peaks
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak_anno <- annotatePeak(DESeq_report, tssRegion=c(-2000, 2000), TxDb=txdb, annoDb="org.Hs.eg.db")
peak_anno <- as.GRanges(peak_anno)
promoters <- peak_anno[which(peak_anno$annotation == "Promoter (<=1kb)" | peak_anno$annotation == "Promoter (1-2kb)"),]
