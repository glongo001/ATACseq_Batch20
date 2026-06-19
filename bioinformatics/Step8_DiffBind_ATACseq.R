#set working directory
setwd("/share/crsp/lab/blumberg/share/")

#load libraries
library(DiffBind) #differential binding analysis
library(DESeq2) #differential analysis

#sample information
bed_files <- list.files("broadPeaks/bedFiles", pattern = "\\.bed$", full.names = TRUE)
bam_files <- list.files("dedup_bam", pattern = "\\.dedup\\.bam$", full.names = TRUE)

#extract basename
bed_ids <- sub("\\.bed$", "", basename(bed_files))
bam_ids <- sub("\\.dedup\\.bam$", "", basename(bam_files))

#match BAMs to BEDs
bam_map <- setNames(bam_files, bam_ids)
matched_bams <- bam_map[bed_ids]

# check for mismatches
if(any(is.na(matched_bams))) {
  stop("Mismatch between BED and BAM files")
}

#extract sample info to get conditions
sample_info <- data.frame(
  sampleID = bed_ids,
  sex = substr(bed_ids, 4, 4),
  treatment = ifelse(grepl("DMSO", bed_ids), "DMSO", "TBT"),
  week = ifelse(grepl("w1", bed_ids), "w1", "w5")
)

#condition including sex, treatment and week
sample_info$condition <- paste0(sample_info$sex, sample_info$treatment, sample_info$week)
#iterate through sample list based on condition
replicates <- as.numeric(ave(bed_ids, sample_info$condition, FUN = seq_along))

# create empty DBA object
ATAC.bam <- NULL

# manually add each peak file + BAM
for (i in seq_along(bed_files)) {
  ATAC.bam <- dba.peakset(
    ATAC.bam,
    peaks = bed_files[i],
    sampID = bed_ids[i],
    condition = sample_info$condition[i],
    replicate = replicates[i],
    peak.caller = "bed",
    peak.format = "bed",
    bamReads = matched_bams[i]
  )
}

#inspect
dba.show(ATAC.bam)

#read count at each peak, create consensus peakset, access with dba.peakset or dba.peakcounts
ATAC.count <- dba.count(ATAC.bam)
#save object
saveRDS(ATAC.count, "bioinformatics/diffbind_counted.rds")

#plot affinity score-based correlation graphs
png("bioinformatics/plots/correlation_affinity.png", width = 800, height = 600)
plot(ATAC.count)
dev.off()

#establish contrasts based on condition
ATAC.contrast <- dba.contrast(ATAC.count, categories = DBA_CONDITION)
#show
dba.show(ATAC.contrast, bContrasts = TRUE)

#save image
save.image(file = "bioinformatics/DiffBind_ATACseq_Batch20.RData")

#Differential Analysis
#run analysis, get differential binding
ATAC.diff <- dba.analyze(ATAC.contrast, method = DBA_DESEQ2, bParallel = FALSE)
dba.show(ATAC.diff, bContrasts=TRUE)
#save object
saveRDS(ATAC.diff, "bioinformatics/diffbind_analyzed.rds")

#plot contrasted correlations
png("bioinformatics/plots/contrasted_correlations.png", width = 800, height = 600)
plot(ATAC.diff)
dev.off()

#get contrast ids
contrasts <- dba.show(ATAC.diff, bContrasts = TRUE)
idx <- which(
  (contrasts$Group == "fTBTw1" & contrasts$Group2 == "mDMSOw1") |
  (contrasts$Group == "mDMSOw1" & contrasts$Group2 == "fTBTw1")
)
#get differential binding sites and write to table
if(length(idx) == 1) {
  ATAC.DB <- tryCatch(dba.report(ATAC.diff, contrast = idx), error = function(e) NULL)
  if (!is.null(ATAC.DB)) {
    write.table(ATAC.DB, file = "bioinformatics/ATAC_DiffPeaks.txt", quote = FALSE, sep = "\t")
  }
}
#male week 1
idx.malew1 <- which(
  (contrasts$Group == "mTBTw1" & contrasts$Group2 == "mDMSOw1") |
  (contrasts$Group == "mDMSOw1" & contrasts$Group2 == "mTBTw1")
)
if(length(idx.malew1) == 1) {
  ATAC.DB.malew1 <- tryCatch(dba.report(ATAC.diff, contrast = idx.malew1), error = function(e) NULL)
  if (!is.null(ATAC.DB.malew1)) {
    write.table(ATAC.DB.malew1, file = "bioinformatics/ATAC_DiffPeaks_malew1.txt", quote = FALSE, sep = "\t")
  }
}
#male week 5
idx.malew5 <- which(
  (contrasts$Group == "mTBTw5" & contrasts$Group2 == "mDMSOw5") |
  (contrasts$Group == "mDMSOw5" & contrasts$Group2 == "mTBTw5")
)
if(length(idx.malew5) == 1) {
  ATAC.DB.malew5 <- tryCatch(dba.report(ATAC.diff, contrast = idx.malew5), error = function(e) NULL)
  if (!is.null(ATAC.DB.malew5)) {
    write.table(ATAC.DB.malew5, file = "bioinformatics/ATAC_DiffPeaks_malew5.txt", quote = FALSE, sep = "\t")
  }
}
#female week 1
idx.femalew1 <- which(
  (contrasts$Group == "fTBTw1" & contrasts$Group2 == "fDMSOw1") |
  (contrasts$Group == "fDMSOw1" & contrasts$Group2 == "fTBTw1")
)
if(length(idx.femalew1) == 1) {
  ATAC.DB.femalew1 <- tryCatch(dba.report(ATAC.diff, contrast = idx.femalew1), error = function(e) NULL)
  if (!is.null(ATAC.DB.femalew1)) {
    write.table(ATAC.DB.femalew1, file = "bioinformatics/ATAC_DiffPeaks_femalew1.txt", quote = FALSE, sep = "\t")
  }
}
#female week 5
idx.femalew5 <- which(
  (contrasts$Group == "fTBTw5" & contrasts$Group2 == "fDMSOw5") |
  (contrasts$Group == "fDMSOw5" & contrasts$Group2 == "fTBTw5")
)
if(length(idx.femalew5) == 1) {
  ATAC.DB.femalew5 <- tryCatch(dba.report(ATAC.diff, contrast = idx.femalew5), error = function(e) NULL)
  if (!is.null(ATAC.DB.femalew5)) {
    write.table(ATAC.DB.femalew5, file = "bioinformatics/ATAC_DiffPeaks_femalew5.txt", quote = FALSE, sep = "\t")
  }
}

#Plots
#function to plot safely, skip plots instead of crashing
safe_plot <- function(filename, expr) {
  png(filename, width = 800, height = 600)
  try(expr)
  dev.off()
}

#Venn Diagram
#minimum overlap of 0.33 makes parameters more lenient so for example if there are 3 replicates peak is kept if it's in at least 1 replicate
#ATACseq peaks be noisy
#subset before creating consensus object
w1_mask <- ATAC.bam$masks$fDMSOw1 |
  ATAC.bam$masks$fTBTw1 |
  ATAC.bam$masks$mDMSOw1 |
  ATAC.bam$masks$mTBTw1
ATAC.w1 <- dba(ATAC.bam, mask = w1_mask)
#consensus object week 1
ATAC.consensus.w1 <- dba.peakset(ATAC.w1, consensus = DBA_CONDITION, minOverlap = 0.33)
#plot
safe_plot("bioinformatics/plots/Venndiagram_w1.png", dba.plotVenn(ATAC.consensus.w1, ATAC.consensus.w1$masks$Consensus, main = "ATAC Overlap Peaks"))

#Venn diagram w5
w5_mask <- ATAC.bam$masks$fDMSOw5 |
  ATAC.bam$masks$fTBTw5 |
  ATAC.bam$masks$mDMSOw5 |
  ATAC.bam$masks$mTBTw5
ATAC.w5 <- dba(ATAC.bam, mask = w5_mask)
ATAC.consensus.w5 <- dba.peakset(ATAC.w5, consensus = DBA_CONDITION, minOverlap = 0.33)
#plot
safe_plot("bioinformatics/plots/Venndiagram_w5.png", dba.plotVenn(ATAC.consensus.w5, ATAC.consensus.w5$masks$Consensus, main = "ATAC Overlap Peaks"))

#MA plots to get accessibility, if >0 more accessible in TBT, <0 more accessible in DMSO
#male week 1
if(length(idx.malew1) == 1) {
  safe_plot("bioinformatics/plots/MAplot_malew1.png", dba.plotMA(ATAC.diff, contrast = idx.malew1))
}
#male week 5
if(length(idx.malew5) == 1) {
  safe_plot("bioinformatics/plots/MAplot_ATAC_malew5.png", dba.plotMA(ATAC.diff, contrast = idx.malew5))
}
#female week 1
if(length(idx.femalew1) == 1) {
  safe_plot("bioinformatics/plots/MAplot_ATAC_femalew1.png", dba.plotMA(ATAC.diff, contrast = idx.femalew1))
}
#female week 5
if(length(idx.femalew5) == 1) {
  safe_plot("bioinformatics/plots/MAplot_ATAC_femalew5.png", dba.plotMA(ATAC.diff, contrast = idx.femalew5))
}

#PCA plots, samples cluster by condition, sex and week
#if PC1 separates the treatment effect dominates, if PC2 sex effect dominates
png("bioinformatics/plots/PCAplot_ATAC.png", width = 800, height = 600)
dba.plotPCA(ATAC.diff, DBA_CONDITION, label = DBA_CONDITION)
dev.off()

#boxplots, shows count per sample and p-values (FDR) after deseq2 normalization
#medians should be aligned
#male week 1
if(length(idx.malew1) == 1) {
  safe_plot("bioinformatics/plots/boxplot_malew1.png", pvals.ATAC <- dba.plotBox(ATAC.diff, contrast = idx.malew1))
}
#male week 5
if(length(idx.malew5) == 1) {
  safe_plot("bioinformatics/plots/boxplot_malew5.png", pvals.ATAC <- dba.plotBox(ATAC.diff, contrast = idx.malew5))
}
#female week 1
if(length(idx.femalew1) == 1) {
  safe_plot("bioinformatics/plots/boxplot_femalew1.png", pvals.ATAC <- dba.plotBox(ATAC.diff, contrast = idx.femalew1))
}
#female week 5
if(length(idx.femalew5) == 1) {
  safe_plot("bioinformatics/plots/boxplot_femalew5.png", pvals.ATAC <- dba.plotBox(ATAC.diff, contrast = idx.femalew5))
}

#heatmaps, shows only significant differential peaks, should have clear separation between DMSO and TBT
#correlation values, based on peak counts across samples, for covariance/correlation analysis
#male week 1
if(length(idx.malew1) == 1) {
  safe_plot("bioinformatics/plots/heatmap_ATAC_differential_malew1.png", corvals.ATAC <- dba.plotHeatmap(ATAC.diff, contrast = idx.malew1, correlations = FALSE))
}
#male week 5
if(length(idx.malew5) == 1) {
  safe_plot("bioinformatics/plots/heatmap_ATAC_differential_malew5.png", corvals.ATAC <- dba.plotHeatmap(ATAC.diff, contrast = idx.malew5, correlations = FALSE))
}
#female week 1
if(length(idx.femalew1) == 1) {
  safe_plot("bioinformatics/plots/heatmap_ATAC_differential_femalew1.png", corvals.ATAC <- dba.plotHeatmap(ATAC.diff, contrast = idx.femalew1, correlations = FALSE))
}
#female week 5
if(length(idx.femalew5) == 1) {
  safe_plot("bioinformatics/plots/heatmap_ATAC_differential_femalew5.png", corvals.ATAC <- dba.plotHeatmap(ATAC.diff, contrast = idx.femalew5, correlations = FALSE))
}

#more heatmaps (unsupervised)
#for all replicates should cluster tightly and conditions should separate
#pearson correlation measures linear relationships
png("bioinformatics/plots/heatmap_ATAC_pearson.png", width = 800, height = 600)
try({dba.plotHeatmap(ATAC.count, distMethod = "pearson")})
dev.off()
#spearman correlation is rank-based and is more robust to outliers
png("bioinformatics/plots/heatmap_ATAC_spearman.png", width = 800, height = 600)
try({dba.plotHeatmap(ATAC.count, distMethod = "spearman")})
dev.off()

save.image(file = "bioinformatics/DiffBind_ATACseq_Batch20.RData")
