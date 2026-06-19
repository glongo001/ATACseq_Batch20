#set working directory
setwd("/share/crsp/lab/blumberg/share/")

#load libraries
library(TxDb.Mmusculus.UCSC.mm10.knownGene) #mm10 genome
library(org.Mm.eg.db) #mouse gene annotations
library(ChIPseeker) #peak annotation
library(clusterProfiler) #GO enrichment

#CHIPseeker
#load genome annotation
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#get promoter annotation
promoter <- getPromoters(TxDb = txdb, upstream = 10000, downstream = 10000)
#upload peak file
if(file.exists("bioinformatics/ATAC_DiffPeaks_malew1.txt")) {
  peak <- readPeakFile("bioinformatics/ATAC_DiffPeaks_malew1.txt", header=TRUE, sep="\t")
} else {
  stop("Peak file not found")
}
#peak annotation, match each peak to genomic features, compares peak coordinates to peak annotations
peak.annot <- annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")

#pie chart, shows genomic distribution
#40% promoter means strong transcriptional regulation, mostly intergenic means enhancer driven regulation
png("bioinformatics/plots/piechart_malew1.png", width = 800, height = 600)
plotAnnoPie(peak.annot)
dev.off()

#upset plot, shows peak overlaps between annotation categories
png("bioinformatics/plots/upsetplot_malew1.png", width = 800, height = 600)
upsetplot(peak.annot)
dev.off()

#coverage plot, shows tss enrichment
#peak at 0 means strong prmoter enrichment, flat means no TSS preference
png("bioinformatics/plots/coverageplot_malew1.png", width = 800, height = 1000)
covplot(peak, weightCol = "Fold")
dev.off()

#heatmap with tags
#tagMatrix transforms genomic coordinates into signal aligned to TSS
tagMatrix <- getTagMatrix(peak, windows = promoter)
#tagHeatmap, shows each row is one peak, each column is position relative to TSS
#bright band at center is peaks enriched at promoters, spread out is enhancer-like behavior
png("bioinformatics/plots/heatmaptags_ATAC.png", width = 800, height = 600)
try({tagHeatmap(tagMatrix)})
dev.off()

#save image
save.image(file = "bioinformatics/DiffBind_ATACseq_Batch20.RData")

#GO analysis, save as table
#remove NA and duplicates
genes <- unique(na.omit(peak.annot@anno$geneId))
#take gene list from peaks, finds enriched biological processes
if(length(genes) > 0) {
  bp <- enrichGO(genes, OrgDb = 'org.Mm.eg.db', ont = "BP", readable = TRUE)
  write.table(bp, file = "bioinformatics/DiffPeaksGO.txt", sep = "\t", quote = FALSE)
}

#output
write.table(as.data.frame(peak.annot), file = "bioinformatics/DiffPeaksAnnot.txt", sep = "\t")

#save final workspace
save.image(file = "bioinformatics/DiffBind_ATACseq_Batch20.RData")
