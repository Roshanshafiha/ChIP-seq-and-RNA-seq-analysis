#R version- 4.0.3

## loading packages

library(ChIPseeker)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)

data_1="ENCFF418TUX.bed.gz"
data_2="ENCFF791SNR.bed.gz"

files<-list("ENCFF418TUX.bed.gz","ENCFF791SNR.bed.gz")

names(files)[1] <- "ENCFF418TUX"
names(files)[2] <- "ENCFF791SNR"

#obtain the bed files in narrowPeak format (data_1)

file1_narrowPeak = data_1
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

gr1_narrowPeak <- import(file1_narrowPeak, format = "BED",
                        extraCols = extraCols_narrowPeak)

covplot(gr1_narrowPeak, weightCol="score")

#obtain the bed files in narrowPeak format (data_2)

file2_narrowPeak = data_2

gr2_narrowPeak <- import(file2_narrowPeak, format = "BED",
                         extraCols = extraCols_narrowPeak)

covplot(gr2_narrowPeak, weightCol="score")

#TSS regions for both the dataset

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))

#TSS regions for each plots

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")

#ChIP peak annotation comparison 

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnnoList)


genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")

dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")


#venn diagram that shows the overlap regions

genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)


#load the chipseq data this way

bed_file <- get_demo_file(format = "bed")


files<-read.table(gzfile("ENCFF418TUX.bed.gz"))
new_file<-import(files)

cols <- c("chrom","chrom-start","chrom-end","name","score","strand","signalValue","pValue","qValue","peak")
colnames(files) <- cols

covplot(peak, weightCol="V5")



