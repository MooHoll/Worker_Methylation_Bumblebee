### Alternative exon inclusion 

### Guide for DEXSeq: http://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.pdf

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Transcription/diff_exon_analysis")

#source("https://bioconductor.org/biocLite.R")
#biocLite("DEXSeq")
library(DEXSeq)
library(readr)
library(ggplot2)
library(ggrepel)

# Ran following to identify the place the python scripts from DEXSeq are kept and then 
# copied these scipts onto ALICE2: dexseq_count.py dexseq_prepare_annotation.py
pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )

# Before running the below script need to use the python scripts on ALICE2 to prepare files
# and count the number of reads over each exon: SCRIPT = python_dexseq.sh

### Moved j8_24 into another directory and then ran this script again, should run with 17 files in list
countFiles = list.files("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Transcription/diff_exon_analysis", pattern=".txt", full.names=TRUE)

flattenedFile = list.files("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Transcription/diff_exon_analysis", pattern="gff", full.names=TRUE)
sampleTable<- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Transcription/diff_exon_analysis/sample_information.csv")

# Also need to remove the j824 row from this
sampleTable<-subset(sampleTable, sampleTable$sample!='j8_24')

sampleTable<-as.data.frame(sampleTable)
rownames(sampleTable)<-c(sampleTable$sample)

dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

# rlog transform the counts
rld = rlog(dxd)

# PCA plot, it's a mess (colony separates)
data1 = plotPCA(rld, intgroup = c("condition"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))
ggplot(data1, aes(PC1, PC2, color=condition)) + geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #coord_fixed()+
  geom_text_repel(aes(label=name), size=6,show.legend=FALSE, point.padding = 0.5, box.padding = 0.25) +
  scale_colour_manual("", breaks=c("non_repro","repro"),
                      values = c("black","grey60"),
                      labels=c("Sterile","Reproductive"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20))



dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
plotDispEsts( dxd )

# This tests two models, full model includes the 'condition': R/NR and reduced leaves it out,
# then it see which fits better
dxd = testForDEU( dxd )

dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

# This makes a useable results table and draws an MA plot
dxr = DEXSeqResults( dxd )
head(dxr)
mcols(dxr)$description

# Shrink the results to improve the MA-plot
dxr1 <- lfcShrink(dxd, contrast=c("condition","repro","non_repro"), res=dxr)
plotMA( dxr1, cex=0.8, ylim=c(-5,5))

# Number of exonic regions significant with FDR of 10%
table ( dxr$padj < 0.1 ) # 83

significant_exons<-as.data.frame(subset(dxr, padj<0.1))
write.csv(significant_exons, file="significant_exon_useage_withoutj824.csv", row.names = FALSE)

# How many genes does this affect
table ( tapply( dxr$padj < 0.1, dxr$groupID, any ) ) # 57
gene_list<-unique(significant_exons$groupID)

#### Graphs ##########

# Make a plot that shows the exon with different significant useage 
plotDEXSeq( dxr, "LOC100644055", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)

# This graph also includes the annoated transcript models for that gene
plotDEXSeq( dxr, "LOC100650436", displayTranscripts=TRUE, legend=TRUE,
            cex.axis=1.2, cex=1.3, lwd=2 )

# Can also show the counts for the different samples
plotDEXSeq( dxr, "LOC100644055", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr, "LOC100650436", expression=T, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

# Graph that removes overall expression change levels 
# THIS COULD BE THE MOST APPROPRIATE GRAPH
pdf('graphs_without_j824.pdf') 
for (i in gene_list) {
plotDEXSeq( dxr, i, expression=FALSE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
}
dev.off()
