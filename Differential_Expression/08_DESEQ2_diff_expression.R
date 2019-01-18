##### Differential Expression Analysis, Repro and NonRepro Workers: DESeq2
# see https://www.bioconductor.org/help/workflows/rnaseqGene/#aligning-reads-to-a-reference-genome

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Transcription/diff_exp_all_genes")

library(DESeq2)
library(reshape2)
library(pheatmap)
library(genefilter)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(PoiClaClu)
library(ReporteRsjars)
library(ggbeeswarm)
library(readr)
library(ggrepel)

#################################################################################

# Count table made in .R script "count_table_diff_exp.R"
countdata <- read.csv("counts_diff_exp.csv")
head(countdata)
# Remove these two genes as massive outliers in samples j8_30, need to rm for all samples
# for further analysis to work 
countdata1<-countdata[!(countdata$geneID== "LOC100652289"),]
countdata<-countdata1[!(countdata1$geneID=="LOC105666893"),]
# exclude theis sample as possibly a classifying mistake and it should be non_repro but labelled as a repro
countdata = subset(countdata, !sample=="j8_24") 

# subset out count table to get metadata (rows of metadata must == columns of count data, so number of samples must match)
coldata<-countdata[,c(2,3,4)]
head(coldata)
coldata1<-coldata[!duplicated(coldata), ]
nrow(coldata1)

# Make count table matrix
countmatr<-countdata[,c(2,5,6)]
countmatr<-reshape2::melt(countmatr)
countmatr1<-dcast(countmatr, geneID ~ sample)
head(countmatr1)
countmatr1<-countmatr1[-c(1:5),]# Remove ambiguous genes etc
row.names(countmatr1)<-(countmatr1$geneID)
countmatr1<-countmatr1[,-1]
countmatr<-as.matrix(countmatr1)
head(countmatr)

# make a DESeq object taking into account colony and spliting by repro status
dds= DESeqDataSetFromMatrix(countData = countmatr, colData = coldata1, design = ~ colony + status)

#################################################################################

# remove features with low counts
dds = dds[ rowMeans(counts(dds)) > 10, ] 
nrow(dds) # 9494 genes left

# rlog transform counts
rld = rlog(dds, blind=FALSE)

#################################################################################

# PCA plot with reproductive status coloured
data = plotPCA(rld, intgroup = c("status"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=status)) + geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()+
  geom_text_repel(aes(label=name), size=6,show.legend=FALSE, point.padding = 0.5, box.padding = 0.25) +
  scale_colour_manual("", breaks=c("non_repro","repro"),
                    values = c("red","blue"),
                    labels=c("Sterile","Reproductive"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
      axis.title=element_text(size=20),
      legend.text=element_text(size=20))
  

# PCA plot with colony coloured
data = plotPCA(rld, intgroup = c("colony"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=colony)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()+
  scale_colour_manual("", breaks=c("j1","j5","j8"),
                      values = c("red","green3","blue"),
                      labels=c("J1","J5","J8"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20))


# two 1st samples plotted against each other to check consistency (for rlog and log2) 
par( mfrow = c( 1, 2 ) )
dds = estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,c(2,4)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
plot(assay(rld)[,c(2,4)],
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
par( mfrow = c( 1, 1 ) )

#################################################################################

# estimate size factors = normalize for library size
dds = DESeq2::estimateSizeFactors(dds)
dds = estimateDispersions(dds)
plotDispEsts(dds, xlab= "Mean of Normalised Counts",
             ylab= "Dispersion", cex=1.0, cex.lab=1.45, cex.axis=1.45)

# check sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )
colnames(sampleDistMatrix) <- c("j1_17_sterile","j1_3_repro","j1_5_sterile","j1_6_repro","j1_7_repro","j1_8_sterile",
                                "j5_13_repro","j5_15_sterile","j7_33_repro","j5_34_sterile","j5_6_repro","j5_7_sterile",
                                "j8_12_repro","j8_18_sterile","j8_23_sterile","j8_29_sterile","j8_30_repro")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         show_colnames = T,
         fontsize = 18)

# check sample distances using the poisson method
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- c("j1_17_sterile","j1_3_repro","j1_5_sterile","j1_6_repro","j1_7_repro","j1_8_sterile",
                                   "j5_13_repro","j5_15_sterile","j7_33_repro","j5_34_sterile","j5_6_repro","j5_7_sterile",
                                     "j8_12_repro","j8_18_sterile","j8_23_sterile","j8_29_sterile","j8_30_repro")
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         show_colnames = T,
         fontsize = 18)

#################################################################################

# differential expression (used Benjamini-Hochberg adjustment and a p-val of <0.1)
dds = DESeq2::DESeq(dds, parallel=TRUE)
resultsNames(dds)

# comparisons
res=results(dds, contrast=c("status","repro","non_repro"))
summary(res) # up-reg = 282 genes and down-reg = 242 genes, outliers = 8

#shrink log2 fold change to make comparison easier
res_plot<-lfcShrink(dds, contrast=c("status","repro","non_repro"), res=res) 

#distribution of coefficents of the model
plotMA(res_plot, ylim=c(-5,5),cex=1.0, cex.lab=1.45, cex.axis=1.45, xlab='Mean of Normalised Counts',
       ylab='Log Fold Change')
#plotMA(res,  ylim=c(-5,5))

# plot of p-vals excluding genes with very small counts
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", xlab= "P-value", ylab="Frequency",
     cex.axis=1.45, cex.lab=1.45, main="")

resOrd=res[order(res$log2FoldChange),]
head(resOrd)
resOrd1<-as.data.frame(resOrd)
resOrd1$gene<-rownames(resOrd1)
rownames(resOrd1)<-c()

resOrd_significant<-subset(resOrd1, padj<0.05)
nrow(resOrd_significant)#334 (including j8_24 110)
#write.csv(as.data.frame(resOrd_significant), file="no_j824_and_without_2_genes_differentially_expressed_genes.csv")
#write.csv(as.data.frame(resOrd), file="no_j824_and_without_2_genes_all.csv")

upreg_repro_sig<-subset(resOrd_significant, log2FoldChange>0 )
nrow(upreg_repro_sig)#172 (including j8_24 56)
#write.csv(upreg_repro_sig, file="upreg_in_repro.csv")

upreg_nonrepro_sig<-subset(resOrd_significant, log2FoldChange<0 )
nrow(upreg_nonrepro_sig)#162 (including j8_24 54)
#write.csv(upreg_nonrepro_sig, file="upreg_in_nonrepro.csv")

# list of genes that are differential between repr and nonrepr 
n=100
topdiff_down = head(c(1:nrow(res))[order(res$log2FoldChange)],n) # rownames(res)
topdiff_up = tail(c(1:nrow(res))[order(res$log2FoldChange)],n)
topdiff = union(topdiff_down,topdiff_up)

################################################################################# 

# heatmap of top differential genes in abdomen of repr vs nonrepr workers
mat = assay(rld)[ topdiff, ]
mat = mat - rowMeans(mat)
df2 = as.data.frame(colData(rld)[,c("status"),drop=FALSE])
colnames(df2)<-c("Status")
df2$Status<-gsub("repro", "Reproductive", df2$Status)
df2$Status<-gsub("non_Reproductive", "Sterile", df2$Status)
pheatmap(mat, annotation_col=df2,
         show_rownames = F,
         fontsize = 16)

#################################################################################

# plot expression levels of some genes

# most downregulated gene in Repro vs Nonrepro
topGene = rownames(res)[which.min(res$log2FoldChange)] 
topGene #LOC100651038
data = plotCounts(dds, gene=topGene, intgroup=c("status"), returnData=TRUE)
ggplot(data, aes(x=status, y=count, fill=status)) + 
  geom_boxplot(outlier.color="red", position=position_dodge(width=0.7), width=0.5) + 
  xlab("Reproductive Status") + ylab("Normalized read count") + 
  ggtitle(topGene) + scale_y_log10() +
  theme(axis.text.x = element_blank())

# most upregulated in Repro vs nonrepro
topGene = rownames(res)[which.max(res$log2FoldChange)] 
topGene #LOC100650436
data = plotCounts(dds, gene=topGene, intgroup=c("status"), returnData=TRUE)
ggplot(data, aes(x=status, y=count, fill=status)) + 
  geom_boxplot(outlier.color="red", position=position_dodge(width=0.7), width=0.5) + 
  xlab("Reproductive Status") + ylab("Normalized read count") + 
  ggtitle(topGene) + scale_y_log10() +
  theme(axis.text.x = element_blank())

# 12 most downregulated genes in repro vs nonrepro
n=12
selGenes = head(rownames(res)[order(res$log2FoldChange)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("status"), returnData=TRUE))))
ggplot(data, aes(x=status, y=count, fill=status)) + 
  geom_boxplot(outlier.color="black", position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Reproductive Status") + ylab("Normalized read count") + 
  scale_y_log10() + 
  ggtitle("Top Downregulated Genes in Reproductive vs Sterile")+ 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))+
  scale_fill_manual("", breaks=c("non_repro","repro"),
                    values = c("grey","grey42"),
                    labels=c("Sterile","Reproductive"))


# 12 most upregulated genes in repro vs nonrepro
n=12
selGenes = head(rownames(res)[order(-res$log2FoldChange)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("status"), returnData=TRUE))))
ggplot(data, aes(x=status, y=count, fill=status)) + 
  geom_boxplot(outlier.color="black", position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Reproductive Status") + ylab("Normalized read count") + 
  scale_y_log10() + 
  ggtitle("Top Upregulated Genes in Repro vs Sterile") + 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))+
  scale_fill_manual("", breaks=c("non_repro","repro"),
                    values = c("grey","grey42"),
                    labels=c("Sterile","Reproductive"))

#################################################################################

# Obtain fpm (fragments per million mapped reads) values (could be useful for analysis with meth)
fpm_values<-fpm(dds)
head(fpm_values)
nrow(fpm_values) # 9494 genes from start analysis
fpm_values<-as.data.frame(fpm_values)
fpm_values$gene<-rownames(fpm_values)
rownames(fpm_values)<-c()
fpm_values$repro_fpm_mean<-apply(fpm_values[,c(2,4,5,7,9,11,13,17)],1,mean)
fpm_values$nonrepro_fpm_mean<-apply(fpm_values[,c(1,3,6,8,10,12,14,15,16)],1,mean)
write.csv(as.data.frame(fpm_values), file="all_genes_fpm_values.csv")

head(res)
signif_res<-subset(res, padj<0.05)
signif_res<- as.data.frame(signif_res)
signif_res$gene<-rownames(signif_res)
nrow(signif_res) #334

fpm_diff_exp_genes<-merge(signif_res, fpm_values, by = "gene")
nrow(fpm_diff_exp_genes)
head(fpm_diff_exp_genes)
write.csv(as.data.frame(fpm_diff_exp_genes), file="diff_exp_genes_fpm_values.csv")

### Also get fpkm values (normalised FPM's)

# Import gene length of every gene in Bter_1.0 
Bter_gene_length <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Transcription/diff_exp_all_genes/Bter_gene_length.csv")
Bter_gene_length<-Bter_gene_length[,c(2,6)]
Bter_gene_length<-as.data.frame(Bter_gene_length)
rownames(Bter_gene_length)<-Bter_gene_length$gene
Bter_gene_length <- Bter_gene_length[ order(row.names(Bter_gene_length)), ]

# Add gene length to the correct genes
Bter_gene_length <- Bter_gene_length[match(rownames(dds), Bter_gene_length$gene),]
Bter_gene_length <- Bter_gene_length[,2]
mcols(dds)$basepairs<-Bter_gene_length

# Obtain fpkm values (not all genes had a length imported so some lost)
fpkm_values<-fpkm(dds)
fpkm_values_with_na<-apply(fpkm_values, 1, function(x){any(is.na(x))})
fpkm_values <- fpkm_values[!fpkm_values_with_na,]

nrow(fpkm_values) #8196 genes left

fpkm_values<-as.data.frame(fpkm_values)
fpkm_values$gene<-rownames(fpkm_values)
rownames(fpkm_values)<-c()
fpkm_values$repro_fpkm_mean<-apply(fpkm_values[,c(2,4,5,7,9,11,13,17)],1,mean)
fpkm_values$nonrepro_fpkm_mean<-apply(fpkm_values[,c(1,3,6,8,10,12,14,15,16)],1,mean)
write.csv(as.data.frame(fpkm_values), file="all_genes_fpkm_values.csv")

fpkm_diff_exp_genes<-merge(signif_res, fpkm_values, by = "gene")
nrow(fpkm_diff_exp_genes)
head(fpkm_diff_exp_genes)
write.csv(as.data.frame(fpkm_diff_exp_genes), file="diff_exp_genes_fpkm_values.csv")
 
