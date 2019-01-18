####### Methylation deciles in differentially expressed genes ########

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/Binned methylation")

library(Hmisc)
library(ggplot2)
library(readr)
library(ggsignif)


# Read in files
diffexpLogFC <- read.csv("../no_j824_and_without_2_genes_all.csv", header=T)
total_meth <- read.csv("../methylation_means_per_gene.csv", header=T)

# Change all log2FoldChange values into posative values for plotting purposes, we don't
# care whether genes are upreg/downreg in R or NR just that they are differentially expressed,
# the further from 0 the FC is the more diff exp
diffexpLogFC$log2FoldChange<-abs(diffexpLogFC$log2FoldChange)


# Take mean meth value of repro and nonrepro
total_meth$mean_meth <- rowMeans(subset(total_meth, select = c(repro, nonrepro)), na.rm = TRUE)


# Merge meth and diff exp files
diffexp_meth <- merge(total_meth,diffexpLogFC,by="gene")


# Bin the mean meth column into 10 bins
diffexp_meth$mean_meth_bin <- as.numeric(cut2(diffexp_meth$mean_meth, g=10))


# Make a mean log2FoldChange (gene expression difference) for each methylation bin
final_data<-as.data.frame(aggregate(diffexp_meth$log2FoldChange, by=list(diffexp_meth$mean_meth_bin), mean))
colnames(final_data)<-c("meth_bin","logFC")
plot(final_data$logFC~final_data$meth_bin)


# Get a list of un-methylated genes to plot an extra bar, other than the bins
all_genes<-as.data.frame(diffexpLogFC[,1])
colnames(all_genes)<-c("gene")
# The total_meth file contains all genes with some level of methylation, so get a list
# of the diff exp genes that don't appear in the total_meth ergo have no methylation
non_methylated<-as.data.frame(all_genes[!(all_genes$gene %in% total_meth$gene),1])
colnames(non_methylated)<-c("gene")


# Make a new dataframe containing the information needed for the unmethylated genes
un<-merge(non_methylated, diffexpLogFC, by="gene")
un<-un[,c(1,3)]
list_unmethylated_genes<-as.data.frame(unique(un$gene))
colnames(list_unmethylated_genes)<-"gene"
final<-merge(list_unmethylated_genes,diffexpLogFC, by="gene")
final$mean_meth_bin<-0


# Cut down the final tables before merging 
new_diffexp_meth<-diffexp_meth[,c(1,13,18)]
final<-final[,c(1,3,8)]

final1<-rbind(final,new_diffexp_meth)


#### Define summary function (ref:http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


# Make a summary of the data for plotting
summary1<-summarySE(final1, measurevar="log2FoldChange", groupvars=c("mean_meth_bin"))
summary1<-summary1[1:11,]
summary1$mean_meth_bin<-factor(summary1$mean_meth_bin)

# plot binned mean meth against diff exp (super exciting, this matches Glastad et al. 2016)
ggplot(summary1, aes(x=mean_meth_bin, y=log2FoldChange))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=log2FoldChange-ci, ymax=log2FoldChange+ci),
                width=.2,
                position = position_dodge(.9))+
  xlab("DNA Methylation Decile")+
  ylab("Differential Expression (Absolute log2FoldChange)")+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20))


# Linear model to see if methylation bin can be predicted by fold change in expression
library(car)
head(final1)

model1<-lm(mean_meth_bin~log2FoldChange, data=final1)
summary(model1)
# log2FoldChange pval<2e-16

ggplot(data=final1, aes(x=mean_meth_bin))+
  geom_smooth(method="lm",aes(y=log2FoldChange))+
  xlab("Binned Methylation")+
  ylab("Differential Expression (Absolute log2FoldChange)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))


# ANOVA on the data before summary
final1$mean_meth_bin<-as.factor(final1$mean_meth_bin)
fit<-aov(log2FoldChange~mean_meth_bin, data=final1)
summary(fit)
# sig-difference between methylation bins in terms of diff expression level
# p-val<2e-16


# Post-hoc tukey test to see where the variation is
tukey.test<-TukeyHSD(fit)
tukey.test
## Bin 0, 3, 9 and 10 are significantly different from all others

