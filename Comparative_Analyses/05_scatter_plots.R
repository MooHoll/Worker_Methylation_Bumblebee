######## Differential expression (FPKM) and methylation data correlation ##########
# (produces scatter plots with a linear regression line)

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/Scatter graphs")

library(ggplot2)
library(readr)

# Read in the fpm gene expression data and the methylation data
all_genes_fpkm_values <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/all_genes_fpkm_values.csv")
methylation <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/methylation_means_per_gene.csv")

# merge diffexp and meth datasets
diffexp_and_meth <- merge(all_genes_fpkm_values, methylation, by='gene') #7907 genes

# Log the expression FPKM values for repro and nonrepro
diffexp_and_meth$logmean_nonrep_fpkm <- log(diffexp_and_meth$nonrepro_fpkm_mean)
diffexp_and_meth$logmean_rep_fpkm <- log(diffexp_and_meth$repro_fpkm_mean)

### Fitting linear model to data (differential expression and methylation level)

# Removing all NA/NAN/Inf 
finiteElements<-which(is.finite(diffexp_and_meth$repro*diffexp_and_meth$logmean_rep_fpkm*diffexp_and_meth$nonrepro*diffexp_and_meth$logmean_nonrep_fpkm))
finiteData<-diffexp_and_meth[finiteElements,] #7902 genes left

# Model for repro
fit_repro<-lm(repro~logmean_rep_fpkm,data=finiteData, na.action=na.exclude)
summary(fit_repro)
anova(fit_repro)

# Model for non-repro
fit_nonrepro<-lm(nonrepro~logmean_nonrep_fpkm,data=finiteData, na.action=na.exclude)
summary(fit_nonrepro)
anova(fit_nonrepro)

# Correlation and linerar model graphs
ggplot(data=finiteData, aes(x=logmean_nonrep_fpkm, y=nonrepro))+
  geom_point(colour="grey40")+
  # Add's a linear regression line to the data
  geom_smooth(method='lm',colour='black')+
  # Using lm modles above obtain R^2 value
  annotate("text",x=7,y=75, label="paste(italic(R) ^ 2, \" = 0.134\")",size=10,parse=T)+
  ylab("Percentage of CpGs Methylated")+
  xlab("Expression Level (log(FPKM))")+
  theme_bw()+
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=28),
        plot.title=element_text(size = 28))+
  scale_x_continuous(breaks = round(seq(min(-6), max(9), by = 2.5),1))+
  ylim(0,100)

ggplot(data=finiteData, aes(x=logmean_rep_fpkm, y=repro))+
  geom_point(colour="grey40")+
  # Add's a linear regression line to the data
  geom_smooth(method='lm', colour='black')+
  # Using lm modles above obtain R^2 value
  annotate("text",x=7,y=75, label="paste(italic(R) ^ 2, \" = 0.127\")",size=10,parse=T)+
  ylab("Percentage of CpGs Methylated")+
  xlab("Expression Level (log(FPKM))")+
  theme_bw()+
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=28),
        plot.title=element_text(size = 28))+
  scale_x_continuous(breaks = round(seq(min(-6), max(9), by = 2.5),1))+
  scale_y_continuous(breaks = round(seq(min(0), max(100), by= 25),1))




# Also try for merged repro and non-repro data
head(finiteData)
finiteData$methylation_both<-(finiteData$repro+finiteData$nonrepro)/2
finiteData$logmean_both_fpkm<-(finiteData$logmean_nonrep_fpkm+finiteData$logmean_rep_fpkm)/2

# Model these data (this is non-significant!)
fit_both<-lm(methylation_both~logmean_both_fpkm,data=finiteData, na.action=na.exclude)
summary(fit_both)
anova(fit_both)

ggplot(data=finiteData, aes(x=logmean_both_fpkm, y=methylation_both))+
  geom_point(colour="grey40")+
  # Add's a linear regression line to the data
  geom_smooth(method='lm',colour='black')+
  # Using lm modles above obtain R^2 value
  annotate("text",x=7,y=60, label="paste(italic(R) ^ 2, \" = 0.1316\")",size=10,parse=T)+
  ylab("Percentage of CpGs Methylated")+
  xlab("Expression Level (log(FPKM))")+
  theme_bw()+
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=28),
        plot.title=element_text(size = 28))+
  scale_x_continuous(breaks = round(seq(min(finiteData$logmean_both_fpkm), 
                                        max(finiteData$logmean_both_fpkm), by = 2.5),1))
