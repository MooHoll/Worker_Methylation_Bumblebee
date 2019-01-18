### Making a line graph of proportion methylation per gene and expression level ###

# Takes output from 2.proportion_meth_bases.R

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/Line graph")
library(ggplot2)
library(Hmisc)
library(reshape2)
library(car)


# Read in methylation proportion file and gene expression fpkm file
methylation<-read.csv(file="prop_meth_bases_per_gene_all_contexts_5%_reads_meth.csv", header = T)
gene_exp<-read.csv(file = "../all_genes_fpkm_values.csv", header=T)

# Sort out the meth file for easier use downstream with old code
meth_repro<-subset(methylation, methylation$staus=="repro")
meth_nonrepro<-subset(methylation, methylation$staus=="nonrepro")

colnames(meth_repro)<-c(".","gene","mean_meth_r","prop_meth_r","status","context")
colnames(meth_nonrepro)<-c(".","gene","mean_meth_nr","prop_meth_nr","status","context")

meth_repro<-meth_repro[,-c(1,3,5)]
meth_nonrepro<-meth_nonrepro[,-c(1,3,5)]

methylation<-merge(meth_nonrepro, meth_repro, by=c("gene","context"))

# Make one dataframe 
all_data<-merge(gene_exp,methylation, by="gene")
all_data<-all_data[,c(1,20,21,22,23,24)]


# Log transform the FPKM values
all_data$log_repro_fpkm<-log10(all_data$repro_fpkm_mean)
all_data$log_nonrepro_fpkm<-log10(all_data$nonrepro_fpkm_mean)


# Binning the exression data for counts (makes 100 bins for each)
all_data$repro_bins <- as.numeric(cut2(all_data$log_repro_fpkm, g=100))
all_data$nonrepro_bins<-as.numeric(cut2(all_data$log_nonrepro_fpkm, g=100))


# Make dataframes for each context (helps with older code)
cpg<-subset(all_data, all_data$context=="CpG")
chg<-subset(all_data, all_data$context=="CHG")
chh<-subset(all_data, all_data$context=="CHH")


# Make an average methylation value for all the genes that occur in each bin (do per context)
repro_FPKM_cpg<-as.data.frame(aggregate(cpg$prop_meth_r, by=list(cpg$repro_bins), mean))
colnames(repro_FPKM_cpg)<-c("expression_bin","r_mean_methylation")
plot(repro_FPKM_cpg$r_mean_methylation~repro_FPKM_cpg$expression_bin)

nonrepro_FPKM_cpg<-as.data.frame(aggregate(cpg$prop_meth_nr, by=list(cpg$nonrepro_bins), mean))
colnames(nonrepro_FPKM_cpg)<-c("expression_bin","nr_mean_methylation")
plot(nonrepro_FPKM_cpg$nr_mean_methylation~nonrepro_FPKM_cpg$expression_bin)



repro_FPKM_chg<-as.data.frame(aggregate(chg$prop_meth_r, by=list(chg$repro_bins), mean))
colnames(repro_FPKM_chg)<-c("expression_bin","r_mean_methylation")
plot(repro_FPKM_chg$r_mean_methylation~repro_FPKM_chg$expression_bin)

nonrepro_FPKM_chg<-as.data.frame(aggregate(chg$prop_meth_nr, by=list(chg$nonrepro_bins), mean))
colnames(nonrepro_FPKM_chg)<-c("expression_bin","nr_mean_methylation")
plot(nonrepro_FPKM_chg$nr_mean_methylation~nonrepro_FPKM_chg$expression_bin)



repro_FPKM_chh<-as.data.frame(aggregate(chh$prop_meth_r, by=list(chh$repro_bins), mean))
colnames(repro_FPKM_chh)<-c("expression_bin","r_mean_methylation")
plot(repro_FPKM_chh$r_mean_methylation~repro_FPKM_chh$expression_bin)

nonrepro_FPKM_chh<-as.data.frame(aggregate(chh$prop_meth_nr, by=list(chh$nonrepro_bins), mean))
colnames(nonrepro_FPKM_chh)<-c("expression_bin","nr_mean_methylation")
plot(nonrepro_FPKM_chh$nr_mean_methylation~nonrepro_FPKM_chh$expression_bin)


# Make a graph showing the trend using loess method
final_data_cpg<-merge(repro_FPKM_cpg, nonrepro_FPKM_cpg, by="expression_bin")
final_data_chg<-merge(repro_FPKM_chg, nonrepro_FPKM_chg, by="expression_bin")
final_data_chh<-merge(repro_FPKM_chh, nonrepro_FPKM_chh, by="expression_bin")

ggplot(data=final_data_cpg, aes(x=expression_bin))+
  geom_smooth(method="loess", size=5, aes(y=r_mean_methylation, colour="r_mean_methylation"#, linetype= "r_mean_methylation"
                                  ))+
  geom_smooth(method="loess", size=5, aes(y=nr_mean_methylation, colour="nr_mean_methylation"#,linetype = "nr_mean_methylation"
                                       ))+
  xlab("Expression Rank (Low to High)")+
  ylab("Proportion of Methylated Bases per Gene")+
  #scale_linetype_manual("", breaks=c("nr_mean_methylation","r_mean_methylation"),
   #                     values = c(1,3),
    #                    labels=c("Sterile","Reproductive"))+
  scale_colour_manual("", breaks=c("nr_mean_methylation","r_mean_methylation"),
                      values = c("dodgerblue","darkorchid"),
                      labels=c("Sterile","Reproductive"))+
  theme_bw() +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=28),
        legend.text=element_text(size=28))



############################################
# Statistical Analysis from Eamonn
############################################

##############
#CpG
data1<-melt(final_data_cpg,id=c("expression_bin"))
colnames(data1)<-c("expression_bin","treatment","mean_methylation")

model1<-lm(mean_methylation~treatment*expression_bin, data=data1)
#plot(model1)
model2<-lm(mean_methylation~treatment+expression_bin, data=data1)
#plot(model2)
anova(model1,model2)
#Analysis of Variance Table

#Model 1: mean_methylation ~ treatment * expression_bin
#Model 2: mean_methylation ~ treatment + expression_bin
#Res.Df        RSS Df   Sum of Sq     F Pr(>F)
#1    196 4.1806e-05                            
#2    197 4.1807e-05 -1 -2.1535e-10 0.001 0.9747

summary.lm(model2)
#Call:
#  lm(formula = mean_methylation ~ treatment + expression_bin, data = data1)

#Residuals:
#  Min         1Q     Median         3Q        Max 
#-1.588e-03 -3.871e-04  2.240e-06  3.660e-04  8.458e-04 

# Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)                   1.085e-03  7.328e-05  14.800   <2e-16 ***
#  treatmentnr_mean_methylation -1.794e-05  6.515e-05  -0.275    0.783    
#  expression_bin                3.045e-05  1.128e-06  26.979   <2e-16 ***

#Residual standard error: 0.0004607 on 197 degrees of freedom
#Multiple R-squared:  0.787,	Adjusted R-squared:  0.7849 
#F-statistic:   364 on 2 and 197 DF,  p-value: < 2.2e-16

drop1(model2, test="F")
#Single term deletions

#Model:
#  mean_methylation ~ treatment + expression_bin
#Df  Sum of Sq        RSS     AIC  F value Pr(>F)    
#<none>                       4.1807e-05 -3070.2                    
#treatment       1 1.6000e-08 4.1823e-05 -3072.1   0.0758 0.7833    
#expression_bin  1 1.5447e-04 1.9627e-04 -2762.9 727.8712 <2e-16 ***

##############
# CHG
data1<-melt(final_data_chg,id=c("expression_bin"))
colnames(data1)<-c("expression_bin","treatment","mean_methylation")

model1<-lm(mean_methylation~treatment*expression_bin, data=data1)
#plot(model1)
model2<-lm(mean_methylation~treatment+expression_bin, data=data1)
#plot(model2)
anova(model1,model2)
#Analysis of Variance Table
#Model 1: mean_methylation ~ treatment * expression_bin
#Model 2: mean_methylation ~ treatment + expression_bin
#Res.Df        RSS Df   Sum of Sq      F  Pr(>F)  
#1    196 3.4065e-07                                
#2    197 3.4581e-07 -1 -5.1575e-09 2.9675 0.08653 .


summary.lm(model2)
#Call:
#  lm(formula = mean_methylation ~ treatment + expression_bin, data = data1)

#Residuals:
#  Min         1Q     Median         3Q        Max 
#-1.276e-04 -2.822e-05 -2.850e-07  3.185e-05  1.287e-04 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                  4.102e-04  6.665e-06  61.555  < 2e-16 ***
#  treatmentnr_mean_methylation 3.391e-05  5.925e-06   5.724 3.83e-08 ***
#  expression_bin               2.062e-06  1.026e-07  20.089  < 2e-16 ***

#Residual standard error: 4.19e-05 on 197 degrees of freedom
#Multiple R-squared:  0.6889,	Adjusted R-squared:  0.6858 
#F-statistic: 218.2 on 2 and 197 DF,  p-value: < 2.2e-16

drop1(model2, test="F")
#Single term deletions

#Model:
#  mean_methylation ~ treatment + expression_bin
#Df  Sum of Sq        RSS     AIC F value    Pr(>F)    
#<none>                       3.4581e-07 -4029.1                      
#treatment       1 5.7510e-08 4.0332e-07 -4000.4  32.763 3.835e-08 ***
#  expression_bin  1 7.0838e-07 1.0542e-06 -3808.2 403.550 < 2.2e-16 ***

##############
# CHH
data1<-melt(final_data_chh,id=c("expression_bin"))
colnames(data1)<-c("expression_bin","treatment","mean_methylation")

model1<-lm(mean_methylation~treatment*expression_bin, data=data1)
#plot(model1)
model2<-lm(mean_methylation~treatment+expression_bin, data=data1)
#plot(model2)
anova(model1,model2)
#Analysis of Variance Table
#Model 1: mean_methylation ~ treatment * expression_bin
#Model 2: mean_methylation ~ treatment + expression_bin
#Res.Df        RSS Df   Sum of Sq      F Pr(>F)
#1    196 8.5658e-06                             
#2    197 8.6668e-06 -1 -1.0101e-07 2.3113   0.13


summary.lm(model2)
#Call:
#  lm(formula = mean_methylation ~ treatment + expression_bin, data = data1)

#Residuals:
#  Min         1Q     Median         3Q        Max 
#-6.416e-04 -1.337e-04  2.414e-05  1.544e-04  4.294e-04 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                  2.136e-03  3.336e-05  64.027  < 2e-16 ***
#  treatmentnr_mean_methylation 1.765e-04  2.966e-05   5.951  1.2e-08 ***
#  expression_bin               1.038e-05  5.138e-07  20.207  < 2e-16 ***

#Residual standard error: 0.0002097 on 197 degrees of freedom
#Multiple R-squared:  0.6925,	Adjusted R-squared:  0.6894 
#F-statistic: 221.9 on 2 and 197 DF,  p-value: < 2.2e-16

drop1(model2, test="F")
#Single term deletions

#Model:
#  mean_methylation ~ treatment + expression_bin
#Df  Sum of Sq        RSS     AIC F value  Pr(>F)    
#<none>                       8.6668e-06 -3384.9                    
#treatment       1 1.5581e-06 1.0225e-05 -3353.8  35.415 1.2e-08 ***
#  expression_bin  1 1.7963e-05 2.6630e-05 -3162.4 408.320 < 2e-16 ***

