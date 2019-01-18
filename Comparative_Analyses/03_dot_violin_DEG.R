################################################################################
## Making Bar plots for methylation levels in diff exp genes/isoforms and non-diff exp genes/isoforms ###
################################################################################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/Bar plots")

library(ggplot2)
library(ggsignif)
library(ggthemes)
library(plyr)
library(reshape2)
library(car)


# Read in files 
methylation<-read.csv(file="../methylation_means_per_gene.csv", header=T)
methylation$overall_mean<-(methylation$nonrepro+methylation$repro)/2
methylation<-methylation[,c(2,9:11)]

diff_exp_genes<-read.csv(file="../diff_exp_genes_fpkm_values.csv", header=T)
diff_exp_genes<-as.data.frame(diff_exp_genes[,2])
colnames(diff_exp_genes)<-c("gene")


# Make a list of non-differentially expressed genes 
meth_genes<-as.data.frame(methylation[,1])
colnames(meth_genes)<-c("gene")
non_DEG<-as.data.frame(meth_genes[!(meth_genes$gene %in% diff_exp_genes$gene),1])
colnames(non_DEG)<-c("gene")
rm(meth_genes)


# Merge to get average methylation in each expression list per gene
DEG_with_meth<-merge(diff_exp_genes, methylation, by="gene")
DEG_with_meth$name<-"DEG"

non_DEG_with_meth<-merge(non_DEG, methylation, by="gene")
non_DEG_with_meth<-non_DEG_with_meth[!is.na(non_DEG_with_meth$overall_mean),]
non_DEG_with_meth$name<-"non_DEG"


# Get methylation level from each gene set
DEG_mean<- mean(DEG_with_meth$overall_mean) #5.19
non_DEG_mean<-mean(non_DEG_with_meth$overall_mean) #8.25


# Get methylation level from each gene set for repro and non-repro
repro_DEG_mean<-mean(DEG_with_meth$repro) #5.19
nonrepro_DEG_mean<-mean(DEG_with_meth$nonrepro) #5.17
repro_nonDEG_mean<-mean(non_DEG_with_meth$repro) #8.26
nonrepro_nonDEG_mean<-mean(non_DEG_with_meth$nonrepro) #8.24


# Need to merge the dataframes
final_data<-rbind(DEG_with_meth,non_DEG_with_meth)
melted_data<-melt(final_data, id=c("gene","name"))
colnames(melted_data)<-c("gene","name","status","methylation")

## Dot plot
dot_data<-subset(melted_data, !status=='overall_mean')
dot_data_DEG<-subset(dot_data, name=="DEG")
DEG_repro<-subset(dot_data_DEG, status=="repro")
DEG_nonrepro<-subset(dot_data_DEG, status=="nonrepro")
dot_data_nonDEG<-subset(dot_data, name=="non_DEG")
nonDEG_repro<-subset(dot_data_nonDEG, status=="repro")
nonDEG_nonrepro<-subset(dot_data_nonDEG, status=="nonrepro")

DEG_repro$name<-"DEG_repro"
DEG_repro<-DEG_repro[,c(2,4)]
DEG_nonrepro$name<-"DEG_nonrepro"
DEG_nonrepro<-DEG_nonrepro[,c(2,4)]
nonDEG_repro$name<-"nonDEG_repro"
nonDEG_repro<-nonDEG_repro[,c(2,4)]
nonDEG_nonrepro$name<-"nonDEG_nonrepro"
nonDEG_nonrepro<-nonDEG_nonrepro[,c(2,4)]

dot_data1<-rbind(DEG_repro,DEG_nonrepro,nonDEG_repro,nonDEG_nonrepro)
dot_data1$status<-dot_data1$name
dot_data1$status<-gsub("nonDEG_","", dot_data1$status)
dot_data1$status<-gsub("DEG_","", dot_data1$status)
dot_data1$log_meth<-log(dot_data1$methylation)

condifence_intervals <- function(x) {
  m <- mean(x)
  sd1 <- sd(x)
  n <- length(x)
  error <- qnorm(0.95)*sd1/sqrt(n)
  ymin <- m-error
  ymax <- m+error
  return(c(y=m,ymin=ymin,ymax=ymax))
}

ggplot(dot_data1, aes(x=name, y=log_meth, fill=status)) + 
  geom_violin(trim = FALSE)+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.025, dotsize = 0.3)+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  xlab("Gene Group")+
  ylab("log(Average Percentage CpGs Methylated)")+
  scale_fill_manual("", breaks=c("repro", "nonrepro"),
                    values = c("grey", "grey35"),
                    labels= c("Repro", "Sterile"))+
  scale_x_discrete(limits=c("DEG_nonrepro", "DEG_repro","nonDEG_nonrepro", "nonDEG_repro"),
                   labels=c("DEG", "DEG", "Non-DEG", "Non-DEG"))+
  theme_bw()+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        legend.text=element_text(size=28))

#ggplot(dot_data1, aes(x=name, y=methylation)) + 
#  geom_violin(trim = FALSE)+
#  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.5, dotsize = 0.5)+
#  stat_summary(fun.data=condifence_intervals, color="blue")

#ggplot(dot_data1, aes(x=name, y=log_meth)) + 
#  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03, dotsize = 0.3)+
#  stat_summary(fun.data=mean_sdl, fun.args = list(mult=2), 
#               geom="pointrange", color="red")

# Do a t-test on the two groups
t.test(final_data$overall_mean~final_data$name) 

#Welch Two Sample t-test

#data:  final_data$overall_mean by final_data$name
#t = -5.3777, df = 295.56, p-value = 1.535e-07
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -4.179825 -1.940162
#sample estimates:
#  mean in group DEG mean in group non_DEG 
#5.192519              8.252513 

## ---------------------------------------------------
# Pretty sure I shouldn't be doing a t-test as doesn't take repro status into account
## ---------------------------------------------------

# Make a variable which is DEG/nonDEG
head(dot_data1)
dot_data1$name <- gsub("DEG.*","DEG", dot_data1$name)

# ...looking online a factorial anova?
# (http://rcompanion.org/handbook/G_09.html)

interaction.plot(x.factor     = dot_data1$name,
                 trace.factor = dot_data1$status, 
                 response     = dot_data1$methylation, 
                 fun = mean,
                 type="b",
                 col=c("black","red"),  
                 pch=c(19, 17),
                 leg.bty = "o")

model = lm(methylation ~ name*status,
           data = dot_data1)
Anova(model,
      type = "II") # there is no interaction effect so type II appropriate?

#Anova Table (Type II tests)

#Response: methylation
#Sum Sq    Df F value    Pr(>F)    
#name           5075     1 45.7922 1.348e-11 ***
#status            2     1  0.0206    0.8859    
#name:status       0     1  0.0005    0.9818    
#Residuals   2442314 22038

# If the above is the right thing to do then following Eamonn's previous code ...


model1<-lm(methylation~name*status, data=dot_data1)
model2<-lm(methylation~name+status, data=dot_data1)

anova(model1,model2) # This tell us the interaction isn't significant so proceed to type II anova

#Analysis of Variance Table

#Model 1: methylation ~ name * status
#Model 2: methylation ~ name + status
#Res.Df     RSS Df Sum of Sq     F Pr(>F)
#1  22038 2442314                          
#2  22039 2442314 -1 -0.057641 5e-04 0.9818

summary.lm(model2) # Get the same result as above code

#Call:
 # lm(formula = methylation ~ name + status, data = dot_data1)

#Residuals:
#  Min     1Q Median     3Q    Max 
#-8.263 -7.833 -5.552  5.784 91.737 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  5.18234    0.45204  11.464  < 2e-16 ***
#  namenonDEG   3.05999    0.45218   6.767 1.35e-11 ***
#  statusrepro  0.02036    0.14181   0.144    0.886    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 10.53 on 22039 degrees of freedom
#Multiple R-squared:  0.002074,	Adjusted R-squared:  0.001984 
#F-statistic: 22.91 on 2 and 22039 DF,  p-value: 1.153e-10

drop1(model2, test="F") 

#Single term deletions

#Model:
#  methylation ~ name + status
#Df Sum of Sq     RSS    AIC F value    Pr(>F)    
#<none>              2442314 103774                      
#name    1    5074.8 2447389 103818 45.7942 1.346e-11 ***
#  status  1       2.3 2442317 103772  0.0206    0.8859    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
