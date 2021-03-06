### Making exon methylation level graphs ###


setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/methylation_over_exons")

library(ggplot2)
library(readr)
library(reshape2)

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


# Read in the exon information file, this was created by the script: mern_meth_gff_annotation.R
# grep was then used to take out only rows containing exon as a feature and only CpG positions (less than 200 non-CpG were removed)
# this file contains methylation values for EVERY CpG in the genome that overlaps an exon 
exon_info <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/methylation_over_exons/exon_information_CpG_only.csv", 
                                               "\t", escape_double = FALSE, col_names = FALSE, 
                                               trim_ws = TRUE)

# Cut down the file to useable information
exon_info<-exon_info[,-c(1,2,3,5,6,7,10)]
colnames(exon_info)<-c("meth_proportion","geneID","exon_number","status")

# Make a new column contianing geneID and exon_number as a unique identifier per exon
exon_info$exon_id<-paste(exon_info$geneID, "_",exon_info$exon_number)
exon_info$exon_id <- gsub('\\s+', '', exon_info$exon_id)
exon_info<-exon_info[,-c(2,3)]

# Take the mean methylation level of every exon
mean_exon_meth<-aggregate( meth_proportion ~ exon_id+status, exon_info, mean )
rm(exon_info)


# Read in the list of genes containing diff exon expression generated by DEXseq, 86 total exons diff expressed
diff_exp_exons<- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/methylation_over_exons/significant_exon_useage_withoutj824.csv")

# Cut down to just useable information
diff_exp_exons<-diff_exp_exons[,c(1,2,10)]

# Make columns compatible with previous data frame
diff_exp_exons$featureID <- gsub('E0', '', diff_exp_exons$featureID)
diff_exp_exons$featureID<-as.numeric(diff_exp_exons$featureID)

colnames(diff_exp_exons)<-c("geneID","exon_number","Log2FC_repro_to_nonrepro")

diff_exp_exons$exon_id<-paste(diff_exp_exons$geneID, "_",diff_exp_exons$exon_number)
diff_exp_exons$exon_id <- gsub('\\s+', '', diff_exp_exons$exon_id)
diff_exp_exons<-diff_exp_exons[,-c(1,2)]


# Merge everything to get an exon list with methylation level and differential expression status
DEE_with_meth<-merge(diff_exp_exons, mean_exon_meth, by="exon_id")
DEE_with_meth$Log2FC_repro_to_nonrepro<-as.numeric(DEE_with_meth$Log2FC_repro_to_nonrepro)
DEE_with_meth$meth_proportion<-as.numeric(DEE_with_meth$meth_proportion)
DEE_with_meth$status<-as.factor(DEE_with_meth$status)

# Make new column so plotting is easier
DEE_with_meth$exp_and_status[DEE_with_meth$Log2FC_repro_to_nonrepro<0 & DEE_with_meth$status=="repro" ]<-"skipped"
DEE_with_meth$exp_and_status[DEE_with_meth$Log2FC_repro_to_nonrepro>0 & DEE_with_meth$status=="non_repro"]<-"skipped"
DEE_with_meth$exp_and_status[DEE_with_meth$Log2FC_repro_to_nonrepro>0 & DEE_with_meth$status=="repro" ]<-"included"
DEE_with_meth$exp_and_status[DEE_with_meth$Log2FC_repro_to_nonrepro<0 & DEE_with_meth$status=="non_repro"]<-"included"

# Make proportion into percentage for comparison with other graphs
DEE_with_meth$meth_proportion<-DEE_with_meth$meth_proportion*100

summary1<-summarySE(DEE_with_meth, measurevar="meth_proportion", groupvars=c("status","exp_and_status"))

# Bar plot with confidence intervals
ggplot(summary1, aes(x=exp_and_status, y=meth_proportion, fill=status))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=meth_proportion-ci, ymax=meth_proportion+ci),
                width=.2,
                position = position_dodge(.9))+
  xlab("Exon Status")+
  ylab("Average CpGs Methylated (%)")+
 # ggtitle("Methylation Levels of Differentially Expressed Exons")+  
  scale_fill_manual("", breaks=c("non_repro","repro"),
                    values = c("grey","grey42"),
                    labels =c("Sterile","Reproductive"))+
  scale_x_discrete(limits=c("included","skipped"),
                   labels=c("Included","Skipped"))+
  theme_bw()+
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=28),
        legend.text=element_text(size=28))

t.test(DEE_with_meth$meth_proportion~DEE_with_meth$exp_and_status)
# t = 0.10862, df = 139.47, p-value = 0.9137





