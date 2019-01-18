################### Making methylation summary graphs ##################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Methylation/methylation_graphs")

library(plyr)
library(sqldf)
library(readr)
library(ggplot2)


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



# Read in data file
merged_final_meth_prop_with_context <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Methylation/methylation_graphs/merged_final_meth_prop_with_context.csv", 
                                                  "\t", escape_double = FALSE, trim_ws = TRUE)


# Make dataframes for each conext
cpg_data<-subset(merged_final_meth_prop_with_context, context=="CG")
chg_data<-subset(merged_final_meth_prop_with_context, context=="CHG")
chh_data<-subset(merged_final_meth_prop_with_context, context=="CHH")

# Make summary data for plotting
summary_CpG<-summarySE(cpg_data, measurevar = "proportion_meth", groupvars = c("feature","repro_status"))
summary_CHG<-summarySE(chg_data, measurevar = "proportion_meth", groupvars = c("feature","repro_status"))
summary_CHH<-summarySE(chh_data, measurevar = "proportion_meth", groupvars = c("feature","repro_status"))


# Making plots
ggplot(summary_CpG, aes(x=feature, y=proportion_meth, fill=repro_status))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=proportion_meth-se, ymax=proportion_meth+se),
                width=.2,
                position = position_dodge(.9))+
  xlab("Genomic Feature")+
  ylab("Average proportion of Methylated Reads in a CpG Context")+
  scale_fill_manual("", breaks=c("non_repro","repro"),
                    values = c("deepskyblue1","darkorchid"),
                    labels=c("Non-Reproductive", "Reproductive"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))


ggplot(summary_CHG, aes(x=feature, y=proportion_meth, fill=repro_status))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=proportion_meth-se, ymax=proportion_meth+se),
                width=.2,
                position = position_dodge(.9))+
  xlab("Genomic Feature")+
  ylab("Average proportion of Methylated Reads in a CHG Context")+
  scale_fill_manual("", breaks=c("non_repro","repro"),
                    values = c("deepskyblue1","darkorchid"),
                    labels=c("Non-Reproductive", "Reproductive"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))


ggplot(summary_CHH, aes(x=feature, y=proportion_meth, fill=repro_status))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=proportion_meth-se, ymax=proportion_meth+se),
                width=.2,
                position = position_dodge(.9))+
  xlab("Genomic Feature")+
  ylab("Average proportion of Methylated Reads in a CHH Context")+
  scale_fill_manual("", breaks=c("non_repro","repro"),
                    values = c("deepskyblue1","darkorchid"),
                    labels=c("Non-Reproductive", "Reproductive"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

