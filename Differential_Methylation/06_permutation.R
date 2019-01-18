#!/bin/bash

#PBS -N permutation
#PBS -l walltime=25:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
#cd $PBS_O_WORKDIR 

# Load software needed
#module load R/3.5.1

#Rscript permutation_with_common.R




library(methylKit)
library(readr)
library(stringr)
library(rlist)

#covariates <- data.frame(coplony=c("J1","J1","J5","J5","J8","J8"))
nsims=10000
diffs <- rep(NA, nsims)

common_sites <- rep(NA, nsims)
sites_to_check <- read.csv("DMRs_min10percentDiff_qval0.05_MSCfiltering.csv")
sites_to_check$chrBase <- paste(sites_to_check$chr,".",sites_to_check$start, sep="")

for (i in 1:nsims){
# Read in subsetted input files
file.list <- list.files("./", pattern = "*_subsetted_final.txt$")
data1<- lapply(file.list, read.delim, sep="\t", colClasses =
                 c("character","character","integer","character","integer","numeric",
                   "numeric"))
names(data1)<-c("j1nr","j1r","j5nr","j5r","j8nr","j8r")

# Make a column in each one that keeps the coverage, freqC and freqT together
for(j in seq_along(data1)){
  data1[[j]]$merged_info <- paste(data1[[j]]$coverage, 
                                  data1[[j]]$freqC,data1[[j]]$freqT,
                                  sep=" ")
}

# Put all the data in one dataframe
all_data <- list.cbind(data1)

# Take only the columns which contain the merged coverage/freqC/freqT info for each sample
for_shuffling <- all_data[,c("j1nr.merged_info","j1r.merged_info","j5nr.merged_info","j5r.merged_info","j8nr.merged_info","j8r.merged_info")]

# Take each row and shuffle this information between the 6 samples 
for (j in 1:nrow(for_shuffling)){
  for_shuffling[j,] <- sample(for_shuffling[j,])
}

# Add the shuffled information back to the oringinal position information
shuffled <- cbind(all_data[,1:4], for_shuffling)

# Re-make 6 samples. Each position now contains the same methylation information, however
# that information, per site, now belongs to a random sample 
j1nr <- shuffled[,c(1:4,5)]
j1r <- shuffled[,c(1:4,6)]
j5nr <- shuffled[,c(1:4,7)]
j5r <- shuffled[,c(1:4,8)]
j8nr <- shuffled[,c(1:4,9)]
j8r <- shuffled[,c(1:4,10)]

new.list <- list(j1nr, j1r, j5nr, j5r, j8nr, j8r)

# Make the final files, by splitting the merged coverage information back into individial
# columns and writing these new shuffled files out 
for(j in seq_along(new.list)){
  colnames(new.list[[j]]) <- c("chrBase","chr","base","strand","merged_info")
  new.list[[j]]<- cbind(new.list[[j]],(str_split_fixed(new.list[[j]]$merged_info, " ", 3)))
  new.list[[j]]<- new.list[[j]][-5]
  colnames(new.list[[j]])[5:7] <- c("coverage","freqC","freqT")
  myfile <- file.path("./", paste0(j,"_","shuffled.txt"))
  write.table(new.list[[j]], file=myfile, quote=F, sep="\t", row.names=F)
}

# Run a differenial methylation analysis on these shuffled files 
file.listA <- list("1_shuffled.txt","2_shuffled.txt",
                   "3_shuffled.txt","4_shuffled.txt",
                   "5_shuffled.txt","6_shuffled.txt")
raw_data <- methRead(file.listA,
                     sample.id = list("j1nr","j1r","j5nr","j5r","j8nr","j8r"),
                     treatment = c(1,0,1,0,1,0),
                     assembly="bter_1.0", 
                     context="CpG")

meth_all_data <- unite(raw_data)
diff_meth <- calculateDiffMeth(meth_all_data, covariates=covariates, mc.cores=7)
diff_meth_10 <- getMethylDiff(diff_meth, difference=10, qvalue=0.05)

# Count the number of differentially methylated sites and write to this variable 
diffs[i] <- nrow(diff_meth_10)

# Work out number of sites in common with original diff meth sites and write to variable 
diff_meth_data <- getData(diff_meth_10)
diff_meth_data$chrBase <- paste(diff_meth_data$chr,".",diff_meth_data$start, sep="")
check <- merge(sites_to_check, diff_meth_data, by="chrBase")
common_sites[i] <- nrow(check)

}

write.csv(diffs, file="permutation_results_diff_meth_sites_2nd.csv")
write.csv(common_sites, file="permutation_results_common_sites_2nd.csv")

pdf("histogram_test.pdf")
hist(diffs)
abline(v=626, col="red", lwd=2)
dev.off()
