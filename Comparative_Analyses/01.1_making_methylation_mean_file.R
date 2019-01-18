### Making a file with mean methylation values per gene from coverage data from Bismark ######
# From the cytosine2coverage output file from Bismark aliognment software. This script requires a lot of 
# computer power, runs overnight on BERTIE. It produces an output files with information for every 'C'
# in the genome. It includes, it's proportion of methylated reads and a labelled gene ID if applicable.

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/Line graph")

library(plyr)
library(dplyr)
library(sqldf)
library(readr)


# Read in file containing every gene start and end position aalong with gene name (made in script:
# getting_length_of_all_genes.R)
genome_annotation<-read.csv.sql("Bter_1.0_gene_length_without_quotes.csv",
                           sql ="select * from file", sep=",",header = T)
genome_annotation<-genome_annotation[,-c(1,6)]


# Read in the coverage files from Bismark, these are produced by running the 'coverage2cytosine' command on
# .cov output files from Bismark alignment
j1nr_all_Cs <- read_delim("./geneID_and_SMP_files/raw_coverage_files/subset_j1nr_all_Cs.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
j1r_all_Cs <- read_delim("./geneID_and_SMP_files/raw_coverage_files/subset_j1r_all_Cs.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
j5nr_all_Cs <- read_delim("./geneID_and_SMP_files/raw_coverage_files/subset_j5nr_all_Cs.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
j5r_all_Cs <- read_delim("./geneID_and_SMP_files/raw_coverage_files/subset_j5r_all_Cs.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
j8nr_all_Cs <- read_delim("./geneID_and_SMP_files/raw_coverage_files/subset_j8nr_all_Cs.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
j8r_all_Cs <- read_delim("./geneID_and_SMP_files/raw_coverage_files/subset_j8r_all_Cs.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


# Cut down file sizes to make merge quicker later 
j1nr_all_Cs<-j1nr_all_Cs[,-c(3,7)]
j1r_all_Cs<-j1r_all_Cs[,-c(3,7)]
j5nr_all_Cs<-j5nr_all_Cs[,-c(3,7)]
j5r_all_Cs<-j5r_all_Cs[,-c(3,7)]
j8nr_all_Cs<-j8nr_all_Cs[,-c(3,7)]
j8r_all_Cs<-j8r_all_Cs[,-c(3,7)]


# Add column names
colnames(j1nr_all_Cs)<-c("scaffold","SMP","meth_reads","unmeth_reads","context")
colnames(j1r_all_Cs)<-c("scaffold","SMP","meth_reads","unmeth_reads","context")
colnames(j5nr_all_Cs)<-c("scaffold","SMP","meth_reads","unmeth_reads","context")
colnames(j5r_all_Cs)<-c("scaffold","SMP","meth_reads","unmeth_reads","context")
colnames(j8nr_all_Cs)<-c("scaffold","SMP","meth_reads","unmeth_reads","context")
colnames(j8r_all_Cs)<-c("scaffold","SMP","meth_reads","unmeth_reads","context")


# Merge files to get control/treatment file
non_repro<-merge(j1nr_all_Cs, j5nr_all_Cs, by=c("scaffold","SMP"))
non_repro<-merge(non_repro, j8nr_all_Cs, by=c("scaffold","SMP"))
repro<-merge(j1r_all_Cs, j5r_all_Cs, by=c("scaffold","SMP"))
repro<-merge(repro, j8r_all_Cs, by=c("scaffold","SMP"))


# Make a total column of reads in the two files for methylated and unmethylated 
non_repro$unmeth_reads<-non_repro$unmeth_reads.x+non_repro$unmeth_reads.y+non_repro$unmeth_reads
non_repro$meth_reads<-non_repro$meth_reads.x+non_repro$meth_reads.y+non_repro$meth_reads
repro$unmeth_reads<-repro$unmeth_reads.x+repro$unmeth_reads.y+repro$unmeth_reads
repro$meth_reads<-repro$meth_reads.x+repro$meth_reads.y+repro$meth_reads


# Cut down files
non_repro<-non_repro[,c(1,2,9,10,11)]
repro<-repro[,c(1,2,9,10,11)]


# Remove any rows where the methylated and unmethylated reads = 0
non_repro_reads<-subset(non_repro, !unmeth_reads==0 | !meth_reads==0)
repro_reads<-subset(repro, !unmeth_reads==0|!meth_reads==0)


# Add row of total reads and filter on minimum coverage per site
repro_reads$total_reads<-repro_reads$unmeth_reads+repro_reads$meth_reads
non_repro_reads$total_reads<-non_repro_reads$unmeth_reads+non_repro_reads$meth_reads

repro_reads_coverage<-subset(repro_reads, total_reads>5) 
nonrepro_reads_coverage<-subset(non_repro_reads, total_reads>5) 


# Proportion of reads that are methylated 
repro_reads_coverage$proportion_meth<-repro_reads_coverage$meth_reads/repro_reads_coverage$total_reads
nonrepro_reads_coverage$proportion_meth<-nonrepro_reads_coverage$meth_reads/nonrepro_reads_coverage$total_reads


# Again cut down files for memory 
repro_reads_coverage<-repro_reads_coverage[,-c(3,4,6)]
nonrepro_reads_coverage<-nonrepro_reads_coverage[,-c(3,4,6)]


# Now annotate the files with gene information
output_repro <- sqldf("SELECT repro_reads_coverage.scaffold,
      repro_reads_coverage.SMP,
      repro_reads_coverage.context,
      repro_reads_coverage.proportion_meth,
      genome_annotation.gene,
      genome_annotation.start,
      genome_annotation.end,
      genome_annotation.scaffold
      FROM repro_reads_coverage AS repro_reads_coverage
      LEFT JOIN genome_annotation AS genome_annotation 
      ON repro_reads_coverage.scaffold = genome_annotation.scaffold
      AND (repro_reads_coverage.SMP >= genome_annotation.start AND repro_reads_coverage.SMP <= genome_annotation.end)") 

output_nonrepro <- sqldf("SELECT nonrepro_reads_coverage.scaffold,
      nonrepro_reads_coverage.SMP,
      nonrepro_reads_coverage.context,
      nonrepro_reads_coverage.proportion_meth,
      genome_annotation.gene,
      genome_annotation.start,
      genome_annotation.end,
      genome_annotation.scaffold
      FROM nonrepro_reads_coverage AS nonrepro_reads_coverage
      LEFT JOIN genome_annotation AS genome_annotation 
      ON nonrepro_reads_coverage.scaffold = genome_annotation.scaffold
      AND (nonrepro_reads_coverage.SMP >= genome_annotation.start AND nonrepro_reads_coverage.SMP <= genome_annotation.end)") 


# Remove any rows not annotaed by a gene
output_repro<-subset(output_repro, !output_repro$gene=="NA")
output_nonrepro<-subset(output_nonrepro, !output_nonrepro$gene=="NA")

# Cut down again
output_repro<-output_repro[,-c(6,7,8)]
output_nonrepro<-output_nonrepro[,-c(6,7,8)]

# Remove duplicate rows
output_repro<-output_repro[!duplicated(output_repro),]
output_nonrepro<-output_nonrepro[!duplicated(output_nonrepro),]


# Make just one file for further analysis with repro_status as column
output_repro$repro_status<-"repro"
output_nonrepro$repro_status<-"non_repro"

final_output<-rbind(output_repro, output_nonrepro)
write.table(final_output, file="merged_final_meth_prop_with_context.csv", sep="\t", row.names=F, quote=F)


### After this has run then:
# awk '$4 > 0 {print $0}' merged_final_meth_prop_with_context.csv > meth_SMPs.csv
### this leaves SMPs with any level of methylation, can then separate into three files
### for each C context with grep, then cut down the files leaving only gene ID, gene length, SMP and status
# cut -f 2,5,6,7 <file> > <output>

### NOTE: determining the cut off proportion of reads to make a base 'methylated' should be considered.
### Have ran with >0, 0.2, 0.1 and 0.05 to check (equateing to more than 0, 20%, 10% and 5% of reads being methylated )
