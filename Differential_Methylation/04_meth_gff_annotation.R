### Overlapping geneID and annotation with SMP positions using a databse 
# NOTE: need to run 'coverate2cytosine' from bismark with --CX flag to output all C context's not just CpG
# Also use sed to remove all " from the input annotation file

setwd("/scratch/monoallelic/hm257/repro_methylation/graphs")

install.packages("sqldf", repos= "https://www.stats.bris.ac.uk/R/")
library(sqldf)
library(readr)

# Main file from GFF containing gene with start and end position, had to edit file with sed before use to remove 
# quote marks and change all N/A to NA
genome_annotation<-read.csv.sql("FINAL_BIG_annotation_file_with_introns_no_quotes.csv",
                           sql ="select * from file", sep=",",header = T)
genome_annotation<-genome_annotation[,-c(1,5,8,11,12)]


# Read in the files generated from bismark coverage2cytosine after alignment

j1nr_all_Cs <- read_delim("../cov_files/sorted_j1nr_all_Cs.txt", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
j1r_all_Cs <- read_delim("../cov_files/sorted_j1r_all_Cs.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
j5nr_all_Cs <- read_delim("../cov_files/sorted_j5nr_all_Cs.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
j5r_all_Cs <- read_delim("../cov_files/sorted_j5r_all_Cs.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
j8nr_all_Cs <- read_delim("../cov_files/sorted_j8nr_all_Cs.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
j8r_all_Cs <- read_delim("../cov_files/sorted_j8r_all_Cs.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)


non_repro<-merge(j1nr_all_Cs, j5nr_all_Cs, by=c("X1","X2"))
non_repro<-merge(non_repro, j8nr_all_Cs, by=c("X1","X2"))
non_repro$non_meth_reads<-non_repro$X5.x+non_repro$X5.y+non_repro$X5
non_repro$meth_reads<-non_repro$X4.x+non_repro$X4.y+non_repro$X4
non_repro_reads<-non_repro[,c(1,2,6,18:19)]
colnames(non_repro_reads)<-c("chr","position","context","non_meth_reads","meth_reads")


repro<-merge(j1r_all_Cs, j5r_all_Cs, by=c("X1","X2"))
repro<-merge(repro, j8r_all_Cs, by=c("X1","X2"))
repro$non_meth_reads<-repro$X5.x+repro$X5.y+repro$X5
repro$meth_reads<-repro$X4.x+repro$X4.y+repro$X4
repro_reads<-repro[,c(1,2,6,18:19)]
colnames(repro_reads)<-c("chr","position","context","non_meth_reads","meth_reads")



# Add row of total reads and filter on minimum coverage per site
repro_reads$total_reads<-repro_reads$non_meth_reads+repro_reads$meth_reads
non_repro_reads$total_reads<-non_repro_reads$non_meth_reads+non_repro_reads$meth_reads

repro_reads_coverage<-subset(repro_reads, total_reads>10) #23579656
nonrepro_reads_coverage<-subset(non_repro_reads, total_reads>10) #22983177


# Proportion of reads that are methylated 
repro_reads_coverage$proportion_meth<-repro_reads_coverage$meth_reads/repro_reads_coverage$total_reads
nonrepro_reads_coverage$proportion_meth<-nonrepro_reads_coverage$meth_reads/nonrepro_reads_coverage$total_reads


# Help with the memory 
rm(j1nr_all_Cs, j1r_all_Cs, j5nr_all_Cs, j5r_all_Cs, j8nr_all_Cs, j8r_all_Cs, repro, non_repro, repro_reads, non_repro_reads)


# Now annotate the files with gene information (might be worth subsetting the genome_annotation file and making
# new annotated files for each annotation, i.e one for gene, one for exons etc, then they would be smaller and more
# manageable.)
output_repro <- sqldf("SELECT sample.chr,
      sample.position,
      sample.context,
      sample.proportion_meth,
      ga.chr,
      ga.start,
      ga.end,
      ga.feature,
      ga.geneID,
      ga.exon_number,
      ga.intron_number
      FROM repro_reads_coverage AS sample
      LEFT JOIN genome_annotation AS ga 
      ON sample.chr = ga.chr
      AND (sample.position >= ga.start AND sample.position <= ga.end)") 

output_repro<-subset(output_repro, !output_repro$geneID=="NA")
output_repro<-output_repro[,-c(5)]
output_repro<-output_repro[!duplicated(output_repro),]


output_nonrepro <- sqldf("SELECT sample.chr,
      sample.position,
      sample.context,
      sample.proportion_meth,
      ga.chr,
      ga.start,
      ga.end,
      ga.feature,
      ga.geneID,
      ga.exon_number,
      ga.intron_number
      FROM nonrepro_reads_coverage AS sample
      LEFT JOIN genome_annotation AS ga 
      ON sample.chr = ga.chr
      AND (sample.position >= ga.start AND sample.position <= ga.end)") 

output_nonrepro<-subset(output_nonrepro, !output_nonrepro$geneID=="NA")
output_nonrepro<-output_nonrepro[,-c(5)]
output_nonrepro<-output_nonrepro[!duplicated(output_nonrepro),]


# Make just one file for further analysis with repro_status as column
output_repro$repro_status<-"repro"
output_nonrepro$repro_status<-"non_repro"

final_output<-rbind(output_repro, output_nonrepro)
write.table(final_output, file="merged_final_meth_prop_with_context.csv", sep="\t", row.names=F, quote=F)


