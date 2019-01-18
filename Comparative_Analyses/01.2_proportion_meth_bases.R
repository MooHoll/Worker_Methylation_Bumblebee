############ Proportion of Methylated Bases per Gene ###########

# Takes edited output from 1.making_methylation_mean_file.R 
# normalises the methylation levels for each 'C' context by taking into account the length of the gene

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/Line graph")

library(readr)

# From 1.making_methylation_mean_file.R cont...

# Read in file giving length of all genes, made in script: getting_length_of_all_genes.R
all_genes<-read.csv(file= "Bter_1.0_gene_length.csv")
all_genes<-all_genes[,c(2,6)]

# Read in SMP files
chg_info<-read.delim("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/Line graph/geneID_and_SMP_files/mirko_meth_CHG_only_SMPs_0.05_cut.csv", header=FALSE)
colnames(chg_info)<-c("meth_point", "geneID", "length","status")
chh_info<-read.delim("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/Line graph/geneID_and_SMP_files/mirko_meth_CHH_only_SMPs_0.05_cut.csv", header=FALSE)
colnames(chh_info)<-c("meth_point", "geneID", "length","status")
cg_info<-read.delim("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/Line graph/geneID_and_SMP_files/mirko_meth_CpG_only_SMPs_0.05_cut.csv", header=FALSE)
colnames(cg_info)<-c("meth_point", "geneID", "length","status")

# Cut up the dataframes into repro and nonrepro for easier counting
chg_info_repro<-subset(chg_info, chg_info$status=='repro')
chg_info_nonrepro<-subset(chg_info, chg_info$status=='non_repro')
chh_info_repro<-subset(chh_info, chh_info$status=='repro')
chh_info_nonrepro<-subset(chh_info, chh_info$status=='non_repro')
cg_info_repro<-subset(cg_info, cg_info$status=='repro')
cg_info_nonrepro<-subset(cg_info, cg_info$status=='non_repro')

# Count up the number of occurances of each gene name to give the total number of methylated C's per gene
table1<-lapply(chg_info_repro[-1], table)
gene_counts_chg_repro=data.frame(item=names(unlist(table1)),count=unlist(table1)[],
                       stringsAsFactors=FALSE)
rownames(gene_counts_chg_repro)=c()
colnames(gene_counts_chg_repro)=c("gene","meth_count")

table1<-lapply(chg_info_nonrepro[-1], table)
gene_counts_chg_nonrepro=data.frame(item=names(unlist(table1)),count=unlist(table1)[],
                                 stringsAsFactors=FALSE)
rownames(gene_counts_chg_nonrepro)=c()
colnames(gene_counts_chg_nonrepro)=c("gene","meth_count")


table1<-lapply(chh_info_repro[-1], table)
gene_counts_chh_repro=data.frame(item=names(unlist(table1)),count=unlist(table1)[],
                                 stringsAsFactors=FALSE)
rownames(gene_counts_chh_repro)=c()
colnames(gene_counts_chh_repro)=c("gene","meth_count")

table1<-lapply(chh_info_nonrepro[-1], table)
gene_counts_chh_nonrepro=data.frame(item=names(unlist(table1)),count=unlist(table1)[],
                                    stringsAsFactors=FALSE)
rownames(gene_counts_chh_nonrepro)=c()
colnames(gene_counts_chh_nonrepro)=c("gene","meth_count")


table1<-lapply(cg_info_repro[-1], table)
gene_counts_cg_repro=data.frame(item=names(unlist(table1)),count=unlist(table1)[],
                                 stringsAsFactors=FALSE)
rownames(gene_counts_cg_repro)=c()
colnames(gene_counts_cg_repro)=c("gene","meth_count")

table1<-lapply(cg_info_nonrepro[-1], table)
gene_counts_cg_nonrepro=data.frame(item=names(unlist(table1)),count=unlist(table1)[],
                                    stringsAsFactors=FALSE)
rownames(gene_counts_cg_nonrepro)=c()
colnames(gene_counts_cg_nonrepro)=c("gene","meth_count")


# Remove additional text not needed
gene_counts_cg_nonrepro$gene <- gsub('geneID.', '', gene_counts_cg_nonrepro$gene)
gene_counts_chg_nonrepro$gene <- gsub('geneID.', '', gene_counts_chg_nonrepro$gene)
gene_counts_chh_nonrepro$gene <- gsub('geneID.', '', gene_counts_chh_nonrepro$gene)
gene_counts_cg_repro$gene <- gsub('geneID.', '', gene_counts_cg_repro$gene)
gene_counts_chg_repro$gene <- gsub('geneID.', '', gene_counts_chg_repro$gene)
gene_counts_chh_repro$gene <- gsub('geneID.', '', gene_counts_chh_repro$gene)


# Diving the meth count by 3 to get an average as at the moment it's all repro/nonrepro together
gene_counts_cg_nonrepro$mean_meth<-(gene_counts_cg_nonrepro$meth_count)/3
gene_counts_chg_nonrepro$mean_meth<-(gene_counts_chg_nonrepro$meth_count)/3
gene_counts_chh_nonrepro$mean_meth<-(gene_counts_chh_nonrepro$meth_count)/3
gene_counts_cg_repro$mean_meth<-(gene_counts_cg_repro$meth_count)/3
gene_counts_chg_repro$mean_meth<-(gene_counts_chg_repro$meth_count)/3
gene_counts_chh_repro$mean_meth<-(gene_counts_chh_repro$meth_count)/3

# Adding the gene length to the dataframes
cg_nonrepro_length<-merge(all_genes, gene_counts_cg_nonrepro, by="gene")
chg_nonrepro_length<-merge(all_genes, gene_counts_chg_nonrepro, by="gene")
chh_nonrepro_length<-merge(all_genes, gene_counts_chh_nonrepro, by="gene")
cg_repro_length<-merge(all_genes, gene_counts_cg_repro, by="gene")
chg_repro_length<-merge(all_genes, gene_counts_chg_repro, by="gene")
chh_repro_length<-merge(all_genes, gene_counts_chh_repro, by="gene")

# Adding proportion of methylation per gene column
cg_nonrepro_length$proportion_meth_bases<-1-(cg_nonrepro_length$total_bases-cg_nonrepro_length$mean_meth)/cg_nonrepro_length$total_bases
chg_nonrepro_length$proportion_meth_bases<-1-(chg_nonrepro_length$total_bases-chg_nonrepro_length$mean_meth)/chg_nonrepro_length$total_bases
chh_nonrepro_length$proportion_meth_bases<-1-(chh_nonrepro_length$total_bases-chh_nonrepro_length$mean_meth)/chh_nonrepro_length$total_bases
cg_repro_length$proportion_meth_bases<-1-(cg_repro_length$total_bases-cg_repro_length$mean_meth)/cg_repro_length$total_bases
chg_repro_length$proportion_meth_bases<-1-(chg_repro_length$total_bases-chg_repro_length$mean_meth)/chg_repro_length$total_bases
chh_repro_length$proportion_meth_bases<-1-(chh_repro_length$total_bases-chh_repro_length$mean_meth)/chh_repro_length$total_bases

# Cut out unwanted info and add context's and repro status as columns
cg_nonrepro_length$staus<-"nonrepro"
chg_nonrepro_length$staus<-"nonrepro"
chh_nonrepro_length$staus<-"nonrepro"
cg_repro_length$staus<-"repro"
chg_repro_length$staus<-"repro"
chh_repro_length$staus<-"repro"

cg_nonrepro_length$context<-"CpG"
cg_repro_length$context<-"CpG"
chg_nonrepro_length$context<-"CHG"
chg_repro_length$context<-"CHG"
chh_nonrepro_length$context<-"CHH"
chh_repro_length$context<-"CHH"

cg_nonrepro_length<-cg_nonrepro_length[,-c(2,3)]
chg_nonrepro_length<-chg_nonrepro_length[,-c(2,3)]
chh_nonrepro_length<-chh_nonrepro_length[,-c(2,3)]
cg_repro_length<-cg_repro_length[,-c(2,3)]
chg_repro_length<-chg_repro_length[,-c(2,3)]
chh_repro_length<-chh_repro_length[,-c(2,3)]


# Make one output file for future use
methylation<-rbind(cg_nonrepro_length,chg_nonrepro_length,chh_nonrepro_length,cg_repro_length,chg_repro_length,chh_repro_length)
write.csv(methylation, file="prop_meth_bases_per_gene_all_contexts_5%_reads_meth.csv")
