#--------------- EAMONN'S  -------------------------

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Transcription/GO_analysis")

#Get seqinr package. This lets you read the fasta files
#install.packages("seqinr") 
library("seqinr")


#Get the fasta file in and make it useable
genome<-read.fasta (file = "final_GFF_bter_1.0.fasta", as.string = TRUE, strip.desc=TRUE)
seq.data<-as.data.frame(do.call(rbind, genome))
seq.data<-data.frame(as.character(names(genome)),seq.data)
colnames(seq.data)=c("Genes","Sequences")

#Input your list of genes
query <- read.table("./diff_exon_GO/diff_exon_genes.csv", quote="\"", comment.char="")
colnames(query)<-"Genes"


#Selecting the node
output <- merge(query,seq.data,by="Genes")


#Getting it out of R
write.csv(output, file="diff_exon_genes_with_fasta.csv")
#Import this (not open) into excel
write.fasta(as.list(output$Sequences), output$Genes, open = "w", nbchar = 100,file.out="diff_exon_genes.fasta")
#The output should be a fasta file with multiple sequences
