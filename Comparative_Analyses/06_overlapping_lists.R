### --------------------------------------
## Overlapping gene and GO lists
### --------------------------------------

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Meth_vs_Exp/LOC_GO_overlap")

library(readr)
library(rowr)

### --------------------------------------
## Gene Lists and Hypergeometric
### --------------------------------------

# Hypergeometric test:
# Population size: 11030 (grep gff file)
# Number of successes in population: methylation number of genes
# Sample size: number of exon/diff expressed genes
# Number of successes in sample: overlapping number 

# lower-tail = T checks if under-represented compared to chance
# lower-tail = F checks if over-represented compared to chance

# Methylation gene lists
diff_meth_genes <- read_csv("gene_lists/diff_meth_genes.csv") #478
gene_hyper_repro_all <- read_csv("gene_lists/gene_hyper_repro_all.csv") #253
gene_hyper_sterile_all <- read_csv("gene_lists/gene_hyper_sterile_all.csv") #276

# Expression gene lists
upreg_in_nonrepro <- read_csv("gene_lists/upreg_in_nonrepro.csv") #162
upreg_in_repro <- read_csv("gene_lists/upreg_in_repro.csv") #172
diff_exp <- rbind(upreg_in_nonrepro, upreg_in_repro) #334

# Exon gene list
diff_exon_gene_list <- read_csv("gene_lists/diff_exon_gene_list.csv") #59

### Number overlapping and get lists

merge1 <- merge(diff_meth_genes, diff_exp, by="geneID") #6
phyper(6, 276, 11030, 478, lower.tail = F) # p-val: 0.9506053
phyper(6, 276, 11030, 478, lower.tail = T) # p-val: 0.04939468 *** JUST LESS THAN BY CHANCE
colnames(merge1) <- "DMG_DEG_Genes"
merge2 <- merge(diff_meth_genes, upreg_in_nonrepro, by="geneID") #3
merge3 <- merge(diff_meth_genes, upreg_in_repro, by="geneID") #3
merge4 <- merge(diff_meth_genes, diff_exon_gene_list, by="geneID") #4
phyper(4, 59, 11030, 478, lower.tail = F) # p-val: 0.1096807
phyper(4, 59, 11030, 478, lower.tail = T) # p-val: 0.8903193
colnames(merge4) <- "DMG_DEE_Genes"

check <- merge(merge1, merge4, by="geneID") #0

merge5 <- merge(gene_hyper_repro_all, diff_exp, by="geneID") #1
merge6 <- merge(gene_hyper_repro_all, upreg_in_nonrepro, by="geneID") #0
merge7 <- merge(gene_hyper_repro_all, upreg_in_repro, by="geneID") #1
merge8 <- merge(gene_hyper_repro_all, diff_exon_gene_list, by="geneID") #0

merge9 <- merge(gene_hyper_sterile_all, diff_exp, by="geneID") #5
merge10 <- merge(gene_hyper_sterile_all, upreg_in_nonrepro, by="geneID") #3
merge11 <- merge(gene_hyper_sterile_all, upreg_in_repro, by="geneID") #2
merge12 <- merge(gene_hyper_sterile_all, diff_exon_gene_list, by="geneID") #0

merge13 <- merge(diff_exp, diff_exon_gene_list, by="geneID") #0

merge14 <- merge(gene_hyper_repro_all, gene_hyper_sterile_all, by="geneID") #51
phyper(51, 59, 11030, 478, lower.tail = F) # p-val: 1.632667e-64 *** AS EXPECTED
phyper(51, 59, 11030, 478, lower.tail = T) # p-val: 1
colnames(merge14) <- "DEG_DEE_Genes"

### --------------------------------------
## GO Lists
### --------------------------------------

# Hypergeometric test:
# Population size: 77573 (taken from number of overall GO terms defined in Eamonn's Fisher's Exact test script, same number found in Alun’s GO file made from the genome)
# Number of successes in population: methylation number of GO terms
# Sample size: number of exon/diff expressed GO terms
# Number of successes in sample: overlapping number 
 
diff_meth_GO_list <- read_csv("GO_lists/diff_meth_GO_list.csv") #532
diff_exp_GO_list <- read_csv("GO_lists/diff_exp_GO_list.csv") #101
diff_exon_GO_list <- read_csv("GO_lists/diff_exon_GO_list.csv") #9


merge15 <- merge(diff_meth_GO_list, diff_exp_GO_list, by="ids") #9
phyper(9, 101, 77573, 532, lower.tail = F) # p-val: 2.298605e-09 ** MORE THAN BY CHANCE
phyper(9, 101, 77573, 532, lower.tail = T) # p-val: 1
colnames(merge15) <- "DMG_DEG_GOs"

merge16 <- merge(diff_meth_GO_list, diff_exon_GO_list, by="ids") #2
phyper(2, 9, 77573, 532, lower.tail = F) # p-val: 2.611836e-05 ** AS ABOVE MORE THAN BY CHANCE
phyper(2, 9, 77573, 532, lower.tail = T) # p-val: 0.9999739
colnames(merge16) <- "DMG_DEE_GOs"

merge17 <- merge(diff_exp_GO_list, diff_exon_GO_list, by="ids") #8
phyper(8, 9, 77573, 101, lower.tail = F) # p-val: 7.445825e-27 ** AS ABOVE MORE THAN BY CHANCE
phyper(8, 9, 77573, 101, lower.tail = T) # p-val: 1
colnames(merge17) <- "DEG_DEE_GOs"


### --------------------------------------
## Write out signif file
### --------------------------------------

final <- cbind.fill(merge1, merge4, merge14, merge15, merge16, merge17, fill="")
write.csv(final, file="overlapping_genes_and_GOs.csv", row.names = F, quote=F)
