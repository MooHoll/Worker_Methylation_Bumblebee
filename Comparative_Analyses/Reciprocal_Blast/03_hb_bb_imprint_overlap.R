##########################################
# Looking for common diff meth genes HB and BB
##########################################

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Methylation/honeybee_overlap")

library(readr)
library(tidyverse)
library(rowr)
library(VennDiagram)
library(stringr)

# Dataframe with BB LOC and HB GB corresponding gene codes
orthologs <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Honeybee/reciprocal_blast_BB_HB.csv")
head(orthologs)
# Filter so only working with 6953 genes but they are reliable matches
orthologs <- subset(orthologs, orthologs$is_unique_hit == "yes")
orthologs <- subset(orthologs, orthologs$match_type == "matches_in_both_blasts")

# Gene lists from Lyko
honeybee_gene_lists <- read_csv("lyko_diff_meth_genes.csv")
head(honeybee_gene_lists)
length(na.omit(unique(honeybee_gene_lists$GLEAN))) #549 genes


# Relevant gene lists from our BB analysis
bumblebee_gene_lists <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Methylation/NEW_diff_meth_methylkit/cpg/list_DMGs_MSCfilter.csv", 
                              col_names = FALSE)
head(bumblebee_gene_lists)
length(na.omit(unique(bumblebee_gene_lists$X1))) #478 diff meth genes


# ------------------------------------------
# Looking for overlap
# ------------------------------------------
overlap <- cbind.fill(honeybee_gene_lists$GLEAN, bumblebee_gene_lists$X1,
                         fill = NA)
overlap <- overlap[rowSums(is.na(overlap)) != ncol(overlap), ]
colnames(overlap) <- c("hb_gene", "bb_gene")
head(overlap)

# Convert the bumblbee gene IDs to the matching HB ortholog so can compare
merge1 <- as.data.frame(overlap$bb_gene)
colnames(merge1) <- "BB"
merge1 <- merge(merge1, orthologs, by = "BB") #417/478 with orthologs
merge1 <- as.data.frame(merge1$HB)
colnames(merge1) <- "BB_orthologes"

# Convert the honeybee gene IDs to the matching BB ortholog 
#merge2 <- as.data.frame(overlap$hb_gene)
#colnames(merge2) <- "HB"
#merge2 <- merge(merge2, orthologs, by = "HB") #X/561 with orthologs
#merge2 <- as.data.frame(merge2$BB)
#colnames(merge2) <- "HB_orthologes"

# Put both gene sets together (honeybee diff meth and hb ortholog of bb diff meth)
bumblbee_lists_HB <- cbind.fill(merge1$BB_orthologes, overlap$hb_gene,
                                     fill = NA)
colnames(bumblbee_lists_HB) <- c("BB_orthologs", "HB_diff_meth")
bumblbee_lists_HB <- as.data.frame(apply(bumblbee_lists_HB,
                                              2,function(x)gsub('\\s+', '',x)))

# For venn diagram need each overlapping dataset as a vector
BB_orthologes <- na.omit(bumblbee_lists_HB$BB_orthologs)
HB_diff_meth <- na.omit(bumblbee_lists_HB$HB_diff_meth)


venn.plot <- venn.diagram(
  x = list(BB_orthologes,HB_diff_meth),
  filename = NULL,
  category = c("A.mel Orthologs of B.ter Diff Meth", 
               "A.mel Diff Meth"),
  fill = c("dodgerblue", "seagreen3"),
  cat.col = c("dodgerblue", "seagreen3"),
  cat.cex = 1,
  margin = 0.05,
  cex = 3
)
grid.newpage()
grid.draw(venn.plot)


### No genes overlapping 