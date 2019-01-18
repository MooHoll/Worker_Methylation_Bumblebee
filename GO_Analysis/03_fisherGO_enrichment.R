### Eamonn's Script

## Identify enriched GO terms in a dataset compare to the GO terms of the whole
## transcriptome. 

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/Transcription/GO_analysis")

annot_GOs_finaltrans <- read.delim("annot_GOs_finaltrans.txt", header=F)
amelmethGO <- read.csv("repro_upreg_blast2go_go_propagation_20180116_1528.txt",sep='\t')#Replace this with your GO terms
amelmethGO <- amelmethGO[ which(amelmethGO$Level > 8),] #getting rid of the higher level GO terms (uniformative). Not sure why they are there

#Renaming the GO columns in each dataset so they can be identified
colnames(annot_GOs_finaltrans)[3] <- "transcriptome"
colnames(amelmethGO)[2] <- "meth"

#Getting rid of all the copies caused by snps being close together.
#amelmethGO<-unique(amelmethGO_temp[2:5])


#Counting the GO terms in each group
transcriptome <- summary(annot_GOs_finaltrans$transcriptome,maxsum = 10000)
dftrans<-data.frame(ids=names(transcriptome), nums=transcriptome)

meth <- summary(amelmethGO$meth,maxsum = 10000)
dfmeth<-data.frame(ids=names(meth), nums=meth)


#Combining the two lists of GO terms
contframe <- merge(dftrans, dfmeth, by = "ids", all=TRUE)
#contframe<-merge(dftrans,dfcolour, by ="ids", all=TRUE)
contframe[is.na(contframe)] <- 0
colnames(contframe)[2] <- "transcriptome"
colnames(contframe)[3] <- "methylation"

#Need to turn it into a matrix
contmatrix<-data.matrix(contframe)

#Need to input these values manually later (Lazy I know)
myvars <- c("methylation", "transcriptome")
totals <- contframe[myvars]
colSums(totals)


################################################################
## Fisher testR


FISHER<-function(contmatrix) {
         
## Prepare a two-dimensional contingency table
contingency.table <- data.frame(matrix(nrow=2, ncol=2))
rownames(contingency.table) <- c("GO term", "Non target terms")
colnames(contingency.table) <- c("Methylation", "Transcriptome")

## Assign the values one by one to make sure we put them in the right
## place (this is not necessary, we could enter the 4 values in a
## single instruction).

contingency.table["GO term", "Methylation"] <- contmatrix[['methylation']] ## Number of marked genes in the selection
contingency.table["GO term", "Transcriptome"] <- contmatrix[['transcriptome']] ## Number of non-marked genes in the selection
contingency.table["Non target terms", "Methylation"] <- 813-contmatrix[['methylation']] ## Number of marked genes outside of the selection
contingency.table["Non target terms", "Transcriptome"] <- 77573-contmatrix[['transcriptome']] ## Number of non-marked genes in the selection




## Run Fisher's exact test
ftest.result <- fisher.test(contingency.table, alternative="greater") 
return(ftest.result$p.value)
}



################################################################

#Adding p values and FDRs to contmatrix and creating a dataframe
p_value<-apply(contmatrix, 1,  FISHER)
p_value <- data.frame(matrix(unlist(p_value), nrow=6648, byrow=T)) #nrow manually taken from the P-value value in environment
results<-cbind(contframe,p_value)
colnames(results)[4] <- "p_value"
results$FDR<-p.adjust(results$p_value, method = "BH")

#Various subsettings and output
enriched <- results[results$ids %in% dfmeth$ids,] 
penriched<-enriched[which(enriched$p_value < 0.05),]
sigenrich<-enriched [ which(enriched$FDR < 0.05),]
write.csv(sigenrich, file = "repro_upreg_sigenrich.csv")
write.csv(penriched, file = "repro_upreg_penriched.csv")

##All you need for REVIGO is the list of GO terms and the p-values