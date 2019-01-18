### Making a count table for differential expression analysis
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/1. Bombus epigenetics/differential_expression_reproductive_status/counts_to_Bter1.0")

library(readr)

# samples (NOTE: count files must be in working directory)
file.list = list.files(("~/Dropbox/PhD Folder/Virtual Lab Book/1. Bombus epigenetics/differential_expression_reproductive_status/counts_to_Bter1.0"),pattern="*.counts")
samples = as.character(sapply(file.list,function(s) gsub("_.counts","", s)))
samples
length(samples) # 18

# sample info
info = read.csv("count_table_metadata.csv")
head(info,18)

# put two together and rbind everything together in a data frame
# takes about 30 mins to run on my mac
out = do.call(rbind, lapply(1:length(samples),
       function(i) {
cnts = read_delim(file.list[[i]], "\t", escape_double = F, col_names = F, trim_ws = T) # counts
colnames(cnts) = c("geneID","count")
sampinfo = info[info$sample==samples[[i]],] # sample info
suppressWarnings(data.frame(sampinfo, cnts))
} ))
head(out)
nrow(out) 
out[out$geneID=="LOC100645068",]$sample # test, woop there are 18 levels

write.csv(out, "counts_diff_exp.csv")
