#!/usr/bin/Rscript
#title      : deseq2_dea.r
#description: differential expression analysis by DESeq2
#author     : Zhuofan Zhang[modified the version of Huamei Li's]
#date       : 2019/5/19

# load packages
library(DESeq2)

# data preprocess
database <- read.table(file = "E:/temp/merge.counts", sep = "\t", header = T, row.names = 1)
database <- round(as.matrix(database))


# parse team condition and factor information
condition <- factor(c('treated', "untreated", "untreated", "untreated","treated","treated"), levels=c("untreated", "treated"))
coldata <- data.frame(row.names = colnames(database), condition)
dds <- DESeqDataSetFromMatrix(countData=database, colData=coldata, design=~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ] 

# set threshvalue and selected differ expression genes
dds <- DESeq(dds)
res  <- results(dds)
res <- res[order(res$padj), ]
diff_res <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))

# write results into specified file
file_name <- paste0("Blast", "_vs_", 'pHSC', "_DESeq2.csv")
dir.create('E:/temp/rna_seq analysis/diff', showWarnings = FALSE)
output_file <- file.path('E:/temp/rna_seq analysis/diff', file_name)
write.table(diff_res, file = output_file, row.names = TRUE, sep=",", quote=FALSE, col.names=NA)

