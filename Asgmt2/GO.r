#!/usr/bin/Rscript
#title      : GO.r
#description: GO annotation
#author     : Zhuofan Zhang
#date       : 2019/5/20

# load packages
library(clusterProfiler)

# load data
data <- read.csv("E:/temp/rna_seq analysis/diff/pHSC_vs_Blast_DESeq2.csv")

# Get genes list of gene-name/ensembl/entrezID
gene.list <- bitr(data$X,fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb="org.Hs.eg.db")

# enrichGO analysis

################## ALL RESULT ##################
# ego.all <- enrichGO(gene = gene.list$ENTREZID,
#                     OrgDb = org.Hs.eg.db,
#                     ont = "ALL",
#                     pAdjustMethod = "BH",
#                     pvalueCutoff = 1,
#                     qvalueCutoff = 1,
#                     keyType = 'ENTREZID',
#                     readable = TRUE)
################################################

################## MF RESULT ##################
ego.MF <- enrichGO(gene = gene.list$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   keyType = 'ENTREZID',
                   readable = FALSE)

# Visualization
# dotplot(ego.MF,title="enrichGO: p-value=0.05,q-value=0.2")
barplot(ego.MF,showCategory=20,title="enrichGO: p-value=0.05,q-value=0.2")
# 'topGO'/'Rgraphviz' package needed:
# plotGOgraph(ego.MF)
###############################################



