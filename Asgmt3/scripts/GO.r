#!/usr/bin/Rscript
#title      : GO.r
#description: GO annotation
#author     : Zhuofan Zhang
#date       : 2019/6/15

# load packages
library(clusterProfiler)

# load data
Up.data <- read.csv("E:/Programming/bioinfo_dm/Asgmt3/Integrative-result/Up_dea_genes.csv")
Down.data <- read.csv("E:/Programming/bioinfo_dm/Asgmt3/Integrative-result/Down_dea_genes.csv")
# Get genes list of gene-name/ensembl/entrezID
Up.gene.list <- bitr(Up.data$Up.dea.gene,fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb="org.Hs.eg.db")
Down.gene.list <- bitr(Down.data$Up.dea.gene,fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb="org.Hs.eg.db")
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
Up.ego.MF <- enrichGO(gene = Up.gene.list$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   keyType = 'ENTREZID',
                   readable = FALSE)

Down.ego.MF <- enrichGO(gene = Down.gene.list$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   keyType = 'ENTREZID',
                   readable = FALSE)
# Visualization
# dotplot(ego.MF,title="enrichGO: p-value=0.05,q-value=0.2")
barplot(Up.ego.MF,showCategory=20,title="Upgenes-enrichGO: p-value=0.05,q-value=0.2")
barplot(Down.ego.MF,showCategory=20,title="Downgenes-enrichGO: p-value=0.05,q-value=0.2")
# 'topGO'/'Rgraphviz' package needed:
# plotGOgraph(ego.MF)
###############################################



