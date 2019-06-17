#!/usr/bin/Rscript
#title      : GO.r
#description: GO annotation
#author     : Zhuofan Zhang
#date       : 2019/6/16

# load packages
library(clusterProfiler)

# load data
Up.data <- read.csv("E:/Programming/bioinfo_dm/Asgmt3/Integrative-result/Up_dea_genes.csv")
Down.data <- read.csv("E:/Programming/bioinfo_dm/Asgmt3/Integrative-result/Down_dea_genes.csv")
All.data <- read.csv("E:/Programming/bioinfo_dm/Asgmt3/Integrative-result/All_dea_genes.csv")
# Get genes list of gene-name/ensembl/entrezID
Up.gene.list <- bitr(Up.data$Up.dea.gene,fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb="org.Hs.eg.db")
Down.gene.list <- bitr(Down.data$Down.dea.gene,fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb="org.Hs.eg.db")
All.gene.list <- bitr(All.data$All.Dea.gene,fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb="org.Hs.eg.db")

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

################## GO RESULT ##################
Up.ego.ALL <- enrichGO(gene = Up.gene.list$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   keyType = 'ENTREZID',
                   readable = FALSE)

Down.ego.ALL <- enrichGO(gene = Down.gene.list$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   keyType = 'ENTREZID',
                   readable = FALSE)

All.ego.ALL <- enrichGO(gene = All.gene.list$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   keyType = 'ENTREZID',
                   readable = FALSE)
# Visualization
# dotplot(ego.MF,title="enrichGO: p-value=0.05,q-value=0.2")
barplot(Up.ego.ALL,showCategory=20,title="Upgenes-enrichGO: p-value=0.05,q-value=0.2")
barplot(Down.ego.ALL,showCategory=20,title="Downgenes-enrichGO: p-value=0.05,q-value=0.2")
barplot(All.ego.ALL,showCategory=20,title="DEAgenes-enrichGO: p-value=0.05,q-value=0.2")

# 'topGO'/'Rgraphviz' package needed:
# plotGOgraph(ego.MF)
###############################################

################## KEGG RESULT ##################
Up.kegg <- enrichKEGG(gene = Up.gene.list$ENTREZID,
                 organism="hsa",
                 keyType="kegg",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 pAdjustMethod = 'BH', 
                 minGSSize = 10,
                 maxGSSize = 500,
                 use_internal_data = FALSE)

Down.kegg <- enrichKEGG(gene = Down.gene.list$ENTREZID,
                 organism="hsa",
                 keyType="kegg",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 pAdjustMethod = 'BH', 
                 minGSSize = 10,
                 maxGSSize = 500,
                 use_internal_data = FALSE)

All.kegg <- enrichKEGG(gene = All.gene.list$ENTREZID,
                 organism="hsa",
                 keyType="kegg",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 pAdjustMethod = 'BH', 
                 minGSSize = 10,
                 maxGSSize = 500,
                 use_internal_data = FALSE)

barplot(Up.kegg,showCategory=20,title="Upgenes-KEGG: p-value=0.05,q-value=0.2")
barplot(Down.kegg,showCategory=20,title="Downgenes-KEGG: p-value=0.05,q-value=0.2")
barplot(All.kegg,showCategory=20,title="DEAgenes-KEGG: p-value=0.05,q-value=0.2")
#################################################




