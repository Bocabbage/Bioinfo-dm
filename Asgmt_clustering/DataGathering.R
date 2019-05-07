package.loader <- function()
{
    # Use the mirror in China
    local({
        BioC <- getOption("BioC_mirror");
        BioC["BioC_mirror"]<-"https://mirrors.ustc.edu.cn/bioc/";
        options(BioC_mirror=BioC)
    })

    # packages #
    # browseVignettes("GEOquery")
    library(GEOquery)
    #library(SummarizedExperiment)

    print("The libraries have been loaded.")
}

setwd("E:/Programming/Dataset/GEO")
package.loader()

# Download data #
gse <- getGEO(GEO="GSE33532",destdir="./GSE33532",getGPL=T)[[1]]
gpl <- getGEO(filename = "./GSE33532/GPL570.soft")
# Get the data matrix #
# metadata(se)
exprset <- exprs(gse)
exprset <- as.data.frame(exprset)   # essential step!

pdata <- pData(gse)
group_list = as.character(pdata[, 39])

{
    # Tag the labels to the data #

    # IA.exprs <- exprset[,grep('1A',group_list)]
    # IB.exprs <- exprset[,grep('1B',group_list)]
    # IIA.exprs <- exprset[,grep('2A',group_list)]
    # IIB.exprs <- exprset[,grep('2B',group_list)]
    # na.exprs <- exprset[,grep('normal lung',as.character(pdata[,39]))]

    # IA.exprs <- rbind('1A',IA.exprs)
    # IB.exprs <- rbind('1B',IB.exprs)
    # IIA.exprs <- rbind('2A',IIA.exprs)
    # IIB.exprs <- rbind('2B',IIB.exprs)
    # na.exprs <- rbind('na',na.exprs)

    A.exprs <- exprset[,grep('adeno',group_list)]
    B.exprs <- exprset[,grep('squamous',group_list)]
    C.exprs <- exprset[,grep('mixed',group_list)]
    na.exprs <- exprset[,grep('normal lung',group_list)]

    A.exprs <- rbind('adeno',A.exprs)
    B.exprs <- rbind('squamous',B.exprs)
    C.exprs <- rbind('mixed',C.exprs)
    na.exprs <- rbind('normal lung',na.exprs)

    # re-bind the data together
    # exprset <- cbind(IA.exprs,IB.exprs,IIA.exprs,IIB.exprs,na.exprs)
    exprset <- cbind(A.exprs,B.exprs,C.exprs,na.exprs)
    exprset$ID <- rownames(exprset)
    gene.entrez.id <- Table(gpl)[,c('ID','Gene Symbol')]
    exprdata <- merge(x=exprset,y=gene.entrez.id,by='ID',all.x=T)
    exprdata$ID <- NULL
}

# write.csv(x=exprdata,file='./GSE33532/GSE33532.csv')
write.csv(x=exprdata,file='./GSE33532/hGSE33532.csv')