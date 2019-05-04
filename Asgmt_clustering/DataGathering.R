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

# Get the data matrix #
# metadata(se)
exprset <- exprs(gse)
pdata <- pData(gse)
group_list = as.character(pdata[, 41])

{
    # Tag the labels to the data #

    IA.exprs <- exprset[,grep('1A',group_list)]
    IB.exprs <- exprset[,grep('1B',group_list)]
    IIA.exprs <- exprset[,grep('2A',group_list)]
    IIB.exprs <- exprset[,grep('2B',group_list)]
    na.exprs <- exprset[,grep('normal lung',as.character(pdata[,39]))]

    IA.exprs <- rbind('1A',IA.exprs)
    IB.exprs <- rbind('1B',IB.exprs)
    IIA.exprs <- rbind('2A',IIA.exprs)
    IIB.exprs <- rbind('2B',IIB.exprs)
    na.exprs <- rbind('na',na.exprs)

    # re-bind the data together
    exprdata <- cbind(IA.exprs,IB.exprs,IIA.exprs,IIB.exprs,na.exprs)
}

write.csv(x=exprdata,file='./GSE33532/GSE33532.csv')
