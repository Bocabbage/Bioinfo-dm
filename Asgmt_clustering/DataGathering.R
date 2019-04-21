package.loader <- function()
{
    # Use the mirror in China
    local({BioC <- getOption("BioC_mirror");BioC["BioC_mirror"]
        <-"https://mirrors.ustc.edu.cn/bioc/";options(BioC_mirror=BioC)})

    # packages #
    # browseVignettes("GEOquery")
    # browseVignettes("SummarizedExperiment")
    library(GEOquery)
    library(SummarizedExperiment)

    print("The libraries have been loaded.")

}

setwd("E:/Programming/Dataset/GEO")
package.loader()

# Download data #
gse <- getGEO(GEO="GSE33532",destdir="./GSE33532")[[1]]
se <- as(gse,"SummarizedExperiment")

# Get the data matrix #
# metadata(se)
data <- assays(se)$exprs
# dim(data)
write.csv(x = data, file='./GSE33532/GSE33532.csv')

