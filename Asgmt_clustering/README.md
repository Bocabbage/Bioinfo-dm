# Assignment-1: Do clustering for micro-array data

## 1.Target

* Apply different cluster methods for micro-array Data
* Compare the result of different methods

## 2.Data-Source

* Database: [GEO](https://www.ncbi.nlm.nih.gov/geo/)
* Data: [GSE33532](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33532)
        (25906 genes, 100 samples)

## 3.Algorithms

* [x]K-Means
* [x]DBSCAN
* [x]Hierarchical clustering
* [ ]Spectral clustering

Meet problem that 'SpectralClustering' take about 99% of my RAM(12GB), also spends a long time to run but 
meet 'out of memory'. I just take 100 samples of the data and run them.

![image](figs/4_clustering_04-22.png)

And here is the time-cost performance of different algorithms:

![image](figs/performance.png)

