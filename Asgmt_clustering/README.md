# Assignment-1: Do clustering for micro-array data

## 1.Target

* Apply different cluster methods for micro-array Data, including clustering the samples/genes.
* Compare the result of different methods

## 2.Data-Source & Instruments

### DataSource
* Database: [GEO](https://www.ncbi.nlm.nih.gov/geo/)
* Data: [GSE33532](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33532)
        (25906 probes, 100 samples)

### Instruments
* sklearn(python-package)
* matplotlib(python-package)
* Bioconductor(R-packages)


## 3.Preprocessing

### Data Clean

The raw-data need 2 steps of preproccessing.

##### Prob to ID/Gene Name
The raw data of the chip is a matrix of [probes\*samples], but what we actually analyze is [genes\*samples].
So we need to transform the probe-ids to gene-names. I use 'Bioconductor' package to download the GPL file for finishing this.
Also, It's important to add tags for each samples considering our clustering targets.I have try 2 kinds of tags:
Use Tumor-stage tags & Histology tags.

#### Ignore unannotated data and Deduplication
After transformation, you may find some probes without annotation and we need to delete them.
Also,It's a feature of micro-array data that one gene might be detected by more than one probe. So we need to deduplication.
I use python to do this.

### Standardized and Decomposition
* Standard-Scaler
* Principal Component Analysis(PCA).[para: n_components=0.9]


## 4.Algorithms

* [x] K-Means
* [x] DBSCAN
* [x] Hierarchical clustering
* [x] Gaussian Mixture

Meet problem that 'SpectralClustering' take about 99% of my RAM(12GB), also spending a long time to run but 
meets 'out of memory'. I just take 100 samples of the data and run them.

## 5.Results


