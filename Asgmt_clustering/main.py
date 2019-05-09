import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from mpl_toolkits.mplot3d import Axes3D # essential

# Preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import Normalizer
from sklearn.decomposition import PCA

# Clustering Module
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from sklearn.cluster import SpectralClustering
from sklearn import metrics
from scipy.stats import pearsonr

def pearson_affinity(M):
    result =  1 - np.array([[pearsonr(a,b)[0] for a in M] for b in M])
    print(result.shape)
    return result

# Hyper Parameters
IS_SAMPLE_CLUSTER = False

PCA_COMPONENTS = 0.9
N_CLUSTERS = 2
DBSCAN_EPS = 1
DBSCAN_MINS = 1
AGG_LINKAGE = "ward"#"average"
AGG_AFFINITY = "euclidean"#pearson_affinity
LOG_TRANSFORM = False
DATA_FILE = 'E:\\Programming\\Dataset\\GEO\\GSE33532\\hGSE33532.csv'
SIGNIF_GENES = 'E:\\Programming\\Dataset\\GEO\\GSE33532\\SignificantGenes.csv'

# Set random seed
random_state = 500

# clustering dict
his = {"adeno":1,"squamous":2,"mixed":3,"normal lung":0}
stage = {"1A":1,"1B":2,"2A":3,"2B":4,"na":0}
# Load data
filename = DATA_FILE
with open(filename,newline='') as csvfile:
    raw_data = list(csv.reader(csvfile))
    #sample_List = raw_data[0][1:-1]      # Sample-codes
    Labels = raw_data[1][1:-1]           # tumor-stage-labels
    #Labels = [his[x] for x in Labels]
    #Labels = [stage[x] for x in Labels]
    # remove the unannotated data
    clean_data = list(filter(lambda line:line[-1]!='',raw_data[2:]))
    clean_data = [x[1:] for x in clean_data]
    clean_data_frame = pd.DataFrame(clean_data)
    clean_data_frame.columns = raw_data[0][1:]
    clean_data_frame = clean_data_frame.drop_duplicates(['Gene Symbol'])
    clean_data = np.array(clean_data_frame).tolist()
    if not IS_SAMPLE_CLUSTER:
        Gene_list = [x[-1] for x in clean_data]
    Array = np.array([x[:-1] for x in clean_data],dtype='f')    # be careful to the dtype!
    print(np.array(Array).shape)
    csvfile.close()
    # We can finally get a data-matrix whose dimension is (14533*100)

# Load Significant-Genes labels from analysis of GEO2R
if not IS_SAMPLE_CLUSTER:
    filename = SIGNIF_GENES
    with open(filename,newline='') as csvfile:
        temp = list(csv.reader(csvfile))
        Signif_genes = [x[1] for x in temp[1:]]
        #print(Signif_genes[:10],len(Signif_genes))
        #print(Gene_list[:10],len(Gene_list))
        csvfile.close()

# Min-Msubfig scaling
if LOG_TRANSFORM:
    from sklearn.preprocessing import FunctionTransformer
    Array = FunctionTransformer(np.log).fit_transform(Array)
#Array = Normalizer().fit_transform(Array)
Array = StandardScaler().fit_transform(Array)
if IS_SAMPLE_CLUSTER:
    Array = Array.T


score_funcs = [
    metrics.adjusted_rand_score,
    metrics.v_measure_score,
    metrics.adjusted_mutual_info_score,
    metrics.mutual_info_score,
]


# PCA
pca = PCA(n_components=PCA_COMPONENTS,svd_solver='full')
Array_pca = pca.fit_transform(Array)
dedimension = Array_pca.shape
print('shape after PCA:',dedimension)

#K-Mean:
Kmeans_pre = KMeans(n_clusters=N_CLUSTERS,random_state=random_state).fit_predict(Array_pca)

# DBSCAN:
DBSCAN_pre = DBSCAN(eps=DBSCAN_EPS,min_samples=DBSCAN_MINS).fit_predict(Array_pca)
#DBSCAN_pre = DBSCAN().fit_predict(Array_pca)

# Hierarchical clustering
Agglom_pre = AgglomerativeClustering(n_clusters=N_CLUSTERS,
                                     affinity=AGG_AFFINITY,
                                     linkage=AGG_LINKAGE).fit_predict(Array_pca)
# Gaussian Mixture
Gmm_pre = GaussianMixture(n_components=N_CLUSTERS,random_state=random_state).fit_predict(Array_pca)

# Spectral Clustering
#Spect_pre = SpectralClustering(n_clusters=N_CLUSTERS).fit_predict(Array_pca)

#### Score for clustering for Samples ####
if IS_SAMPLE_CLUSTER:
    for score_func in score_funcs:
        print(score_func.__name__,'KMeans:',score_func(Kmeans_pre,Labels))
        print(score_func.__name__,'DBSCAN:',score_func(DBSCAN_pre,Labels))
        print(score_func.__name__,'Agglom:',score_func(Agglom_pre,Labels))
        print(score_func.__name__,'GMM:',score_func(Gmm_pre,Labels))
#       print(score_func.__name__,'Spectral:',score_func(Spect_pre,Labels))
else:
    SIGNIF_Labels = [0]*Array.shape[0]
    for i in range(Array.shape[0]):
        if Gene_list[i] in Signif_genes:
            SIGNIF_Labels[i] = 1
    for score_func in score_funcs:
        print(score_func.__name__,'KMeans:',score_func(Kmeans_pre,SIGNIF_Labels))
        print(score_func.__name__,'DBSCAN:',score_func(DBSCAN_pre,SIGNIF_Labels))
        print(score_func.__name__,'Agglom:',score_func(Agglom_pre,SIGNIF_Labels))
        print(score_func.__name__,'GMM:',score_func(Gmm_pre,SIGNIF_Labels))
##########################################

fig = plt.figure(1)

# KMeans result
# subfig = fig.add_subplot(2,2,1)
# subfig.scatter(Array_pca[:,0],Array_pca[:,1],c=Labels)
# subfig.set_title("True Labels")
subfig = fig.add_subplot(2,2,1)
subfig.scatter(Array_pca[:,0],Array_pca[:,1],c=Kmeans_pre)
subfig.set_title("KMeans: n_comp = %d"% N_CLUSTERS)
# subfig = fig.add_subplot(2,2,1,projection='3d')
# subfig.scatter(Array_pca[:,0], Array_pca[:,1],Array_pca[:,2],c=Kmeans_pre)
# subfig.set_title("KMeans: n_comp = %d"% N_CLUSTERS)

# DBSCAN result
subfig = fig.add_subplot(2,2,2)
subfig.scatter(Array_pca[:,0], Array_pca[:,1],c=DBSCAN_pre)
subfig.set_title("DBSCAN: eps= %d,min_samples= %d,n_comp = %d"% 
             (DBSCAN_EPS,DBSCAN_MINS,len(set(DBSCAN_pre))))
# subfig = fig.add_subplot(2,2,2,projection='3d')
# subfig.scatter(Array_pca[:,0], Array_pca[:,1],Array_pca[:,2],c=DBSCAN_pre)
# subfig.set_title("DBSCAN: eps= %d,min_samples= %d,n_comp = %d"% 
#              (DBSCAN_EPS,DBSCAN_MINS,len(set(DBSCAN_pre))))

# Hierarchical result
subfig = fig.add_subplot(2,2,3)
subfig.scatter(Array_pca[:,0], Array_pca[:,1],c=Agglom_pre)
subfig.set_title("Ward Linkage: n_comp = %d"% N_CLUSTERS)
# subfig = fig.add_subplot(2,2,3,projection='3d')
# subfig.scatter(Array_pca[:,0], Array_pca[:,1],Array_pca[:,2],c=Agglom_pre)
# subfig.set_title("Ward Linkage: n_comp = %d"% N_CLUSTERS)

# Gaussian Mixture result
subfig = fig.add_subplot(2,2,4)
subfig.scatter(Array_pca[:,0], Array_pca[:,1],c=Gmm_pre)
subfig.set_title("Gaussian Mixture: n_comp = %d"% N_CLUSTERS)

plt.show()

#### Clustering for Genes : Visualization ####
if not IS_SAMPLE_CLUSTER:
    fig = plt.figure(2)

    subfig1=fig.add_subplot(1,2,1)
    subfig2=fig.add_subplot(1,2,2)

    subfigs = [subfig1,subfig2]
    cluster1 = []
    cluster2 = []
    clusters = [cluster1,cluster2]
    for i in range(len(Kmeans_pre)):
        if len(cluster1) >=20 and len(cluster2) >=20:
            break
        if Kmeans_pre[i]==0 and len(cluster1) >=20:
            continue
        elif Kmeans_pre[i]==1 and len(cluster2) >=20:
            continue
        clusters[Kmeans_pre[i]].append(i)

    i = 0
    for cluster in clusters:
        for x in cluster:
            subfigs[i].plot(range(50),Array[x][-50:])
        i += 1

    subfig1.set_title("Normal Genes")
    subfig2.set_title("Differentially Expressed Genes")

    plt.show()
##############################################