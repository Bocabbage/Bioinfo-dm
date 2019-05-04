# 更新时间：2019/4/21(未完成)
#          2019/4/22(增加层次聚类；谱聚类小样本test；未完成)
#          2019/5/4(增加性能评估；未完成)

import csv
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D # essential

# Preprocessing
from sklearn.preprocessing import StandardScaler
#from sklearn.preprocessing import MinMsubfigScaler
from sklearn.decomposition import PCA

# Clustering Module
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from sklearn.mixture import GaussianMixture
#from sklearn.cluster import SpectralClustering
from sklearn import metrics

# Hyper Parameters
PCA_COMPONENTS_LIST = range(2,100)
#PCA_COMPONENTS = 0.9
N_CLUSTERS = 5
DBSCAN_EPS = 1
DBSCAN_MINS = 5

# Set random seed
random_state = 170

# Load data
filename = 'E:\\Programming\\Dataset\\GEO\\GSE33532\\GSE33532.csv'
with open(filename,newline='') as csvfile:
    rawdata = list(csv.reader(csvfile))
    probs_List = rawdata[0][1:]
    Gene_List = [x[0] for x in rawdata[1:]]
    # Will get a data matrix whose shape is (100,25906) and 5-classes labels
    Array = np.array([x[1:] for x in rawdata[2:]]).reshape(len(rawdata[0])-1,len(rawdata)-2)
    Labels = rawdata[1][1:]


# Min-Msubfig scaling
Array = StandardScaler().fit_transform(Array)

score_funcs = [
    metrics.adjusted_rand_score,
    metrics.v_measure_score,
    metrics.adjusted_mutual_info_score,
    metrics.mutual_info_score,
]


# Scores_Score lists
Scores_Kmeans = np.zeros([len(PCA_COMPONENTS_LIST),len(score_funcs)])
Scores_DBSCAN = np.zeros([len(PCA_COMPONENTS_LIST),len(score_funcs)])
Scores_Aggloom = np.zeros([len(PCA_COMPONENTS_LIST),len(score_funcs)])
Scores_Gmm = np.zeros([len(PCA_COMPONENTS_LIST),len(score_funcs)])

i=0
j=0

for PCA_COMPONENTS in PCA_COMPONENTS_LIST:
    j=0
    # Use PCA
    pca = PCA(n_components=PCA_COMPONENTS,svd_solver='full')
    Array_pca = pca.fit_transform(Array)
    #print(Array_pca.shape)

    # K-Mean:
    Kmeans_pre = KMeans(n_clusters=N_CLUSTERS,random_state=random_state).fit_predict(Array_pca)

    # DBSCAN:
    DBSCAN_pre = DBSCAN(eps=DBSCAN_EPS,min_samples=DBSCAN_MINS).fit_predict(Array_pca) 

    # Hierarchical clustering
    Agglom_pre = AgglomerativeClustering(n_clusters=N_CLUSTERS,
                                         affinity="euclidean",
                                         linkage="ward").fit_predict(Array_pca)
    # Gaussian Mixture
    Gmm_pre = GaussianMixture(n_components=N_CLUSTERS).fit_predict(Array_pca)

    for score_func in score_funcs:
        Scores_Kmeans[i,j] = score_func(Kmeans_pre,Labels)
        Scores_DBSCAN[i,j] = score_func(DBSCAN_pre,Labels)
        Scores_Aggloom[i,j] = score_func(Agglom_pre,Labels)
        Scores_Gmm[i,j] = score_func(Gmm_pre,Labels)
        j+=1
    i+=1

fig_score = plt.figure(1)
subfig = fig_score.add_subplot(2,2,1)
for i in range(0,len(score_funcs)):
    subfig.plot(PCA_COMPONENTS_LIST,Scores_Kmeans[:,i],label=score_funcs[i].__name__)
subfig.set_title("KMeans Scores-score")
subfig.legend()


subfig = fig_score.add_subplot(2,2,2)
for i in range(0,len(score_funcs)):
    subfig.plot(PCA_COMPONENTS_LIST,Scores_DBSCAN[:,i],label=score_funcs[i].__name__)
subfig.set_title("DBSCAN Scores-score")
subfig.legend()


subfig = fig_score.add_subplot(2,2,3)
for i in range(0,len(score_funcs)):
    subfig.plot(PCA_COMPONENTS_LIST,Scores_Aggloom[:,i],label=score_funcs[i].__name__)
subfig.set_title("Aggloom Scores-score")
subfig.legend()


subfig = fig_score.add_subplot(2,2,4)
for i in range(0,len(score_funcs)):
    subfig.plot(PCA_COMPONENTS_LIST,Scores_Gmm[:,i],label=score_funcs[i].__name__)
subfig.set_title("Gmm Scores-score")
subfig.legend()


plt.show()


# Spectral_clustering:
# SpecClus_pre = SpectralClustering(
#                 n_clusters=N_CLUSTERS,
#                 assign_labels='discretize',
#                 random_state=random_state).fit_predict(Array_pca[0:100,:])



# PCA
# pca = PCA(n_components=PCA_COMPONENTS,svd_solver='full')
# Array_pca = pca.fit_transform(Array)
# print(Array_pca.shape)

# #K-Mean:
# Kmeans_pre = KMeans(n_clusters=N_CLUSTERS,random_state=random_state).fit_predict(Array_pca)

# # DBSCAN:
# DBSCAN_pre = DBSCAN(eps=DBSCAN_EPS,min_samples=DBSCAN_MINS).fit_predict(Array_pca) 

# # Hierarchical clustering
# Agglom_pre = AgglomerativeClustering(n_clusters=N_CLUSTERS,
#                                      affinity="euclidean",
#                                      linkage="ward").fit_predict(Array_pca)
# # Gaussian Mixture
# Gmm_pre = GaussianMixture(n_components=N_CLUSTERS).fit_predict(Array_pca)

# print(metrics.adjusted_rand_score(Kmeans_pre,Labels))
# print(metrics.adjusted_rand_score(DBSCAN_pre,Labels))
# print(metrics.adjusted_rand_score(Agglom_pre,Labels))
# print(metrics.adjusted_rand_score(Gmm_pre,Labels))

# fig = plt.figure()

# # KMeans result
# subfig = fig.add_subplot(2,2,1)
# subfig.scatter(Array_pca[:,0],Array_pca[:,1],c=Kmeans_pre)
# subfig.set_title("KMeans: n_comp = %d"% N_CLUSTERS)
# # subfig = fig.add_subplot(2,2,1,projection='3d')
# # subfig.scatter(Array_pca[:,0], Array_pca[:,1],Array_pca[:,2],c=Kmeans_pre)
# # subfig.set_title("KMeans: n_comp = %d"% N_CLUSTERS)

# # DBSCAN result
# subfig = fig.add_subplot(2,2,2)
# subfig.scatter(Array_pca[:,0], Array_pca[:,1],c=DBSCAN_pre)
# subfig.set_title("DBSCAN: eps= %d,min_samples= %d,n_comp = %d"% 
#              (DBSCAN_EPS,DBSCAN_MINS,len(set(DBSCAN_pre))))
# # subfig = fig.add_subplot(2,2,2,projection='3d')
# # subfig.scatter(Array_pca[:,0], Array_pca[:,1],Array_pca[:,2],c=DBSCAN_pre)
# # subfig.set_title("DBSCAN: eps= %d,min_samples= %d,n_comp = %d"% 
# #              (DBSCAN_EPS,DBSCAN_MINS,len(set(DBSCAN_pre))))

# # Hierarchical result
# subfig = fig.add_subplot(2,2,3)
# subfig.scatter(Array_pca[:,0], Array_pca[:,1],c=Agglom_pre)
# subfig.set_title("Ward Linkage: n_comp = %d"% N_CLUSTERS)
# # subfig = fig.add_subplot(2,2,3,projection='3d')
# # subfig.scatter(Array_pca[:,0], Array_pca[:,1],Array_pca[:,2],c=Agglom_pre)
# # subfig.set_title("Ward Linkage: n_comp = %d"% N_CLUSTERS)

# # Gaussian Mixture result
# subfig = fig.add_subplot(2,2,4)
# subfig.scatter(Array_pca[:,0], Array_pca[:,1],c=Gmm_pre)
# subfig.set_title("Gaussian Mixture: n_comp = %d"% N_CLUSTERS)

# # Spectral-Clustering result
# # subfig = fig.add_subplot(2,2,4,projection='3d')
# # subfig.scatter(Array_pca[0:100,0], Array_pca[0:100,1],Array_pca[0:100,2],c=SpecClus_pre)
# # subfig.set_title("Spectral-Clustering: n_comp=%d(n_samples=100)"% N_CLUSTERS)

# plt.show()