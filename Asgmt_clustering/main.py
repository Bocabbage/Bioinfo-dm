# 更新时间：2019/4/21(未完成)
#          2019/4/22(增加层次聚类；谱聚类小样本test；未完成)

import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # essential
from sklearn.decomposition import PCA

from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import SpectralClustering

# Set random seed
random_state = 170

# Load data
filename = 'E:\\Programming\\Dataset\\GEO\\GSE33532\\GSE33532.csv'
with open(filename,newline='') as csvfile:
    rawdata = list(csv.reader(csvfile))
    GSM_List = rawdata[0][1:]
    Gene_List = [x[0] for x in rawdata[1:]]
    # Will get a data matrix which shape is (25906,100)
    Array = np.array([x[1:] for x in rawdata[1:]]).reshape(len(rawdata)-1,len(rawdata[0])-1)

# Use PCA
pca = PCA(n_components=3,svd_solver='full')
Array_pca = pca.fit_transform(Array)
#print(Array_pca.shape)

# K-Mean:
Kmeans_pre = KMeans(n_clusters=10,random_state=random_state).fit_predict(Array_pca)

# DBSCAN:
DBSCAN_pre = DBSCAN(eps=1,min_samples=1).fit_predict(Array_pca) 

# Hierarchical clustering
Agglom_pre = AgglomerativeClustering(n_clusters=10,
                                     affinity="euclidean",
                                     linkage="ward").fit_predict(Array_pca)

# Spectral_clustering:
SpecClus_pre = SpectralClustering(
                n_clusters=10,
                assign_labels='discretize',
                random_state=random_state).fit_predict(Array_pca[0:100,:])


fig = plt.figure()

# KMeans result
ax = fig.add_subplot(2,2,1,projection='3d')
ax.scatter(Array_pca[:,0], Array_pca[:,1],Array_pca[:,2],c=Kmeans_pre)
ax.set_title("KMeans: n_comp = 10")

# DBSCAN result
ax = fig.add_subplot(2,2,2,projection='3d')
ax.scatter(Array_pca[:,0], Array_pca[:,1],Array_pca[:,2],c=DBSCAN_pre)
ax.set_title("DBSCAN: eps=1,min_samples=1")

# Hierarchical result
ax = fig.add_subplot(2,2,3,projection='3d')
ax.scatter(Array_pca[:,0], Array_pca[:,1],Array_pca[:,2],c=Agglom_pre)
ax.set_title("Ward Linkage: n_comp = 10")

# Spectral-Clustering result
ax = fig.add_subplot(2,2,4,projection='3d')
ax.scatter(Array_pca[0:100,0], Array_pca[0:100,1],Array_pca[0:100,2],c=SpecClus_pre)
ax.set_title("Spectral-Clustering: n_comp=10(n_samples=100)")


plt.show()