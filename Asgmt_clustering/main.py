import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from mpl_toolkits.mplot3d import Axes3D # essential

# Preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA

# Clustering Module
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from sklearn import metrics

# Hyper Parameters
PCA_COMPONENTS = 0.9
N_CLUSTERS = 5
DBSCAN_EPS = 0.1
DBSCAN_MINS = 1
AGG_LINKAGE = "ward"
AGG_AFFINITY = "euclidean"

# Set random seed
random_state = 170

# clustering dict
his = {"adeno":0,"squamous":1,"mixed":2,"normal lung":3}

# Load data
filename = 'E:\\Programming\\Dataset\\GEO\\GSE33532\\hGSE33532.csv'
with open(filename,newline='') as csvfile:
    raw_data = list(csv.reader(csvfile))
    sample_List = raw_data[0][1:-1]      # Sample-codes
    Labels = raw_data[1][1:-1]           # tumor-stage-labels
    # remove the unannotated data
    clean_data = list(filter(lambda line:line[-1]!='' and line[-1][0].isalpha(),raw_data[2:]))
    clean_data = [x[1:] for x in clean_data]
    clean_data_frame = pd.DataFrame(clean_data)
    clean_data_frame.columns = raw_data[0][1:]
    clean_data_frame = clean_data_frame.drop_duplicates(['Gene Symbol'])
    clean_data = np.array(clean_data_frame).tolist()
    #clean_data_frame.to_csv("./test.csv",sep=',')
    Array = np.array([x[:-1] for x in clean_data])
    # We can finally get a data-matrix whose dimension is (14533*100)
    print(np.array(Array).shape)


    

# Min-Msubfig scaling
#Array = MinMaxScaler().fit_transform(Array)
Array = StandardScaler().fit_transform(Array)
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
dedimension = Array_pca.shape[1]
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

#### Score for clustering for Samples ####

for score_func in score_funcs:
    print(score_func.__name__,'KMeans:',score_func(Kmeans_pre,Labels))
    print(score_func.__name__,'DBSCAN:',score_func(DBSCAN_pre,Labels))
    print(score_func.__name__,'Agglom:',score_func(Agglom_pre,Labels))
    print(score_func.__name__,'GMM:',score_func(Gmm_pre,Labels))

##########################################

fig = plt.figure(1)

# KMeans result
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
# fig = plt.figure(2)

# subfig1=fig.add_subplot(3,3,1)
# subfig2=fig.add_subplot(3,3,2)
# subfig3=fig.add_subplot(3,3,3)
# subfig4=fig.add_subplot(3,3,4)
# subfig5=fig.add_subplot(3,3,5)
# subfig6=fig.add_subplot(3,3,6)
# subfig7=fig.add_subplot(3,3,7)
# subfig8=fig.add_subplot(3,3,8)
# subfig9=fig.add_subplot(3,3,9)


# subfigs = [subfig1,subfig2,subfig3,subfig4,
#            subfig5,subfig6,subfig7,subfig8,subfig9]

# i = 0
# for label in Agglom_pre[:100]:
#     subfigs[label].plot(range(0,100),Array[i])
#     i+=1

# plt.show()
##############################################