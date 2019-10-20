import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import sys

# subtract time series for each station by its mean over 24hr
def subtractMean(usage):
    for i in usage.index:
        mean = usage.iloc[i].mean()
        usage.iloc[i] -= mean

def getClusters(usage, cluster_range, n_init, init=None):

    silhouette_avg = []

    for n_clusters in n_clusters_range:
        clustering = KMeans(n_clusters=n_clusters).fit(usage)

        labels = clustering.labels_
        traj = clustering.cluster_centers_
        silhouette_avg.append(silhouette_score(usage, labels))

    print silhouette_avg

# alternatively...generate cluster centroids by averaging
# then run getClusters using these centroids
def precomputeCentroids(usage, cluster_range, num_iterations):
    centroids = {}
    for n_clusters in n_clusters_range:
        centroids[n_clusters] = np.zeros((n_clusters, 24))

        for i in range(num_iterations):
            clustering = KMeans(n_clusters=n_clusters, init='random').fit(usage)
            traj = clustering.cluster_centers_ / num_iterations

            # correlate centroids from one iteration to the next
            if i ==0:
                centroids[n_clusters] += traj
            else:
                for row in range(n_clusters):
                    corrcoef_mat = np.corrcoef(centroids[n_clusters][row], traj)
                    # only look at first row for col index of max correlation
                    corr_max = corrcoef_mat[0][1:].max()
                    col_max = np.where(corrcoef_mat[0][1:] == corr_max)[0][0]
                    centroids[n_clusters][row] += traj[col_max]

    return centroids

def createDataFrame(nameDoc):

    fileOpen = open(nameDoc, 'r')
    matrixData = []

    line = fileOpen.readline()
    line2 = line

    while line:
        row = line.replace("\n", "").split(",")

        for i in range(len(row)):
            row[i] = float(row[i])
        matrixData.append(row)
        line = fileOpen.readline()

    #armamos el header
    header = []
    data = line2.replace("\n", "").split(",")

    for i in range(len(data)):
        element = 'col'+str(i+1)
        header.append(element)

    dataFrameValue = pd.DataFrame(matrixData, columns=header)
    return dataFrameValue

###############################################################################
if __name__ == '__main__':

    #generamos el dataFrame

    usage_all = createDataFrame(sys.argv[1])
    #stations = usage_all[usage_all.columns[0]]
    #usage_all.drop(usage_all.columns[0], axis=1, inplace=True)
    subtractMean(usage_all)


# we evalute the silhouette score for clustering on different parameters
# change n_clusters_range and the data fed into getClusters
    n_clusters_range = [2,3,4,5,6,7]
    #getClusters(usage_all, n_clusters_range, 100)

# precompute cluster centroids
    #centroids = precomputeCentroids(usage_all, n_clusters_range,10)
    getClusters(usage_all, n_clusters_range, 1, init=None)
