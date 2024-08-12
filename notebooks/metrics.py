import numpy as np
from sklearn.metrics import silhouette_score, calinski_harabasz_score, completeness_score, homogeneity_score, adjusted_mutual_info_score
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import linear_sum_assignment


#General metric functions used in the benchmarks.

#Dunn index calculation function.
def dunn_index(X, labels):
    #Compute pairwise distances.
    distances = squareform(pdist(X))
    
    #Find unique cluster labels.
    unique_labels = np.unique(labels)
    
    #Calculate inter-cluster distances (minimum distance between points in different clusters).
    inter_cluster_dists = np.inf
    for i in range(len(unique_labels)):
        for j in range(i + 1, len(unique_labels)):
            cluster_i = X[labels == unique_labels[i]]
            cluster_j = X[labels == unique_labels[j]]
            inter_cluster_dists = min(inter_cluster_dists, np.min(distances[np.ix_(labels == unique_labels[i], labels == unique_labels[j])]))
    
    #Calculate intra-cluster distances (maximum distance within a cluster).
    intra_cluster_dists = 0
    for label in unique_labels:
        cluster = X[labels == label]
        intra_cluster_dists = max(intra_cluster_dists, np.max(pdist(cluster)) if len(cluster) > 1 else 0)
    
    #Calculate Dunn index.
    return inter_cluster_dists / intra_cluster_dists if intra_cluster_dists != 0 else np.inf


#Special accuracy function for clustering.
def acc(y_true, y_pred):
    """
    Calculate clustering accuracy. This acc function calculates clustering accuracy, which is a measure of how well a clustering algorithm assigns
    data points to clusters compared to the ground truth labels. It does so by finding the best matching between the predicted cluster assignments
    (y_pred) and the true labels (y_true), and then calculating the accuracy based on this matching. Here's how the function works:
    1. It first ensures that both y_true and y_pred are integer arrays.
    2. It asserts that the sizes of y_true and y_pred are equal.
    3. It initializes a confusion matrix w where each entry w[i, j] denotes the number of data points that belong to cluster i according to y_pred
       and to cluster j according to y_true.
    4. It computes the optimal one-to-one matching between the predicted clusters and the true clusters using the Hungarian algorithm, implemented
       in linear_assignment from scikit-learn.
    5. It calculates the accuracy by summing the values of the confusion matrix corresponding to the optimal matching and dividing by the total 
       number of data points.
    6. This accuracy measure is specifically tailored for clustering algorithms, where there is no inherent ordering of cluster labels. It 
       accounts for the fact that clusters may be assigned arbitrary labels by the algorithm, and it finds the best matching between these labels and the true labels to compute accuracy.

    Traditional accuracy metrics assume a one-to-one correspondence between predicted and true labels, which is not the case in clustering. 
    Instead, the acc function finds the best matching between clusters and true labels, providing a more appropriate measure of clustering 
    performance. We need to find the best matching between the predicted clusters and the true clusters and then evaluate the performance based on
    this matching. This is what the acc function does by finding the optimal matching using the Hungarian algorithm and computing the accuracy
    based on this matching.
    Arguments
        y: true labels, numpy.array with shape `(n_samples,)`
        y_pred: predicted labels, numpy.array with shape `(n_samples,)`
    Return
        accuracy, in [0,1]
    """
    y_true = y_true.astype(np.int64)
    assert y_pred.size == y_true.size
    D = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    ind = linear_sum_assignment(w.max() - w)
    return sum(w[ind[0], ind[1]]) * 1.0 / y_pred.size


def print_stats(act_mat, true_labels, kmeans_labels):
    #Silhouette score.
    silhouette_avg = silhouette_score(act_mat, kmeans_labels)
    print(f"Silhouette Score: {silhouette_avg}")

    #Calinski-Harabasz index.
    calinski_harabasz = calinski_harabasz_score(act_mat, kmeans_labels)
    print(f"Calinski-Harabasz Index: {calinski_harabasz}")

    #Special accuracy function.
    print(f"Special accuracy: {acc(true_labels, kmeans_labels)}")

    #Completeness score.
    print(f'completeness score: {completeness_score(true_labels, kmeans_labels)}')

    #Homogeneity_score.
    print(f"homogeneity_score: {homogeneity_score(true_labels, kmeans_labels)}")

    #Adjusted_mutual_info_scor
    print(f"adjusted_mutual_info_score: {adjusted_mutual_info_score(true_labels, kmeans_labels)}")