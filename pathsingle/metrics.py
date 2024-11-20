import numpy as np
from sklearn.metrics import silhouette_score, calinski_harabasz_score, completeness_score, homogeneity_score, adjusted_mutual_info_score
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import linear_sum_assignment
from scipy.stats import combine_pvalues
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.notebook import tqdm
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_selection import VarianceThreshold


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
    # Ensure y_true and y_pred are NumPy arrays of integers.
    y_true = np.asarray(y_true, dtype=np.int64)
    y_pred = np.asarray(y_pred, dtype=np.int64)

    assert y_pred.size == y_true.size

    # Determine the size of the confusion matrix and initialize it.
    D = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)

    # Fill the confusion matrix.
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    
    # Perform the assignment using the Hungarian algorithm.
    ind = linear_sum_assignment(w.max() - w)

    # Return accuracy.
    return sum(w[ind[0], ind[1]]) * 1.0 / y_pred.size


def print_stats(act_mat, true_labels, kmeans_labels):
    #Silhouette score.
    silhouette_avg = silhouette_score(act_mat, kmeans_labels, metric='euclidean')
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


def decoupler_feature_importance(scores, reactome, pvals=None):
    print(f"Shape of data: {scores.shape}")
    num_genesets_returned = scores.shape[1]

    # Extract the processed gene sets from the scores DataFrame.
    processed_gene_sets = scores.columns

    # Store the normalized enrichment scores (NES) in the result matrix.
    gsea_results_matrix = scores.T.values

    if pvals is not None:
        # Combine p-values across all cells using Fisher's method.
        combined_pvals = np.zeros(num_genesets_returned)
        for j in range(num_genesets_returned):
            # Replace zero values with a very small number (e.g., 1e-10).
            pvals_col = np.where(pvals.iloc[:, j] == 0, 1e-10, pvals.iloc[:, j])
            combined_pvals[j] = combine_pvalues(pvals_col)[1]

        # Create a DataFrame to store the combined p-values and sort by p-value.
        gene_sets = reactome['geneset'].unique()
        combined_pvals_df = pd.DataFrame({'geneset': processed_gene_sets, 'combined_pval': combined_pvals})
        combined_pvals_df = combined_pvals_df.sort_values('combined_pval')

        # Replace zero p-values with a very small number (e.g., 1e-10).
        combined_pvals_df['combined_pval'] = combined_pvals_df['combined_pval'].replace(0, 1e-10)
        # Replace p-values of 1 with a value slightly less than 1 (e.g., 1 - 1e-10).
        combined_pvals_df['combined_pval'] = combined_pvals_df['combined_pval'].replace(1, 1 - 1e-10)
        # Set a threshold for p-values.
        threshold = 1e-20
        combined_pvals_df['adjusted_pval'] = combined_pvals_df['combined_pval'].apply(lambda x: max(x, threshold))
        combined_pvals_df['pval_rank'] = combined_pvals_df['combined_pval'].rank()

        # Plot the distribution of p-values for the gene sets.
        #plt.figure(figsize=(10, 6))
        #plt.hist(combined_pvals_df['combined_pval'], bins=50, edgecolor='k')
        #plt.yscale('log')
        #plt.xlabel('Gene sets p-value')
        #plt.ylabel('Frequency')
        #plt.title('Distribution of Gene Sets p-values')
        #plt.show()
    else:
        combined_pvals_df = pd.DataFrame({'geneset': processed_gene_sets, 'combined_pval': 0, 'pval_rank': 0})

    # Use additional metric.
    combined_pvals_df['enrichment_score'] = scores.mean().values
    # Rank by combined p-value and enrichment score.
    combined_pvals_df['enrichment_rank'] = combined_pvals_df['enrichment_score'].rank(ascending=False)
    # Create a composite score.
    combined_pvals_df['composite_score'] = combined_pvals_df['pval_rank'] + combined_pvals_df['enrichment_rank']
    # Sort by composite score
    combined_pvals_df = combined_pvals_df.sort_values('composite_score')

    # Display the top 20 gene sets based on combined p-values.
    top_gene_sets = combined_pvals_df.head(20).copy()
    #top_gene_sets['-log10(combined_pval)'] = -np.log10(top_gene_sets['combined_pval'])
    #print(top_gene_sets.head())

    # Plot the top 20 gene sets with regards to combined ranks.
    plt.figure(figsize=(10, 8))
    sns.barplot(x='combined_pval', y='geneset', data=top_gene_sets)
    plt.xlabel('Gene sets p-value)')
    plt.ylabel('Gene Set')
    plt.title('Top 20 Gene Sets by p-value')
    plt.show()

    return top_gene_sets


def feature_importance(scores_df, true_labels):
        # Loop through the features, omitting one at a time and calculating the acc score.
        # Dictionary to store the accuracy results.
        acc_results = {}

        # Loop through each gene set.
        for pathway in tqdm(scores_df.columns):
            # Omit the current pathway.
            pathsingle_results_matrix_omitted = scores_df.drop(columns=[pathway]).values

            selector = VarianceThreshold(threshold=0.01)
            pathsingle_results_matrix_omitted = selector.fit_transform(pathsingle_results_matrix_omitted)

            #Scale the data.
            scaler = MinMaxScaler()
            pathsingle_results_matrix_omitted = scaler.fit_transform(pathsingle_results_matrix_omitted)
            
            # Perform KMeans clustering.
            num_clusters = len(np.unique(true_labels))
            kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(pathsingle_results_matrix_omitted)
            cluster_assignments = kmeans.labels_
            
            # Evaluate clustering results using acc metric.
            clustering_accuracy = acc(true_labels, cluster_assignments)
            
            # Store the accuracy result.
            acc_results[pathway] = clustering_accuracy

        # Sort the gene sets according to the acc result (the best gene set is the one that reduced the results the most).
        sorted_genesets = sorted(acc_results.items(), key=lambda x: x[1])

        # Display the sorted gene sets.
        print("Top 40 gene sets by clustering accuracy:")
        for pathway, accuracy in sorted_genesets[:40]:
            print(f'Pathway set: {pathway}')