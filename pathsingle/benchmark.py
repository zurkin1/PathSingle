from __future__ import annotations
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler
from tqdm.notebook import tqdm
import os
import scanpy.external as sce
from metrics import *
from activity import *
from sklearn.preprocessing import Normalizer # Unit norm. StandardScaler # Normal distribution. MinMaxScaler # [0,1] range.
from scipy import stats
from itertools import chain, repeat
import urllib.request
import anndata as ad
from pathlib import Path
import scprep
import magic


os.environ["LOKY_MAX_CPU_COUNT"] = '4' #Prevents overloading system. Used by SKLearn, decoupler, MAGIC.

def run_method(method_name, method, adata, reactome):
    """Run a given method 30 times, calculate metrics and confidence intervals."""
    metrics = [[] for _ in range(6)] # Initialize 6 empty lists.
    metric_names = ['Silhouette', 'Calinski', 'Special', 'Completeness', 'Homogeneity', 'Adjusted']
    for i in tqdm(range(30)):
        pathway_activity_df = method(adata, reactome)
        #Perform KMeans clustering and plot UMAP.
        kmeans = cluster_with_kmeans(method_name, pathway_activity_df, adata, n_clusters=5) #Change to 10 for the PBMC benchmark.
        scores = calc_stats(pathway_activity_df, true_labels, kmeans.labels_, debug=True)
        print(i)
        for score, metric_list in zip(scores, metrics):
            metric_list.append(score)

    print(f'Results for {method_name}:')
    for name, metric in zip(metric_names, metrics):
        print(f"{name} - mean: {np.mean(metric)}, ci: {stats.t.interval(0.95, len(metric)-1, loc=np.mean(metric), scale=stats.sem(metric))}")

# Running GSEA. We will use the python package [`decoupler`](https://decoupler-py.readthedocs.io/en/latest/) <cite>`badia2022decoupler`</cite> to perform GSEA enrichment tests on our data.
# We use the normalized scores from sc.pp.normalize_total(adata) as a proxy for differential expression (DE) scores, which will significantly speed up the process since we don't have to 
# calculate DE scores for each cell individually.
def run_gsea(adata, reactome):
    #Prepare the result matrix for GSEA scores.
    num_cells = adata.shape[0]
    num_gene_sets = len(reactome['geneset'].unique())
    gsea_results_matrix = np.zeros((num_cells, num_gene_sets))

    #Loop through each cell to run GSEA.
    for cell_index in range(num_cells):
        #Get normalized expression values for the specific cell.
        cell_expr = adata.X[cell_index]
        #Create a DataFrame to hold DE scores.
        de_scores = pd.DataFrame(cell_expr, index=adata.var_names, columns=['scores'])
        #Run GSEA using decoupler.
        _, norm, _ = decoupler.run_gsea(de_scores.T, reactome, source="geneset", target="genesymbol")
        #Store the normalized enrichment scores (NES) in the result matrix.
        gsea_results_matrix[cell_index, :] = norm.iloc[:, 0].values
        #Save the result matrix for later use.
        #np.save('./data/gsea_results_matrix.npy', gsea_results_matrix)
        print(cell_index, end='\r')
    return gsea_results_matrix

def run_progeny(adata, reactome):
    progeny = decoupler.get_progeny(organism='human', top=2000)
    # Convert gene names to consistent format.
    adata.var_names = adata.var_names.str.upper()

    # Print overlap with PROGENy.
    progeny_genes = set(progeny['target'].unique())
    data_genes = set(adata.var_names)
    overlap = progeny_genes.intersection(data_genes)
    print("PROGENy genes:", len(progeny_genes), "Overlap genes:", len(overlap))

    decoupler.run_mlm(mat=adata, net=progeny, source='source', target='target', weight='weight', verbose=False, use_raw=False)
    acts = decoupler.get_acts(adata, obsm_key='mlm_estimate')
    #Convert the pathway activity matrix to a DataFrame.
    return pd.DataFrame(acts.obsm['mlm_estimate'], index=adata.obs_names, columns=acts.var_names)

def run_aucell(adata, reactome):
    decoupler.run_aucell(adata, reactome, source="geneset", target="genesymbol", use_raw=False, verbose=False)
    return adata.obsm["aucell_estimate"]

def run_pathsingle(adata, reactome):
    from sklearn.decomposition import PCA
    
    activity = sc.AnnData(adata.X, obs=adata.obs, var=adata.var)
    calc_activity(activity)
    output_activity = pd.read_csv('./data/output_activity.csv', index_col=0)

    #Scale the data.
    scaler = Normalizer() #For each cell (row), devide each activity by L2 norm of the row (square root of the sum of squares). 
                          #Each row will have length 1. print(np.sqrt(np.sum(X_normalized**2, axis=1)))  # [1. 1.]
    output_activity = scaler.fit_transform(output_activity)
    PCA = PCA(n_components=30, svd_solver='arpack')
    output_activity = PCA.fit_transform(output_activity)
    return output_activity


if __name__ == '__main__':
    # We first download the 68K PBMC data and follow the standard `scanpy` workflow for normalisation of read counts and subsetting on the highly variable genes. For the
    # T-cells data uncomment the following block.
    '''
    adata = sc.datasets.pbmc68k_reduced()
    adata.obs['labels'] = adata.obs.bulk_labels.map({'CD14+ Monocyte':0, 'Dendritic':1, 'CD56+ NK':2, 'CD4+/CD25 T Reg':3, 'CD19+ B':4, 'CD8+ Cytotoxic T':5, 'CD4+/CD45RO+ Memory':6, 'CD8+/CD45RA+ Naive Cytotoxic':7, 'CD4+/CD45RA+/CD25- Naive T':8, 'CD34+':9})
    true_labels = adata.obs.labels
    print(adata)
    '''

    num_splits = 5
    split_files = [f'./data/sc_training_split_{i+1}.h5ad' for i in range(num_splits)]
    splits = [sc.read_h5ad(file) for file in split_files]
    # Concatenate the splits back into a single AnnData object.
    adata = ad.concat(splits, join='outer', merge='same', label='batch', keys=[f'batch_{i+1}' for i in range(num_splits)])
    adata = sc.pp.subsample(adata, fraction=0.1, copy=True) #28697 cells × 15077 genes.
    print(adata)
    true_labels = adata.obs.state.map({'cycling':0, 'effector':1, 'other':2, 'progenitor':3, 'terminal exhausted':4})
    sc.pp.filter_genes(adata, min_cells=1)  # Remove unexpressed genes. Keep genes expressed in at least 1 cell.
    #sc.pp.normalize_total(adata)  # Library size normalization (works on adata.X).
    #sc.pp.sqrt(adata)             # Square root transformation (works on adata.X).
    #adata.raw = adata.copy()      # Copy adata.X plus other objects to adata.raw.
    adata.X = scprep.normalize.library_size_normalize(adata.X) #For each cell (row), divide each expression value by the sum of the row and multiply by the scaling factor (default 1e4).
    #adata.X = scprep.transform.sqrt(adata.X) #For each value x in the expression matrix take √x. Stabilizes variance and reduces outliers.

    # MAGIC imputation.
    print(adata.X.toarray()[:5,:5]) #adata.raw.to_adata().X.toarray()
    magic_op = magic.MAGIC()
    adata.X = magic_op.fit_transform(adata.X)
    adata.X = adata.X.astype(np.float16)
    #sce.pp.magic(adata, name_list='all_genes')
    print(adata.X[:5,:5]) #adata.raw.to_adata().X

    # Retrieving gene sets. Download and read the `gmt` file for the REACTOME pathways annotated in the C2 collection of MSigDB. 
    url = 'https://figshare.com/ndownloader/files/35233771'
    output_path = './data/c2.cp.reactome.v7.5.1.symbols.gmt'
    if not Path(output_path).is_file():
        print(f"Downloading {url} to {output_path}")
        urllib.request.urlretrieve(url, output_path)

    def gmt_to_decoupler(pth: Path) -> pd.DataFrame:
        """Parse a gmt file to a decoupler pathway dataframe."""
        pathways = {}
        with Path(pth).open("r") as f:
            for line in f:
                name, _, *genes = line.strip().split("\t")
                pathways[name] = genes

        return pd.DataFrame.from_records(
            chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
            columns=["geneset", "genesymbol"],
        )

    reactome = gmt_to_decoupler("./data/c2.cp.reactome.v7.5.1.symbols.gmt")
    # Define list of method functions.
    methods = [run_pathsingle] #run_gsea, run_progeny, run_aucell, 

    # Loop through method functions.
    for method_func in methods:
        method_name = method_func.__name__.replace('run_', '')  # Remove 'run_' prefix.
        print(f"\nRunning {method_name.upper()}...")
        run_method(method_name, method_func, adata, reactome)