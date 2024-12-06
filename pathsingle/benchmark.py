from __future__ import annotations
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler
from tqdm.notebook import tqdm
import session_info
import os
import warnings
import scanpy.external as sce
from metrics import *
import scprep
from sklearn.cluster import KMeans
import umap.umap_ as umap
import matplotlib.pyplot as plt
import magic
from activity import *
from pathlib import Path
from sklearn.preprocessing import StandardScaler  # Normal distribution.
from sklearn.preprocessing import MinMaxScaler    # [0,1] range.
from sklearn.preprocessing import Normalizer      # Unit norm.
from sklearn.feature_selection import VarianceThreshold
from scipy import stats
from itertools import chain, repeat
import urllib.request


os.environ["LOKY_MAX_CPU_COUNT"] = '4'
sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))
#Filtering warnings from current version of matplotlib.
warnings.filterwarnings("ignore", message=".*Parameters 'cmap' will be ignored.*", category=UserWarning)
warnings.filterwarnings("ignore", message="Tight layout not applied.*", category=UserWarning)

# We first download the 68K PBMC data and follow the standard `scanpy` workflow for normalisation of read counts and subsetting on the highly variable genes.
adata = sc.datasets.pbmc68k_reduced()
adata.obs['labels'] = adata.obs.bulk_labels.map({'CD14+ Monocyte':0, 'Dendritic':1, 'CD56+ NK':2, 'CD4+/CD25 T Reg':3, 'CD19+ B':4, 'CD8+ Cytotoxic T':5, 'CD4+/CD45RO+ Memory':6, 'CD8+/CD45RA+ Naive Cytotoxic':7, 'CD4+/CD45RA+/CD25- Naive T':8, 'CD34+':9})
true_labels = adata.obs.labels
print(adata)

# Run Magic.
print(adata.raw.to_adata().X.toarray()[:5,:5])
sce.pp.magic(adata, name_list='all_genes')
print(adata.raw.to_adata().X[:5,:5])

# Retrieving gene sets. Download and read the `gmt` file for the REACTOME pathways annotated in the C2 collection of MSigDB. 
os.makedirs('./data', exist_ok=True)
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
'''
print(reactome.head())
geneset	genesymbol
0	REACTOME_INTERLEUKIN_6_SIGNALING	JAK2
1	REACTOME_INTERLEUKIN_6_SIGNALING	TYK2
2	REACTOME_INTERLEUKIN_6_SIGNALING	CBL
3	REACTOME_INTERLEUKIN_6_SIGNALING	STAT1
4	REACTOME_INTERLEUKIN_6_SIGNALING	IL6ST
'''

def run_method(method_name, method):
    """Run a given method 30 times, calculate metrics and confidence intervals."""
    metrics = [[] for _ in range(6)] # Initialize 6 empty lists.
    metric_names = ['Silhouette', 'Calinski', 'Special', 'Completeness', 'Homogeneity', 'Adjusted']
    for _ in tqdm(range(30)):
        pathway_activity_df = method()
        #Perform KMeans clustering and plot UMAP.
        kmeans = cluster_with_kmeans(method_name, pathway_activity_df, adata, n_clusters=10)
        scores = calc_stats(pathway_activity_df, true_labels, kmeans.labels_)
        for score, metric_list in zip(scores, metrics):
            metric_list.append(score)

    print(f'Results for {method_name}:')
    for name, metric in zip(metric_names, metrics):
        print(f"{name} - mean: {np.mean(metric)}, ci: {stats.t.interval(0.95, len(metric)-1, loc=np.mean(metric), scale=stats.sem(metric))}")

# Running GSEA. We will use the python package [`decoupler`](https://decoupler-py.readthedocs.io/en/latest/) <cite>`badia2022decoupler`</cite> to perform GSEA enrichment tests on our data.
# We use the normalized scores from sc.pp.normalize_total(adata) as a proxy for differential expression (DE) scores, which will significantly speed up the process since we don't have to 
# calculate DE scores for each cell individually.
def run_gsea():
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
'''
Results for gsea:
Silhouette - mean: 0.5561046677924425, ci: (0.5534683894232597, 0.5587409461616253)
Calinski - mean: 16294.160463012933, ci: (15962.582829852261, 16625.738096173605)
Special - mean: 0.31719047619047624, ci: (0.31164243899418, 0.3227385133867725)
Completeness - mean: 0.2179805089345317, ci: (0.21609472267907712, 0.21986629518998627)
Homogeneity - mean: 0.24682039680701184, ci: (0.24492730757050962, 0.24871348604351406)
Adjusted - mean: 0.20717204963411157, ci: (0.20538114825921425, 0.20896295100900888)
'''

progeny = decoupler.get_progeny(organism='human', top=2000)

def run_progeny():
    decoupler.run_mlm(mat=adata, net=progeny, source='source', target='target', weight='weight', verbose=False, use_raw=False)
    acts = decoupler.get_acts(adata, obsm_key='mlm_estimate')
    #Convert the pathway activity matrix to a DataFrame.
    return pd.DataFrame(acts.obsm['mlm_estimate'], index=adata.obs_names, columns=acts.var_names)
'''
Results for progeny:
Silhouette - mean: 0.6144946217536926, ci: (0.59571121144838, 0.6332780320590052)
Calinski - mean: 1802.2639849762188, ci: (1755.1765120853065, 1849.351457867131)
Special - mean: 0.6490476190476191, ci: (0.63785434370834, 0.6602408943868981)
Completeness - mean: 0.6347861014556874, ci: (0.6295494352718396, 0.6400227676395351)
Homogeneity - mean: 0.6431875507404569, ci: (0.6408570554777759, 0.6455180460031378)
Adjusted - mean: 0.6269649005116761, ci: (0.6235347482357637, 0.6303950527875885)
'''

def run_aucell():
    decoupler.run_aucell(adata, reactome, source="geneset", target="genesymbol", use_raw=False, verbose=False)
    return adata.obsm["aucell_estimate"]
'''
Results for aucell:
Silhouette - mean: 0.5813433527946472, ci: (0.5586026801822613, 0.6040840254070331)
Calinski - mean: 985.4724865675855, ci: (945.7024049948561, 1025.242568140315)
Special - mean: 0.6157142857142858, ci: (0.6029025957025977, 0.6285259757259738)
Completeness - mean: 0.6249309793986906, ci: (0.616305298690219, 0.6335566601071622)
Homogeneity - mean: 0.6483159001639436, ci: (0.6441170250837653, 0.6525147752441218)
Adjusted - mean: 0.6243653270044098, ci: (0.618072118848946, 0.6306585351598737)
'''

def run_pathsingle():
    from sklearn.decomposition import PCA
    
    activity_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    #magic_op = magic.MAGIC()
    #activity_df = magic_op.fit_transform(activity_df)
    activity = sc.AnnData(activity_df)

    calc_activity(activity)
    output_activity = pd.read_csv('./data/output_interaction_activity.csv', index_col=0)

    #selector = VarianceThreshold(threshold=0.01)
    #output_activity = selector.fit_transform(output_activity)

    #Scale the data.
    #scaler = MinMaxScaler()
    #output_activity = scaler.fit_transform(output_activity)
    PCA = PCA(n_components=50, svd_solver='arpack')
    output_activity = PCA.fit_transform(output_activity)
    return output_activity
'''
Silhouette Score: 0.5577373354652382
Calinski-Harabasz Index: 1055.0147987294065
Special accuracy: 0.6214285714285714
completeness score: 0.6103618744991666
homogeneity_score: 0.6504527943574239
adjusted_mutual_info_score: 0.6176909258981184
'''

# Define list of method functions.
methods = [run_pathsingle] #run_gsea, run_progeny, run_aucell, 

# Loop through method functions.
for method_func in methods:
    method_name = method_func.__name__.replace('run_', '')  # Remove 'run_' prefix
    print(f"\nRunning {method_name.upper()}...")
    run_method(method_name, method_func)