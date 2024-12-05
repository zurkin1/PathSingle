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
from sklearn.preprocessing import MinMaxScaler
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

        # Append all scores at once using list comprehension
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
Silhouette Score: 0.5818618817999396
Calinski-Harabasz Index: 17819.743450336544
Special accuracy: 0.3314285714285714
completeness score: 0.22306129171029385
homogeneity_score: 0.2461920972304762
adjusted_mutual_info_score: 0.2096547513727912
'''

progeny = decoupler.get_progeny(organism='human', top=2000)

def run_progeny():
    decoupler.run_mlm(mat=adata, net=progeny, source='source', target='target', weight='weight', verbose=False, use_raw=False)
    acts = decoupler.get_acts(adata, obsm_key='mlm_estimate')
    #Convert the pathway activity matrix to a DataFrame.
    return pd.DataFrame(acts.obsm['mlm_estimate'], index=adata.obs_names, columns=acts.var_names)
'''
Silhouette - mean: 0.5994470119476318, ci: (0.5837380086204144, 0.6151560152748493)
Calinski - mean: 1809.87665776202, ci: (1773.4737280450795, 1846.2795874789604)
Special - mean: 0.6271428571428571, ci: (0.6114787182419504, 0.6428069960437638)
Completeness - mean: 0.6207960242943827, ci: (0.6128990147889134, 0.6286930337998521)
Homogeneity - mean: 0.6358791241359868, ci: (0.6309284707239026, 0.640829777548071)
Adjusted - mean: 0.6159252478430198, ci: (0.6098366484002805, 0.6220138472857591)
'''

def run_aucell():
    decoupler.run_aucell(adata, reactome, source="geneset", target="genesymbol", use_raw=False, verbose=False)
    return adata.obsm["aucell_estimate"]
'''
Silhouette - mean: 0.5645602345466614, ci: (0.5418287845709704, 0.5872916845223524)
Calinski - mean: 945.4906096018539, ci: (915.5391413815059, 975.4420778222018)
Special - mean: 0.5946666666666667, ci: (0.5827601074933282, 0.6065732258400052)
Completeness - mean: 0.60093050891257, ci: (0.5916599459964447, 0.6102010718286954)
Homogeneity - mean: 0.6342885403046792, ci: (0.6296686799647475, 0.6389084006446109)
Adjusted - mean: 0.6045281693242305, ci: (0.5974978431803457, 0.6115584954681152)
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
    scaler = MinMaxScaler()
    output_activity = scaler.fit_transform(output_activity)
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
Results for pathsingle:
Silhouette - mean: 0.5358855452382626, ci: (0.5232767136639006, 0.5484943768126246)
Calinski - mean: 1090.5438033592732, ci: (1054.1400981262054, 1126.947508592341)
Special - mean: 0.5937619047619047, ci: (0.5835908538349613, 0.6039329556888481)
Completeness - mean: 0.6089631563405739, ci: (0.6041706907373023, 0.6137556219438454)
Homogeneity - mean: 0.6464026295938959, ci: (0.6402176559580726, 0.6525876032297193)
Adjusted - mean: 0.6148974658974095, ci: (0.6098015038476241, 0.6199934279471948)
'''

# Define list of method functions.
methods = [run_gsea, run_progeny, run_aucell, run_pathsingle]

# Loop through method functions.
for method_func in methods:
    method_name = method_func.__name__.replace('run_', '')  # Remove 'run_' prefix
    print(f"\nRunning {method_name.upper()}...")
    run_method(method_name, method_func)