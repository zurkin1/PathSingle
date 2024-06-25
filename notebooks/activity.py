# %%
#Step 1: Prepare Pathway Matrix Before Runtime.
import pandas as pd
import numpy as np
import json
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm


#Load pathway relations data.
pathway_relations_df = pd.read_csv('./data/pathway_relations.csv')

#Ensure 'source' column is properly cleaned and converted to lists of strings.
pathway_relations_df['source'] = pathway_relations_df['source'].fillna('').astype(str).str.lower().str.split('*')

#Create a binary matrix for pathway relations with adjustments for interaction types.
pathway_genes = pathway_relations_df['source']
all_genes = sorted(set(sum(pathway_genes.tolist(), [])))

gene_index = {gene: i for i, gene in enumerate(all_genes)}
pathway_index = {pathway: i for i, pathway in enumerate(pathway_relations_df['pathway'].unique())}

pathway_matrix = np.zeros((len(all_genes), len(pathway_index)))


#Function to determine if an interaction type is inhibitory.
def is_inhibitory(interaction_type):
    inhibitory_keywords = ['inhibition', 'repression', 'dissociation', 'missing interaction', 'dephosphorylation', 'ubiquitination']
    return any(keyword in interaction_type for keyword in inhibitory_keywords)


for _, row in pathway_relations_df.iterrows():
    pathway = row['pathway']
    interaction_type = row['interactiontype']
    adjustment = -1 if is_inhibitory(interaction_type) else 1
    for gene in row['source']:
        if gene in gene_index:
            pathway_matrix[gene_index[gene], pathway_index[pathway]] = adjustment

#Save the pathway matrix and indices
np.save('pathway_matrix.npy', pathway_matrix)
with open('gene_index.json', 'w') as f:
    json.dump(gene_index, f)
with open('pathway_index.json', 'w') as f:
    json.dump(pathway_index, f)

# %%
#Step 2: Load Pathway Matrix and Gene Expression Data at Runtime.
#Load pathway matrix and indices.
def calc_activity_from_adata(adata):
    #On the index we have cell names (adata.obs_names) and the columns are gene names (adata.var_names). We transpose this dataframe before calculating activity.
    df = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names).T

    #Load pathway matrix and indices.
    pathway_matrix = np.load('./data/pathway_matrix.npy')
    with open('./data/gene_index.json', 'r') as f:
        gene_index = json.load(f)
    with open('./data/pathway_index.json', 'r') as f:
        pathway_index = json.load(f)

    #Load gene expression data.
    gene_expression_df = df #pd.read_csv('./data/sample_file.csv')
    gene_expression_df['GeneSymbol'] = gene_expression_df.index.str.lower() #['GeneSymbol'].str.lower()

    #Take intersection of genes present in both datasets.
    intersecting_genes = list(set(gene_expression_df['GeneSymbol']) & set(gene_index.keys()))
    intersecting_indices = [gene_index[gene] for gene in intersecting_genes]

    #Create a reduced pathway matrix.
    reduced_pathway_matrix = pathway_matrix[intersecting_indices, :]
    gene_expression_matrix = gene_expression_df.set_index('GeneSymbol').loc[intersecting_genes].values

    #Convert to PyTorch tensors.
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    reduced_pathway_tensor = torch.tensor(reduced_pathway_matrix, dtype=torch.float32).to(device)
    gene_expression_tensor = torch.tensor(gene_expression_matrix, dtype=torch.float32).to(device)

    #Step 3: Parallelize Multiplication Using Data Loaders.
    #Define a Dataset.
    class PathwayDataset(torch.utils.data.Dataset):
        def __init__(self, pathway_tensor):
            self.pathway_tensor = pathway_tensor

        def __len__(self):
            return self.pathway_tensor.shape[1]

        def __getitem__(self, idx):
            return self.pathway_tensor[:, idx]

    #Create DataLoader.
    batch_size = 1000  #Adjust based on your memory constraints.
    pathway_dataset = PathwayDataset(reduced_pathway_tensor)
    pathway_loader = DataLoader(pathway_dataset, batch_size=batch_size, shuffle=False)

    #Initialize a list to hold results.
    results = []

    #Process the multiplication in batches using DataLoader.
    for pathway_batch in tqdm(pathway_loader):
        pathway_batch = pathway_batch.to(device)
        result_batch = torch.matmul(pathway_batch, gene_expression_tensor) #Multiplying a batch of 1000 rows from the pathway_matrix numpy array map we created, with the gene expression matrix.
        results.append(result_batch)

    #Concatenate the results.
    activity_matrix = torch.cat(results, dim=0)

    #Convert back to CPU and NumPy if needed.
    activity_matrix = activity_matrix.cpu().numpy()

    #Save results to CSV.
    activity_df = pd.DataFrame(activity_matrix, index=list(pathway_index.keys()), columns=gene_expression_df.columns[1:])
    activity_df.to_csv('./data/output_activity.csv')