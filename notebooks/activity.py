# %%
#Step 1: Prepare Pathway Matrix Before Runtime.
import pandas as pd
import numpy as np
import json
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm


MAX_NUMBER_OF_INTERACTIONS = 200


def create_pathway_matrix():
    #Load pathway relations data.
    pathway_relations_df = pd.read_csv('./data/pathway_relations.csv')

    #Ensure 'source' column is properly cleaned and converted to lists of strings.
    pathway_relations_df['source'] = pathway_relations_df['source'].fillna('').astype(str).str.lower().str.split('*')

    #Create a binary matrix for pathway relations with adjustments for interaction types.
    pathway_genes = pathway_relations_df['source']
    all_genes = sorted(set(sum(pathway_genes.tolist(), [])))

    gene_index = {gene: i for i, gene in enumerate(all_genes)}
    pathway_index = {pathway: i for i, pathway in enumerate(pathway_relations_df['pathway'].unique())}

    pathway_matrix = np.zeros((len(all_genes), len(pathway_index), MAX_NUMBER_OF_INTERACTIONS))
    inhibition_flags = np.zeros(len(pathway_index))
    current_interaction_counter = {}


    #Function to determine if an interaction type is inhibitory.
    def is_inhibitory(interaction_type):
        inhibitory_keywords = ['inhibition', 'repression', 'dissociation', 'missing interaction', 'dephosphorylation', 'ubiquitination']
        return any(keyword in interaction_type for keyword in inhibitory_keywords)


    for _, row in pathway_relations_df.iterrows():
        pathway = row['pathway']
        interaction_type = row['interactiontype']
        
        #Initialize current_interaction_counter dictionary.
        interaction_counter = current_interaction_counter.get(pathway, 0)
        if interaction_counter == 0:
            current_interaction_counter[pathway] = 0
        
        #If row['source'] have more than one gene we will use a new pathway_matrix, otherwise we will use the one in index 0 to save memory.
        len_row = len(row['source'])
        for gene in row['source']:
            if gene in gene_index:
                if len_row == 1:
                    pathway_matrix[gene_index[gene], pathway_index[pathway], 0] = 1
                else:
                    pathway_matrix[gene_index[gene], pathway_index[pathway], current_interaction_counter[pathway]] = 1
        
        if len_row > 1:
            current_interaction_counter[pathway] += 1
        
        if is_inhibitory(interaction_type):
            inhibition_flags[pathway_index[pathway]] = 1

    #Save the pathway matrix and indices.
    np.save('./data/pathway_matrix.npy', pathway_matrix)
    np.save('./data/inhibition_flags.npy', inhibition_flags)
    with open('./data/gene_index.json', 'w') as f:
        json.dump(gene_index, f)
    with open('./data/pathway_index.json', 'w') as f:
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
    
    #Create a mapping of gene names to their indices in the gene_expression_tensor.
    gene_to_index = {gene: i for i, gene in enumerate(gene_expression_df['GeneSymbol'])}

    #Take intersection of genes present in both datasets.
    intersecting_genes = list(set(gene_expression_df['GeneSymbol']) & set(gene_index.keys()))
    intersecting_indices = [gene_index[gene] for gene in intersecting_genes]

    #Create a reduced pathway matrix and a reduced gene expression matrix.
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

# %%
import numpy as np
import pandas as pd
import json
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm
import gc
import scanpy as sc


epsilon = 1e-10
#Load gene expression data.
#gene_expression_df = pd.read_csv('./data/sample_file.csv', index_col=0)
#gene_expression_df.index = gene_expression_df.index.map(str.lower)
adata = sc.read_h5ad('./data/sc_trainingmagic.h5ad')

#Convert the gene expression data to a PyTorch tensor.
gene_expression_tensor = torch.tensor(adata.X.T, dtype=torch.float32)  # Transpose to match the desired orientation

#Create a mapping of gene names to their indices in the gene_expression_tensor.
gene_names = adata.var_names.str.lower()
gene_to_index = {gene: i for i, gene in enumerate(gene_names)}

#Create a mapping of gene names to their indices in the gene_expression_tensor.
#gene_to_index = {gene: i for i, gene in enumerate(gene_expression_df.index)}

#Load pathway relations.
pathway_relations = pd.read_csv('./data/pathway_relations.csv')

#Ensure 'source' column is properly cleaned and converted to lists of strings.
pathway_relations['source'] = pathway_relations['source'].fillna('').astype(str).str.lower().str.split('*')

#Parse pathway relations to create a simplified data structure.
pathway_interactions = {}


#Function to determine if an interaction type is inhibitory.
def is_inhibitory(interaction_type):
    inhibitory_keywords = ['inhibition', 'repression', 'dissociation', 'missing interaction', 'dephosphorylation', 'ubiquitination']
    return any(keyword in interaction_type for keyword in inhibitory_keywords)


for _, row in pathway_relations.iterrows():
    pathway = row['pathway']
    sources = row['source']
    inttype = row['interactiontype']
    
    if pathway not in pathway_interactions:
        pathway_interactions[pathway] = []
    
    pathway_interactions[pathway].append((sources,inttype))

#Convert gene expression data to PyTorch tensor.
#gene_expression_tensor = torch.tensor(gene_expression_df.values, dtype=torch.float32)

#Define a Dataset.
class GeneExpressionDataset(torch.utils.data.Dataset):
    def __init__(self, gene_expression_tensor):
        self.gene_expression_tensor = gene_expression_tensor

    def __len__(self):
        return self.gene_expression_tensor.shape[1]

    def __getitem__(self, idx):
        return self.gene_expression_tensor[:, idx]

#Create DataLoader.
batch_size = 5
gene_expression_dataset = GeneExpressionDataset(gene_expression_tensor)
gene_expression_loader = DataLoader(gene_expression_dataset, batch_size=batch_size, shuffle=False)

#Initialize a dictionary to store pathway activities.
pathway_activities = {pathway: [] for pathway in pathway_interactions.keys()}

#Process the multiplication in batches using DataLoader.
for gene_expression_batch in tqdm(gene_expression_loader):
    gene_expression_batch = gene_expression_batch.to('cpu')
    
    for sample in range(0,gene_expression_batch.shape[0]):

        #Calculate activities for each pathway.
        for pathway, interactions in pathway_interactions.items():
                pathway_activity = 0
                interactions_counter = 0
                for interaction in interactions:
                    interaction_activity = 1
                    for gene in interaction[0]:
                        if gene in gene_to_index:
                            interaction_activity *= gene_expression_batch[sample, gene_to_index[gene]]
                    
                    if is_inhibitory(interaction[1]):
                        interaction_activity = 1 - interaction_activity
                    pathway_activity += interaction_activity
                    interactions_counter += 1
                
                #Calculate the mean activity for each pathway and for all cells.                
                pathway_activity /= interactions_counter
            
                pathway_activities[pathway].append(pathway_activity) #.cpu().numpy())

mean_activity_matrix = np.zeros((len(adata.obs_names), len(pathway_interactions)))

for idx, (pathway_name, activities) in enumerate(pathway_activities.items()):
    if activities:
        mean_activity_matrix[:, idx] = activities
    else:
        print(f"No activities for pathway {pathway_name}")

#Create the DataFrame for the activity matrix.
activity_df = pd.DataFrame(mean_activity_matrix, index=adata.obs_names, columns=list(pathway_interactions.keys())).T

#Save results to CSV.
activity_df.to_csv('./data/output_activity.csv')