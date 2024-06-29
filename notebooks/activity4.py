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


epsilon = 1e-10
#Load pathway matrix and indices.
pathway_matrix = np.load('./data/pathway_matrix.npy')
with open('./data/gene_index.json', 'r') as f:
    gene_index = json.load(f)
with open('./data/pathway_index.json', 'r') as f:
    pathway_index = json.load(f)

#Load gene expression data.
gene_expression_df = pd.read_csv('./data/sample_file.csv', index_col=0)
gene_expression_df.index = gene_expression_df.index.map(str.lower)

#Create a mapping of gene names to their indices in the gene_expression_tensor.
gene_to_index = {gene: i for i, gene in enumerate(gene_expression_df.index)}

#Take intersection of genes present in both datasets.
intersecting_genes = list(set(gene_expression_df.index) & set(gene_index.keys()))
intersecting_indices = [gene_index[gene] for gene in intersecting_genes] #These are the indices in the pathway matrix.

#Create a reduced pathway matrix and a reduced gene expression matrix.
reduced_pathway_matrix = pathway_matrix[intersecting_indices, :, :] # genes x pathways x interactions
reduced_gene_expression_matrix = gene_expression_df.loc[intersecting_genes].values # genes x cells

#Convert to PyTorch tensors.
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
reduced_pathway_tensor = torch.tensor(reduced_pathway_matrix, dtype=torch.float32).to(device)
gene_expression_tensor = torch.tensor(reduced_gene_expression_matrix, dtype=torch.float32).to(device)

#Step 3: Parallelize Multiplication Using Data Loader. Define a Dataset.
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

#Initialize a list to hold results and a dictionary to aggregate by pathway.
results = []
pathway_activities = {pathway: [] for pathway in list(pathway_index.keys())}

#Expand pathway_tensor to match the batch dimension of gene_expression_batch.
expanded_pathway_tensor = reduced_pathway_tensor.unsqueeze(1).T  # shape: [num_genes, num_pathways, max_interactions, 1]
print(f'expanded_pathway_tensor shape:', {expanded_pathway_tensor.shape}) #[200, 314, 1, 4865] interactions, pathways, batch, genes.

#Process the multiplication in batches using DataLoader.
for gene_expression_batch in tqdm(gene_expression_loader):
    gene_expression_batch = gene_expression_batch.unsqueeze(0) 
    #print(f'gene_expression_batch shape:', {gene_expression_batch.shape}) #[1, 10, 4865]

    #Perform the element-wise (dot) multiplication, with broadcasting. Gene_expression_batch has to be unsqueezed to match the pathway tensor.
    gene_expression_batch = gene_expression_batch * expanded_pathway_tensor  # shape: [interactions, pathways, batch_size, num_genes]
    #print(f'gene_expression_batch shape: {gene_expression_batch.shape}') #[200, 314, 4, 4865]

    #Multiply along the gene dimension to get the interaction activity for each pathway 
    # (in the previous example, removing the 4865 dimension which represents genes that participates in interactions in pathways).
    #First replace zeros with ones in place in order to neutralize zero value genes in the product.
    #Convert to log space, replacing zeros with epsilon.
    gc.collect()
    gene_expression_batch = torch.where(gene_expression_batch == 0, torch.tensor(epsilon, dtype=gene_expression_batch.dtype, device=gene_expression_batch.device), gene_expression_batch).log()

    #Sum in log space along the desired dimension (dim=3) (instead of prod(dim=3)).
    gene_expression_batch = gene_expression_batch.sum(dim=3)

    #Convert back from log space to get the product.
    gene_expression_batch = gene_expression_batch.exp()
    gc.collect()
    #print(f'gene_expression_batch shape: {gene_expression_batch.shape}') #[200, 314, 4865]
    
    #Aggregate the result by pathway, since we have multiple lines with the same pathway name, representing different interactions.
    for i in range(gene_expression_batch.shape[1]):
        pathway_name = list(pathway_index.keys())[i]
        pathway_activities[pathway_name].extend(gene_expression_batch[:, i, :])

#Calculate the mean activity for each pathway and for all cells.
mean_activity_matrix = np.zeros((len(gene_expression_df.columns), len(pathway_index)))
print(f'mean_activity_matrix shape: {mean_activity_matrix.shape}') #(20, 314) [num_cells, num_pathways]

for idx, (pathway_name, activities) in enumerate(pathway_activities.items()):
    if activities: #len(activities) == <number of cells>
        mean_activity = np.mean(activities, axis=0) #Calculate the mean along the correct axis.
        mean_activity_matrix[:, idx] = mean_activity
    else:
        print(f"No activities for pathway {pathway_name}")

#Create the DataFrame for the activity matrix.
activity_df = pd.DataFrame(mean_activity_matrix, index=gene_expression_df.columns, columns=list(pathway_index.keys())).T

#Save results to CSV.
activity_df.to_csv('./data/output_activity.csv')
# %%
create_pathway_matrix()
# %%
