import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm
import scanpy as sc


def calc_activity(adata):
    #Load gene expression data.
    #gene_expression_df = pd.read_csv('./data/sample_file.csv', index_col=0)
    #gene_expression_df.index = gene_expression_df.index.map(str.lower)
    
    #adata = sc.read_h5ad('./data/sc_trainingmagic.h5ad')
    #adata = sc.pp.subsample(adata, fraction=0.2, copy=True)
    print(adata)

    #Convert the gene expression data to a PyTorch tensor. Transpose to match the desired orientation.
    gene_expression_tensor = torch.tensor(adata.X.T, dtype=torch.float32) #gene_expression_df.values

    #Create a mapping of gene names to their indices in the gene_expression_tensor.
    gene_names = adata.var_names.str.lower() #gene_expression_df.index
    gene_to_index = {gene: i for i, gene in enumerate(gene_names)}

    #Load pathway relations.
    pathway_relations = pd.read_csv('./data/pathway_relations.csv')

    #Ensure 'source' column is properly cleaned and converted to lists of strings.
    pathway_relations['source'] = pathway_relations['source'].fillna('').astype(str).str.lower().str.split('*')

    #Parse pathway relations to create a simplified data structure.
    pathway_interactions = {}


    #Function to determine if an interaction type is inhibitory.
    def is_inhibitory(interaction_type):
        inhibitory_keywords = ['inhibition', 'repression', 'dissociation', 'dephosphorylation', 'ubiquitination'] #'missing interaction'
        return any(keyword in interaction_type for keyword in inhibitory_keywords)


    for _, row in pathway_relations.iterrows():
        pathway = row['pathway']
        sources = row['source']
        inttype = row['interactiontype']
        
        if pathway not in pathway_interactions:
            pathway_interactions[pathway] = []
        
        pathway_interactions[pathway].append((sources,inttype))

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

    # Initialize a list to store dictionaries of interaction activities.
    interaction_dicts = []

    #Process the multiplication in batches using DataLoader.
    for batch_idx, gene_expression_batch in enumerate(tqdm(gene_expression_loader)):
        gene_expression_batch = gene_expression_batch.to('cpu')
        
        for sample_idx in range(0,gene_expression_batch.shape[0]):
            sample_name = adata.obs_names[batch_idx * batch_size + sample_idx]
            interaction_dict = {'sample_name': sample_name}
            
            #Calculate activities for each pathway. interactions is a list of tuples: (['baiap2', 'wasf2', 'wasf3', 'wasf1'], 'activation').
            for pathway, interactions in pathway_interactions.items():
                pathway_activity = 0
                interactions_counter = 0
                for interaction_idx, interaction in enumerate(interactions):
                    interaction_activity = 0
                    for gene in interaction[0]:
                        if gene in gene_to_index:
                            interaction_activity += gene_expression_batch[sample_idx, gene_to_index[gene]]
                    #Once we finish calculating the interaction activity, we check if it is ihibitory.
                    if is_inhibitory(interaction[1]):
                        interaction_activity = -interaction_activity
                    pathway_activity += interaction_activity
                    interactions_counter += 1

                # Add interaction activity to the dictionary.
                interaction_dict[f'interaction_{pathway}_{interaction_idx}'] = interaction_activity
                
                #Once we finish with the pathway, calculate the mean activity for the pathway.               
                pathway_activity /= interactions_counter
                pathway_activities[pathway].append(pathway_activity) #.cpu().numpy())
            
            interaction_dicts.append(interaction_dict)

    # Convert the list of dictionaries to a DataFrame.
    interaction_activities = pd.DataFrame(interaction_dicts)
    # Set the sample name as the index.
    interaction_activities.set_index('sample_name', inplace=True)

    mean_activity_matrix = np.zeros((len(adata.obs_names), len(pathway_interactions)))

    for idx, (pathway_name, activities) in enumerate(pathway_activities.items()):
        if activities:
            mean_activity_matrix[:, idx] = activities
        else:
            print(f"No activities for pathway {pathway_name}")

    #Create the DataFrame for the activity matrix.
    activity_df = pd.DataFrame(mean_activity_matrix, index=adata.obs_names, columns=list(pathway_interactions.keys())).T

    #Save results to CSV.
    activity_df.T.to_csv('./data/output_activity.csv')
    interaction_activities.to_csv('./data/output_interaction_activity.csv')


if __name__ == '__main__':
    cdata = sc.read('./data/emt_magic.csv', delimiter=',', cache=True) #emt_magic
    calc_activity(cdata)