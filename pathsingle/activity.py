import pandas as pd
import numpy as np
import scanpy as sc
from concurrent.futures import ProcessPoolExecutor


#Function to determine if an interaction type is inhibitory.
def is_inhibitory(interaction_type):
    inhibitory_keywords = ['inhibition', 'repression', 'dissociation', 'dephosphorylation', 'ubiquitination'] #'missing interaction'
    return any(keyword in interaction_type for keyword in inhibitory_keywords)

def gaussian_scaling(p, q, sigma=0.5):
    """Calculate the scaled activity of inputs (p) of an interaction by the outputs (q).
    1. Distance term: (p - q)**2 : Measures squared difference between p and q. Always positive. Larger when p and q are far apart.
    2. Scaling factor: exp(-(p-q)²/(2σ²)) : Returns 1.0 when p=q. Decreases exponentially as |p-q| increases. σ controls how quickly scaling drops off.
    3. Final value: p * scaling : When p=q: returns p (scaling=1). When p≠q: reduces p based on distance. Never increases above p.
    """
    scaling = np.exp(-(p - q)**2 / (2*sigma**2))
    return p * scaling

def process_pathway(args):
    """Calculate the activities of all pathways for a given sample."""
    pathway, interactions, gene_expression, gene_to_index = args
    pathway_activity = 0
    interactions_counter = 0
    interaction_activities = {}
    
    #Calculate activities for each interaction. interactions is a list of tuples: (['baiap2', 'wasf2', 'wasf3', 'wasf1'], 'activation', 'activation', ['cyp2b6', 'cyp2j2']).
    for interaction_idx, interaction in enumerate(interactions):
        interaction_activity = 0
        # Calculate input activity.
        for gene in interaction[0]:
            if gene in gene_to_index:
                interaction_activity += gene_expression[gene_to_index[gene]]
        
        # Calculate output activity.
        output_activity = 0
        for gene in interaction[2]:
            if gene in gene_to_index:
                output_activity += gene_expression[gene_to_index[gene]]
        output_activity = max(1e-10, output_activity)
        
        # Calculate the interaction activity using modified cross entropy function.
        #interaction_activity = -interaction_activity * (1 - np.log(output_activity))
        # Calculate interaction activity.
        interaction_activity = gaussian_scaling(interaction_activity, output_activity)
        #Once we finish calculating the interaction activity, we check if it is ihibitory.
        if is_inhibitory(interaction[1]):
            interaction_activity = -interaction_activity
            
        pathway_activity += interaction_activity
        interactions_counter += 1
        
        # Add interaction activity to the dictionary. item() converts tensor to float.
        interaction_activities[f'interaction_{pathway}_{interaction_idx}'] = interaction_activity
    
    return pathway, pathway_activity / interactions_counter, interaction_activities

def calc_activity(adata):
    """Calculate the activity of pathways based on gene expression data."""
    #Load gene expression data.
    #gene_expression_df = pd.read_csv('./data/sample_file.csv', index_col=0)
    #gene_expression_df.index = gene_expression_df.index.map(str.lower)

    #Convert the gene expression data to a PyTorch tensor. Transpose to match the desired orientation.
    gene_expression_tensor = adata.X #gene_expression_df.values (samples, genes)

    #Create a mapping of gene names to their indices in the gene_expression_tensor.
    gene_names = adata.var_names.str.lower() #gene_expression_df.index
    gene_to_index = {gene: i for i, gene in enumerate(gene_names)}

    #Load and parse pathway relations.
    pathway_relations = pd.read_csv('./data/pathway_relations.csv')
    pathway_relations['source'] = pathway_relations['source'].fillna('').astype(str).str.lower().str.split('*')
    pathway_relations['target'] = pathway_relations['target'].fillna('').astype(str).str.lower().str.split('*')

    #Parse pathway relations to create a simplified data structure.
    pathway_interactions = {}
    for _, row in pathway_relations.iterrows():
        pathway = row['pathway']
        sources = row['source']
        targets = row['target']
        inttype = row['interactiontype']
        if pathway not in pathway_interactions:
            pathway_interactions[pathway] = []
        pathway_interactions[pathway].append((sources,inttype,targets)) # For example: (['baiap2', 'wasf2', 'wasf3', 'wasf1'], 'activation', 'activation', ['cyp2b6', 'cyp2j2']).

    #Initialize a dictionary to store pathway activities.
    pathway_activities = {pathway: [] for pathway in pathway_interactions.keys()}

    # Initialize a list to store dictionaries of interaction activities per sample.
    interaction_dicts = []

    with ProcessPoolExecutor(max_workers=40) as executor:
        for sample_idx in range(gene_expression_tensor.shape[0]):
            sample_name = adata.obs_names[sample_idx]
            sample_data = gene_expression_tensor[sample_idx]
            interaction_dict = {'sample_name': sample_name}
            
            # Prepare arguments for parallel processing.
            args = [(pathway, interactions, sample_data, gene_to_index) 
                    for pathway, interactions in pathway_interactions.items()]
        
            # Process pathways in parallel.
            futures = [executor.submit(process_pathway, arg) for arg in args]
            results = [f.result() for f in futures]
            # Collect results.
            for pathway, pathway_activity, interaction_acts in results:
                pathway_activities[pathway].append(pathway_activity)
                interaction_dict.update(interaction_acts)
        
            interaction_dicts.append(interaction_dict)
            print(f"Processed sample {sample_idx+1}/{gene_expression_tensor.shape[0]}", end='\r')

    mean_activity_matrix = np.zeros((gene_expression_tensor.shape[0], len(pathway_interactions))) # (samples, pathways)

    for idx, (pathway_name, activities) in enumerate(pathway_activities.items()):
        if activities:
            mean_activity_matrix[:, idx] = activities
        else:
            print(f"No activities for pathway {pathway_name}")

    #Create the DataFrame for the activity matrix.
    activity_df = pd.DataFrame(mean_activity_matrix, index=adata.obs_names, columns=list(pathway_interactions.keys())).T

    #Save results to CSV.
    activity_df.T.to_csv('./data/output_activity.csv')
    
    # Convert the list of dictionaries to a DataFrame.
    interaction_activities = pd.DataFrame(interaction_dicts)
    # Set the sample name as the index.
    interaction_activities.set_index('sample_name', inplace=True)
    interaction_activities = interaction_activities.astype(np.float16)
    interaction_activities.to_csv('./data/output_interaction_activity.csv')