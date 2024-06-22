import pandas as pd
from tqdm import tqdm


#Load the necessary data.
gene_expression_file = './data/sample_file.csv'
pathway_relations_file = './data/pathway_relations.csv'
output_file = './data/output_activity.csv'

#Load gene expression data.
gene_expression_df = pd.read_csv(gene_expression_file, index_col=0)
gene_expression_df.index = gene_expression_df.index.map(str.lower)


#Split the 'source' and 'target' columns by '*' and lowercase the gene names.
def split_and_lower(x):
    if isinstance(x, str):
        return x.lower().split('*')
    return []


#Load pathway relations data.
pathway_relations_df = pd.read_csv(pathway_relations_file)
pathway_relations_df['source'] = pathway_relations_df['source'].apply(split_and_lower)
pathway_relations_df['target'] = pathway_relations_df['target'].apply(split_and_lower)

#Prepare a dictionary for fast lookup of gene expression values.
gene_expression_dict = gene_expression_df.to_dict(orient='index')


#Calculate activity and consistency of paths.
#Handle molecules (proteins, rnas and compouns).
#Complexes are proteins == basic complexes (built from links to probes) or group of proteins.
#Compounds are assumed to always be present (UDP == 1).
#If a molecule does not have a probability we need to remove the whole interaction from activity and consistency calculations.
#Function to calculate activity for a single interaction and a single sample.
def calculate_interaction_activity(row, sample):
    sources = row['source']
    interaction_type = row['interactiontype']
    
    # Calculate the product of gene values in the source column
    activity = 1
    for gene in sources:
        if gene.startswith('cpd:'):
            activity *= 1
        else:
            gene_value = gene_expression_dict.get(gene, {}).get(sample, 1)
            activity *= gene_value

    # Adjust activity for inhibition
    if 'inhibit' in interaction_type:
        activity = 1 - activity

    return activity

#Calculate activity for each sample and each pathway.
sample_names = gene_expression_df.columns
pathway_activities = []

for pathway in tqdm(pathway_relations_df['pathway'].unique()):
    pathway_data = pathway_relations_df[pathway_relations_df['pathway'] == pathway]
    pathway_activity = [pathway]
    
    for sample in sample_names:
        activities = pathway_data.apply(lambda row: calculate_interaction_activity(row, sample), axis=1)
        mean_activity = activities.mean()
        pathway_activity.append(mean_activity)
    
    pathway_activities.append(pathway_activity)

# Create a DataFrame for the results and save to a CSV file
columns = ['pathway'] + list(sample_names)
output_df = pd.DataFrame(pathway_activities, columns=columns)
output_df.to_csv(output_file, index=False)

print(f"Activity calculations saved to {output_file}")