# %%
import pandas as pd
from tqdm import tqdm
import gc
import numpy as np
#import cProfile
#import pstats
import multiprocessing as mp
import scanpy as sc
import scanpy.external as sce


class path_activity:
    def __init__(self, udp):
        #Load the necessary data.
        self.pathway_relations_file = './data/pathway_relations.csv'
        self.output_file = './data/output_activity.csv'

        #Load gene expression data.
        self.gene_expression_df = udp
        self.gene_expression_df.index = self.gene_expression_df.index.map(str.lower)


        #Split the 'source' and 'target' columns by '*' and lowercase the gene names.
        def split_and_lower(x):
            if isinstance(x, str):
                return x.lower().split('*')
            return []


        #Load pathway relations data.
        self.pathway_relations_df = pd.read_csv(self.pathway_relations_file)
        self.pathway_relations_df['source'] = self.pathway_relations_df['source'].apply(split_and_lower)
        self.pathway_relations_df['target'] = self.pathway_relations_df['target'].apply(split_and_lower)


    # Function to determine if an interaction type is inhibitory
    def is_inhibitory(self, interaction_type):
        inhibitory_keywords = ['inhibition', 'repression', 'dissociation', 'missing interaction', 'dephosphorylation', 'ubiquitination']
        return any(keyword in interaction_type for keyword in inhibitory_keywords)

    #Calculate activity and consistency of paths.
    #Handle molecules (proteins, rnas and compouns).
    #Complexes are proteins == basic complexes (built from links to probes) or group of proteins.
    #Compounds are assumed to always be present (UDP == 1).
    #If a molecule does not have a probability we need to remove the whole interaction from activity and consistency calculations.
    #Reference https://www.kegg.jp/kegg/xml/docs/.
    #Function to calculate activity for a single interaction and a single sample.
    #Calculate the product of gene values in the source column.
    def calculate_interaction_activity(self, sources, interaction_type, sample):
        activity = 1.0
        for gene in sources:
            gene_value = self.gene_expression_df.at[gene, sample] if gene in self.gene_expression_df.index else 0.0
            activity *= gene_value

        #Adjust activity for different relation subtypes.
        if self.is_inhibitory(interaction_type):
           activity = 1 - activity

        return activity


    def process_sample(self, sample):
        gc.collect()
        sample_activities = []
        for pathway in self.pathway_relations_df['pathway'].unique():
            #Define dataframe for pathway.
            pathway_data = self.pathway_relations_df[self.pathway_relations_df['pathway'] == pathway]
            #activities = pathway_data.apply(lambda row: self.calculate_interaction_activity(row['source'], row['interactiontype'], sample), axis=1)
            activities = [
                            self.calculate_interaction_activity(row['source'], row['interactiontype'], sample)
                            for _, row in pathway_data.iterrows()
                        ]
            pathway_activity = np.mean(activities, axis=0)
            sample_activities.append([pathway, sample, pathway_activity])
        return sample_activities


    def calculate_activity(self):
        #Calculate activity for each sample and each pathway.
        sample_names = self.gene_expression_df.columns
        pathway_activities = []
    
        pool = mp.Pool()
        results = [pool.apply_async(self.process_sample, args=(sample,)) for sample in sample_names]
        
        for result in tqdm(results, total=len(sample_names)):
            pathway_activities.extend(result.get(timeout=100))

        #Create a DataFrame for the results and save to a CSV file.        
        pathway_activities_df = pd.DataFrame(pathway_activities, columns=['pathway', 'sample', 'activity'])
        self.output_df = pathway_activities_df.pivot(index='pathway', columns='sample', values='activity')
        self.output_df.to_csv(self.output_file)
        print(f"Activity calculations saved to {self.output_file}")

        
def calc_activity_from_adata(adata):
    #On the index we have cell names (adata.obs_names) and the columns are gene names (adata.var_names). We transpose this dataframe before calculating activity.
    df = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names).T
    activity_obj = path_activity(df)
    activity_obj.calculate_activity()

# %%
if __name__ == '__main__':
    udp = pd.read_csv('./data/sample_file.csv', index_col=0)
    udp.index = udp.index.map(str.lower)
    activity_obj = path_activity(udp)
    activity_obj.calculate_activity()
    '''
    profile_filename = './data/profile_output.prof'
    cProfile.run('activity_obj.calculate_activity()', profile_filename)

    with open('./data/profile_stats.txt', 'w') as f:
        stats = pstats.Stats(profile_filename, stream=f)
        stats.strip_dirs().sort_stats('cumulative').print_stats(10)
    '''
# %%
adata = sc.read_h5ad('./data/sc_trainingmagic.h5ad')
calc_activity_from_adata(adata)