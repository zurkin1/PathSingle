# %%
import pandas as pd
from tqdm import tqdm
import gc
import numpy as np
#import cProfile
#import pstats
from concurrent.futures import ProcessPoolExecutor, as_completed
import scanpy as sc
from scvi.data import AnnDataManager
from scvi.dataloaders import AnnDataLoader


class path_activity:
    def __init__(self, adata):
        #Load the necessary data.
        self.pathway_relations_file = './data/pathway_relations.csv'
        self.output_file = './data/output_activity.csv'

        #Split the 'source' and 'target' columns by '*' and lowercase the gene names.
        def split_and_lower(x):
            if isinstance(x, str):
                return x.lower().split('*')
            return []

        #Load pathway relations data.
        self.pathway_relations_df = pd.read_csv(self.pathway_relations_file)
        self.pathway_relations_df['source'] = self.pathway_relations_df['source'].apply(split_and_lower)
        #self.pathway_relations_df['target'] = self.pathway_relations_df['target'].apply(split_and_lower)

        #Filter gene_expression_df to include only genes from pathway_relations data.
        #genes_in_source = set(gene for sublist in self.pathway_relations_df['source'] for gene in sublist)

        #Load gene expression data.
        #self.gene_expression_df = udp
        #self.gene_expression_df.index = self.gene_expression_df.index.map(str.lower)
        self.adata = adata
        self.adata.var_names = self.adata.var_names.str.lower()

        #Round the values to 5 decimal places.
        #self.gene_expression_df = self.gene_expression_df.round(5)
        #self.adata.X = np.round(self.adata.X, 5)
        sc.pp.scale(self.adata, max_value=1e-5, copy=True)

        #self.gene_expression_df = self.gene_expression_df[self.gene_expression_df.index.isin(genes_in_source)]
        #self.adata = self.adata[:, self.adata.var_names.isin(genes_in_source)]

        #Initialize AnnDataManager.
        self.adata_manager = AnnDataManager(adata=self.adata)
        self.adata_manager.register_fields()

        
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
    #https://kernprof.readthedocs.io/en/latest/
    def calculate_interaction_activity(self, sources, interaction_type, sample_name):
        activity = 1.0
        for gene in sources:
            #gene_value = self.gene_expression_df.at[gene, sample] if gene in self.gene_expression_df.index else 0.0
            if gene in self.adata.var_names:
                gene_idx = self.adata.var_names.get_loc(gene)
                sample_idx = self.adata.obs_names.get_loc(sample_name)
                gene_value = self.adata.X[sample_idx, gene_idx]
            else:
                gene_value = 0.0
            activity *= gene_value

        #Adjust activity for different relation subtypes.
        if self.is_inhibitory(interaction_type):
           activity = 1 - activity

        return activity

    def process_sample_batch(self, batch):
        try:
            for sample_idx in batch:
                sample_name = self.adata.obs_names[sample_idx]
                sample_activities = []
                for pathway in self.pathway_relations_df['pathway'].unique():
                    pathway_data = self.pathway_relations_df[self.pathway_relations_df['pathway'] == pathway]
                    activities = [
                        self.calculate_interaction_activity(row['source'], row['interactiontype'], sample_idx)
                        for _, row in pathway_data.iterrows()
                    ]
                    pathway_activity = np.mean(activities, axis=0)
                    sample_activities.append([pathway, sample_name, pathway_activity])

                df = pd.DataFrame(sample_activities, columns=['pathway', 'sample', 'activity'])
                df.to_csv(f'./data/samples/{sample_name}_activity.csv', index=False)
            return True
        except Exception as e:
            print(f"Error processing batch: {e}")
            return False

    def calculate_activity(self):
        #Calculate activity for each sample and each pathway.
        #sample_names = self.gene_expression_df.columns
        # Initialize the AnnDataLoader
        data_loader = AnnDataLoader(
            adata_manager=self.adata_manager,
            batch_size=128,  # Adjust based on your memory constraints
            shuffle=False
        )

        try:

            success_samples = []
            failures = 0
            successes = 0
        
            #pool = mp.Pool(processes=10)
            #results = [pool.apply_async(self.process_sample, args=(i, sample,)) for i, sample in enumerate(sample_names)]
            
            #Use multiprocessing to process samples.
            with ProcessPoolExecutor(max_workers=20) as executor:
                futures = {executor.submit(self.process_sample_batch, batch): batch for batch in data_loader}

                for future in tqdm(as_completed(futures), total=len(futures)):
                    batch = futures[future]
                    try:
                        success = future.result()
                        if success:
                            success_samples.append(batch)
                            successes += 1
                        else:
                            failures += 1
                    except Exception as e:
                        print(f"Error processing batch: {e}.")
                        failures += 1

            print(f"Successfully processed {successes} batches, {failures} batches failed.")

            #Merge all the individual CSV files into one DataFrame.
            all_dfs = []
            for batch in success_samples:
                for sample in batch:
                    df = pd.read_csv(f'./data/samples/{sample}_activity.csv')
                    all_dfs.append(df)

            #Create a DataFrame for the results and save to a CSV file.        
            pathway_activities_df = pd.concat(all_dfs, ignore_index=True)
            self.output_df = pathway_activities_df.pivot(index='pathway', columns='sample', values='activity')
            self.output_df.to_csv(self.output_file)
            print(f"Activity calculations saved to {self.output_file}")
        
        except Exception as e:
            print(f"Error in calculate_activity: {e}.")

       
def calc_activity_from_adata(adata):
    #On the index we have cell names (adata.obs_names) and the columns are gene names (adata.var_names). We transpose this dataframe before calculating activity.
    #df = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names).T
    activity_obj = path_activity(adata)
    activity_obj.calculate_activity()

# %%
if __name__ == '__main__':
    #udp = pd.read_csv('./data/sample_file.csv', index_col=0)
    #udp.index = udp.index.map(str.lower)
    #activity_obj = path_activity(udp)
    #activity_obj.calculate_activity()
    
    '''
    profile_filename = './data/profile_output.prof'
    cProfile.run('activity_obj.calculate_activity()', profile_filename)

    with open('./data/profile_stats.txt', 'w') as f:
        stats = pstats.Stats(profile_filename, stream=f)
        stats.strip_dirs().sort_stats('cumulative').print_stats(10)
    '''
    adata = sc.read_h5ad('./data/kang_counts_25kmagic.h5ad')
    adata = sc.pp.subsample(adata, fraction=0.4, copy=True)
    print('Finished reading data.')
    calc_activity_from_adata(adata)
# %%
'''
import pandas as pd


udp = pd.read_csv('./data/sample_file.csv', index_col=0)
udp_large = pd.concat([udp] * 50, axis=1)
udp_large.columns = [f"{col}_{i}" for i in range(50) for col in udp.columns]
udp_large.to_csv('./data/sample_file_large.csv')
'''
# %%