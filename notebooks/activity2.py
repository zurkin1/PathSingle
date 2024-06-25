import pandas as pd
from tqdm import tqdm
import numpy as np


tqdm.pandas()

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


    #Function to get the expression value of a gene for a specific sample.
    def get_expression_value(self, gene, sample):
        if gene.startswith('cpd:'):
            return 1.0  #Compounds have a value of 1.
        try:
            return self.gene_expression_df.at[gene, sample]
        except KeyError:
            return 0.0  #Assume missing genes have a value of 1.


    #Function to determine if an interaction type is inhibitory.
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
    def calculate_interaction_activity(self, row, sample):
        sources = row['source']
        interaction_type = row['interactiontype']
        
        #Calculate the product of gene values in the source column.
        activity = 1.0
        for gene in sources:
            if gene.startswith('cpd:'):
                gene_value = 1.0
            else:
                gene_value = self.get_expression_value(gene, sample)
            
            activity *= gene_value

        #Adjust activity for different relation subtypes.
        if self.is_inhibitory(interaction_type):
           activity = 1 - activity

        return activity


    def calculate_activity(self):
        #Calculate activity for each sample and each pathway.
        sample_names = self.gene_expression_df.columns
          

        def calculate_for_pathway_and_sample(pathway_data, sample):
            return pathway_data.apply(lambda row: self.calculate_interaction_activity(row, sample), axis=1)


        def calculate_for_pathway(pathway):
            pathway_data = self.pathway_relations_df[self.pathway_relations_df['pathway'] == pathway]
            activities = sample_names.to_series().apply(lambda sample: calculate_for_pathway_and_sample(pathway_data, sample)).values.tolist()
            return [pathway] + activities
            

        pathway_activities = self.pathway_relations_df['pathway'].unique().tolist()
        results = pd.Series(pathway_activities).progress_apply(calculate_for_pathway)

        #Create a DataFrame for the results and save to a CSV file.
        columns = ['pathway'] + list(sample_names)
        output_df = pd.DataFrame(results.tolist(), columns=columns).T
        output_df.to_csv(self.output_file, index=False)

        print(f"Activity calculations saved to {self.output_file}")


def calc_activity_from_adata(adata):
    #On the index we have cell names (adata.obs_names) and the columns are gene names (adata.var_names). We transpose this dataframe before calculating activity.
    df = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names).T
    activity_obj = path_activity(df)
    activity_obj.calculate_activity()