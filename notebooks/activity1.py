import numpy as np
import pandas as pd
import time
import multiprocessing as mp
import sys
import gc
from scipy.stats import mannwhitneyu as mann
from scipy.stats import f_oneway, levene
import os
from tqdm import tqdm


pathway_file = './data/pathay_elations.csv'

class path_activity:
    def __init__(self, udp):
        print(time.ctime(), 'Init activity object')
        self.chunk_num = 0
        self.orig_paths = pd.read_csv(pathway_file)
        self.udp = udp
        if len(udp) == 0:
            print('Empty UDP.')
            sys.exit()
    
    def calc_kegg_link_to_pr(self, sample):
        udp_sample = self.udp.loc[:, sample].copy().reset_index()
        udp_sample.columns = ['gene', 'pr']
        self.link_to_pr_dict[sample] = udp_sample.set_index('gene').pr.to_dict()

    #Calculate activity and consistency of paths.
    #Handle molecules (proteins, rnas and compouns).
    #Complexes are proteins == basic complexes (built from links to probes) or group of proteins.
    #Compounds are assumed to always be present (UDP == 1).
    def process_samples(self, udp_chunk):
        gc.collect()
        print(self.chunk_num, end="")
        self.chunk_num = self.chunk_num + 1000
        #Table columns are (sampleID, path_id, activity, consistency).
        sample_results = np.empty((0, 6))
        for sample in tqdm(udp_chunk):
            #Calculate UDP of the molecules in all the paths.
            if self.kegg: #New pathway that was added by the new procedure.
               self.calc_kegg_link_to_pr(sample)
            else:
               self.calc_link_to_gene_to_pr(sample) #Single cell is always rnaseq.
            cmplx_to_pr_dict = self.calc_cmplx_to_pr(sample)
            paths = self.orig_paths.copy()

            #If a molecule does not have a probability we need to remove the whole interaction from activity and consistency calculations.
            paths.loc[paths.molType == 'protein', 'pr'] = paths.molLink.apply(lambda x: max([self.link_to_udp(i, sample) for i in str(x).split(',')]))  # Max returns nan if there is at least one nan in the list.
            #Compounds are assumed to always be present.
            paths.loc[paths.molType == 'compound', 'pr'] = 1
            paths.loc[paths.molType == 'rna', 'pr'] = 1
            paths.loc[paths.molType == 'complex', 'pr'] = paths.molNum.apply(lambda x: cmplx_to_pr_dict.get(x, np.nan))

            #Calculate activity and consistency of each interaction. Activity of an interaction is a multiplication of the inputs. Note that some interactions participate in few paths.
            paths.loc[paths.molRole == 'inhibitor', 'pr'] = 1 - paths.pr

            #Save the pr for output molecules, since it also has a UDP. First allocate a new column.
            paths['pr_output'] = 1
            paths['pr_output'] = paths['pr_output'].astype(float)
            paths.loc[paths.molRole == 'output', 'pr_output'] = paths.pr
            
            #Now assign 1 to the output molecules pr. By the next block, interactions with NA, will be removed after multiplication.
            paths.loc[paths.molRole == 'output', 'pr'] = 1
            
            # Handle molecules that we could not find UDP for, we need to remove the whole interaction. Pandas like R is ignoring NA values in gropuby.transform, so instead we zero it.
            # paths = paths.dropna(subset=['pr'])
            paths = paths.fillna({'pr': 0})
            paths['interaction_activity'] = paths.groupby(['path_id', 'intID'])['pr'].transform('prod')
            paths['interaction_consistency'] = paths.interaction_activity * paths.pr_output + (1 - paths.interaction_activity) * (1 - paths.pr_output)

            # Calculate activity and consistency of each path.
            # Both activity and consistency are averages of all interactions activities and consistencies.
            paths['activity'] = paths.groupby('path_id')['interaction_activity'].transform('mean')
            paths['consistency'] = paths.groupby('path_id')['interaction_consistency'].transform('mean')
            paths['sampleID'] = sample
            paths_result = paths[['sampleID', 'path_name', 'path_id', 'activity', 'consistency', 'molRole']].drop_duplicates()
            sample_results = np.concatenate((sample_results, paths_result.values))  # axis=0
        
        results_df = pd.DataFrame(data=sample_results, columns=['sampleID', 'path_name', 'pathID', 'Activity', 'Consistency', 'molRole'])
        results_df.drop(['molRole'], inplace=True, axis=1)
        results_df.drop_duplicates(inplace=True)
        results_df = pd.pivot(results_df, columns='sampleID', index='path_name', values='Activity')
        results_df.to_csv(f'./data/activity/{self.chunk_num}.csv')
        sys.stdout.flush()
        return results_df, paths

    # Read chunks of dataframe columns.
    def chunker_columns(self, size):
        len_df = len(self.udp.columns)
        cols = self.udp.columns
        return [cols[pos:min(pos + size, len_df)] for pos in range(0, len_df, size)]

    # Run activity on parallel.
    def calc_activity_consistency(self):
        gc.collect()
        print(time.ctime(), 'Calculate activity and consistency...')
        #df = pd.DataFrame()
        '''
        pool = mp.Pool()  # Configure number of CPUs processes.
        results = [pool.apply_async(self.process_samples, args=(x,)) for x in self.chunker_columns(1000)]
        for p in results:
            p.get(timeout=100)
            #df = pd.concat([df, p.get()[0]]) # f.get(timeout=100)
            print('.', end="")
            sys.stdout.flush()
        '''
        #df = self.process_samples(self.udp)[0]
        #df['Activity'] = df.Activity #Consistency
        for x in tqdm(self.chunker_columns(1000)):
            self.process_samples(x)
        
        #List all CSV files in the directory.
        csv_files = [file for file in os.listdir('./data/activity/') if file.endswith('.csv')]

        #Check if there are CSV files to process.
        if len(csv_files) == 0:
            print("No CSV files found in the directory.")
            return pd.DataFrame()
        
        #Initialize a flag to keep track of whether we are processing the first file.
        first_file = True
        merged_df = None

        #Iterate over each CSV file.
        for file in csv_files:
            file_path = os.path.join('./data/activity/', file)
    
            #Read the CSV file into a pandas DataFrame.
            df = pd.read_csv(file_path, index_col=0)
            dfT = df.T
            #If it's the first file, initialize the merged_df.
            if first_file:
                merged_df = dfT
                first_file = False
            else:
                #Concatenate subsequent files to the merged_df along the row axis (axis=0).
                merged_df = pd.concat([merged_df, dfT], axis=0)

        #Export the merged DataFrame to a new CSV file.
        merged_df.to_csv('output_activity.csv')
        #Delete all the CSV files from the directory.
        #for file in csv_files:
        #    os.remove(file_path)

        #pool.close()
        self.activity = merged_df
        print(time.ctime(), "Done.")
        return merged_df


    def calc_mann_whitney(self, pathname, file1, file2):
        data1 = pd.read_csv(file1)
        data1 = data1.loc[data1.path_name == pathname]
        data2 = pd.read_csv(file2)
        data2 = data2.loc[data2.path_name == pathname]
        n1 = len(data1)
        n2 = len(data2)
        m_u = n1 * n2 * 0.5
        stat, p = mann(data1.Activity, data2.Activity)
        # sigma_u = np.sqrt(n1*n2*(n1+n2+1)*(10/12))
        print(
            f'Stat: {stat}, P-value: {p}, n1: {n1}, n2: {n2}, m_u(expected under H0): {m_u}')


def calc_activity_from_adata(adata):
    #On the index we have cell names (adata.obs_names) and the columns are gene names (adata.var_names). We transpose this dataframe before calculating activity.
    df = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names).T
    df.index = df.index.map(str.lower)
    activity_obj = path_activity(df)
    return activity_obj.calc_activity_consistency()


if __name__ == '__main__':
    udp = pd.read_csv('./data/output_udp.csv', index_col=0)
    udp.index = udp.index.map(str.lower)
    activity_obj = path_activity(udp)
    activity_obj.calc_activity_consistency()