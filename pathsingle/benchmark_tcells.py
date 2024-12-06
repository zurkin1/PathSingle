# -*- coding: utf-8 -*-
import scanpy as sc
import pandas as pd
from metrics import *
from activity import *
import scprep
import magic
from warnings import filterwarnings
import numpy as np


adata = sc.read_h5ad('./data/sc_training.h5ad')
# Load the data.
num_splits = 5
split_files = [f'./data/sc_training_split_{i+1}.h5ad' for i in range(num_splits)]
splits = [sc.read_h5ad(file) for file in split_files]

# Concatenate the splits back into a single AnnData object.
adata = splits[0].concatenate(*splits[1:], batch_key='batch', batch_categories=[f'batch_{i+1}' for i in range(num_splits)])

#adata = sc.pp.subsample(adata, fraction=0.3, copy=True)
#28697 cells × 15077 genes.
print(adata)
adata.X = adata.layers['rawcounts'].copy()

activity_df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
activity_df = scprep.normalize.library_size_normalize(activity_df)
activity_df = scprep.transform.sqrt(activity_df)

magic_op = magic.MAGIC()
activity_df = magic_op.fit_transform(activity_df)
activity_df = activity_df.astype(np.float16)
activity_df.to_csv('./data/activity_df.csv', index=True)

# Run PathSingle.
activity = sc.read('./data/activity_df.csv', delimiter=',', cache=False)
calc_activity(activity)