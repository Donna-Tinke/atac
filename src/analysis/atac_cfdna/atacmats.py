# %%
import numpy as np
import pandas as pd
import glob
import os
from scipy.stats import pearsonr
import anndata as ad
from scipy.sparse import csr_matrix
from scipy.stats import zscore
import time
from scipy.sparse import issparse
import re

# %%
atac_path = 'data/atac_files/small'

# %%
ann_genes = pd.read_csv('atac_expres/ann_genes.csv')

# %%
ann_genes = ann_genes['feature_name'].astype(str)

# %%
##Eerst loopen over npy files en dictionary maken
npy_file_map = {}
for fname in os.listdir(atac_path):
    if not fname.endswith('.npy'):
        continue
    parts = re.split(r'atac_nearest_cell_types_|_harmony.npy', fname)
    if len(parts) > 3:
        tissue = parts[1] + '_' + parts[2]
    else:
        tissue = parts[1]
    npy_file_map[tissue] = fname


#Dan loopen over csvs en csv waar geen npy voor is skippen
for csv_file in os.listdir(atac_path):
    if not csv_file.endswith('.csv'):
        continue

    tissue_data = pd.read_csv(os.path.join(atac_path, csv_file))

    #Tssue name fetchen
    tissue = csv_file.replace("tissue_subset_", "").replace(".csv", "")

    print(f"\n▶ Processing tissue: {tissue}")

    #npy file matchen
    if tissue not in npy_file_map:
        print(f" No matching .npy file for {tissue}")
        continue

    types_file = npy_file_map[tissue]
    cell_types = np.load(os.path.join(atac_path, types_file))
    cell_types = cell_types.astype(str) 

    #Matching genes filteren
    tissue_data = tissue_data.loc[:, tissue_data.columns.isin(ann_genes)]
    print(f"Filtered matrix shape: {tissue_data.shape}")

    #Group by cell type
    tissue_data['cell_type'] = cell_types
    mean_matrix = tissue_data.groupby('cell_type').mean()

    #Cell type labels saven
    cell_type_df = pd.DataFrame({'cell_type': cell_types})
    cell_type_df.to_csv(f'data/atac_files/{tissue}_cell_types_raw.csv', index=False)

   
    out_path = f'data/atac_files/celltype_mats/{tissue}.csv'
    mean_matrix.to_csv(out_path)
    print(f"Saved mean matrix of {tissue}")