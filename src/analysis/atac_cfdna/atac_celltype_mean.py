# %%
import numpy as np
import pandas as pd
import glob
import os
import re


# %%
atac_path = 'data/atac_files/small'

bed_df = pd.read_csv("TSSClassification/beds/dts_beds/atac_match.bed", sep='\t', header=None)
gene_names_in_bed = set(bed_df[3].unique())

# %%
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

    tissue_data['cell_type'] = cell_types

    # Keep only numeric columns + cell_type
    numeric_cols = tissue_data.select_dtypes(include=[np.number]).columns.tolist()
    group_input = tissue_data[numeric_cols + ['cell_type']]

    # Group and compute mean
    mean_matrix = group_input.groupby('cell_type').mean()

    #Filter for genes that match the bed
    filtered = mean_matrix.reset_index()
    matching_cols = [col for col in filtered.columns if col in gene_names_in_bed or col == 'cell_type']
    filtered = filtered[matching_cols]

    filtered.to_csv(f'data/atac_files/cell_type_mean/{tissue}.csv', index=False)


# %%



