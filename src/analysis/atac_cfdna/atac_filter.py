import numpy as np
import pandas as pd
import glob
import os


atac_path = 'data/atac_files/cell_type_mean'
bed_df = pd.read_csv("TSSClassification/beds/dts_beds/atac_match.bed", sep='\t', header=None)

gene_names_in_bed = set(bed_df[3].unique())

for my_file in os.listdir(atac_path):
    if not my_file.endswith('.csv'):
        continue
    tissue_data = pd.read_csv(os.path.join(atac_path, my_file))
    matching_cols = [col for col in tissue_data.columns if col in gene_names_in_bed]
    tissue_data = tissue_data[matching_cols]

    tissue_data.to_csv(f'data/atac_files/cell_type_mean/{my_file}')