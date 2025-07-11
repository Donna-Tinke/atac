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
##61759 genes
##19376 bed length

# %%
expr_mats_dir = "data/scrna/expr_mats/stom"
bed_df = pd.read_csv("TSSClassification/beds/dts_beds/tabula_match.bed", sep='\t', header=None)

# %%
# cfDNA
ff_path = 'TSSClassification/results/04-15_10-25_cf_tabula'
ff_feature = 'cov_spread'

# %%
file_names = [f for f in os.listdir(ff_path) if f.endswith(".csv")]

# %%
cf_samples = []

# Load each cfDNA file and append to list
for ff_file in file_names:
    cf_file_path = os.path.join(ff_path, ff_file)
    cf_sample = pd.read_csv(cf_file_path, header=0, dtype='object', usecols=[ff_feature])
    cf_sample = cf_sample.astype(float).values.flatten()
    cf_samples.append(cf_sample)
    ## Elke row is een sample

sample_names = [re.match(r"^[^\.]+", f).group(0) for f in file_names]
cf_df = pd.DataFrame(cf_samples, index=sample_names)
cf_z = (cf_df - cf_df.mean(axis=1).values[:, None]) / cf_df.std(axis=1).values[:, None]

# %%
print(cf_z.head())

# %%
for file in os.listdir(expr_mats_dir):
    if file.endswith(".csv"):
        expr_path = os.path.join(expr_mats_dir, file)
        expr_df = pd.read_csv(expr_path, index_col=0)

        # Z-score normalization of expr_df by row
        expr_z = (expr_df - expr_df.mean(axis=1).values[:, None]) / expr_df.std(axis=1).values[:, None]

        # Ensure same number of columns
        if cf_z.shape[1] != expr_z.shape[1]:
            print(f"Skipping {file}: mismatched number of features (columns).")
            continue

        # Compute correlation matrix: dot product of z-scored rows
        correlation_matrix = np.abs(np.dot(cf_z.values, expr_z.values.T) / cf_z.shape[1])

        # Build DataFrame
        corr_df = pd.DataFrame(
            correlation_matrix,
            index=cf_df.index,
            columns=expr_df.index
        )

        # Save output
        out_path = os.path.join("correlation_results/scrna/control/filtered_cells/cov_spread", file)
        corr_df.to_csv(out_path)


