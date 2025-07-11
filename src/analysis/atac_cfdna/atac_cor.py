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
# Bed 
bed = pd.read_csv("TSSClassification/beds/dts_beds/atac_match.bed", sep='\t', header=None)

# %%
# cfDNA
ff_path = 'TSSClassification/results/06-05_13-04_COAD_atac'
feature = 'cov_spread'
file_names = [f for f in os.listdir(ff_path) if f.endswith(".csv")]

# %%
# atac
atac_df = pd.read_csv("tacje/tissue_mean_atac_signals.csv", index_col=0, header=0)

# %%
bed_genes = bed.iloc[:, 3].tolist()
filtered_atac = atac_df.loc[:, atac_df.columns.isin(bed_genes)]

# %%
cf_samples = []

# Load each cfDNA file and append to list
for ff_file in file_names:
    cf_file_path = os.path.join(ff_path, ff_file)
    cf_sample = pd.read_csv(cf_file_path, header=0, dtype='object', usecols=[feature])
    cf_sample = cf_sample.astype(float).values.flatten()
    cf_samples.append(cf_sample)
    ## Elke row is een sample

sample_names = [re.match(r"^[^\.]+", f).group(0) for f in file_names]
cf_df = pd.DataFrame(cf_samples, index=sample_names)

# %%
cf_z = (cf_df - cf_df.mean(axis=1).values[:, None]) / cf_df.std(axis=1).values[:, None]
atac_z = (filtered_atac - filtered_atac.mean(axis=1).values[:, None]) / filtered_atac.std(axis=1).values[:, None]

# Compute pairwise Pearson correlation (dot product of z-scores)
corr_matrix = np.dot(cf_z.values, atac_z.values.T) / cf_z.shape[1]
corr_matrix = np.abs(corr_matrix)

# Convert to DataFrame with proper labels
corr_df = pd.DataFrame(corr_matrix, index=cf_df.index, columns=atac_df.index)

# %%
print(corr_df)

# %%
ranked_df = corr_df.rank(axis=1, ascending=False, method='min')


# %%
corr_df.to_csv("tacje/coad_abs_cor_mat.csv")
ranked_df.to_csv("tacje/coad_ranked_cor_mat.csv")


