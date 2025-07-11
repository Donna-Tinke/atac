# %%
import numpy as np
import pandas as pd
import glob
import os
from scipy.sparse import csr_matrix
import re

# %%
ff_path = 'TSSClassification/results/06-04_12-12_healthy_atac'
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

# %%
sample_ids = [os.path.splitext(f)[0] for f in file_names]                # strip “.csv”
df = pd.DataFrame(cf_samples, index=sample_ids)                           # rows = samples     # name columns cov_spread_1,…  
df.to_csv("data/cfdna/ataccfmatje.csv", index_label='sample_id')

# %%
print(df)

# %%



