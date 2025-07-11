# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

# %%
coad_ff_path = 'TSSClassification/results/04-15_10-29_all_coad'
brca_ff_path = 'TSSClassification/results/04-15_10-32_all_brca'
peak_names = pd.read_csv('brca_correlaties/TCGA_data/filtered_peak_names.csv', header=None).values

# %%
ff_feature = 'cov_x'

# %%
coad_file_names = [f for f in os.listdir(coad_ff_path) if f.endswith(".csv")]
brca_file_names = [f for f in os.listdir(brca_ff_path) if f.endswith(".csv")]

ff_list = []  # ← moved outside the loop

for ff_file in coad_file_names:
    cf_file_path = os.path.join(coad_ff_path, ff_file)
    new_column = pd.read_csv(cf_file_path, header=0, dtype='object', usecols=[ff_feature])
    new_column = new_column.astype(float).values.flatten()
    ff_list.append(new_column)

for ff_file in brca_file_names:
    cf_file_path = os.path.join(brca_ff_path, ff_file)
    new_column = pd.read_csv(cf_file_path, header=0, dtype='object', usecols=[ff_feature])
    new_column = new_column.astype(float).values.flatten()
    ff_list.append(new_column)    

# Now stack all at once
ff_matrix = None
ff_matrix = np.stack(ff_list)

# %%
peak_var = np.empty((len(peak_names), 2), dtype=object)
peak_var[:, 0] = peak_names.flatten()

# %%
for i in range(len(peak_var)):
    peak_var[i, 1] = np.var(ff_matrix[:, i])

# %%
peak_var_df = pd.DataFrame(peak_var)
top_peaks = peak_var_df.sort_values(by=peak_var_df.columns[1], ascending=False).head(10000)
high_var_peaks = top_peaks.values[:,0]

# %%
np.savetxt('brca_correlaties/TCGA_data/high_var_covx_peaks.csv', high_var_peaks, delimiter=",", fmt="%s")


