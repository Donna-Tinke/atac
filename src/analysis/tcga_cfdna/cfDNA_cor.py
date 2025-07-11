# %%
import numpy as np
import pandas as pd
import glob
import os
from scipy.stats import pearsonr

# %%
peak_names = pd.read_csv('brca_correlaties/TCGA_data/filtered_peak_names.csv', header=None).values
ff_feature = 'cov_x'

# %%
ff_path = 'TSSClassification/results/04-15_10-29_all_coad'

# %%
def calc_correlations(ff_path, ff_feature):
        file_names = [f for f in os.listdir(ff_path) if f.endswith(".csv")]
        column_headers = [f.split('.')[0] for f in file_names]

        corr_mat = pd.DataFrame(
            np.nan,
            index= column_headers,
            columns=column_headers
        )

        for ff_file1 in file_names: ##Rows 
            cf_file_path = os.path.join(ff_path, ff_file1)
            cf_sample1 = pd.read_csv(cf_file_path, header=0, dtype='object', usecols=[ff_feature])
            cf_sample1 = cf_sample1.astype(float).values.flatten()

            for ff_file2 in file_names: ##Columns
                cf_file_path = os.path.join(ff_path, ff_file2)
                cf_sample2 = pd.read_csv(cf_file_path, header=0, dtype='object', usecols=[ff_feature])
                cf_sample2 = cf_sample2.astype(float).values.flatten()

                # Calculate the Pearson correlation
                corr, _ = pearsonr(cf_sample1, cf_sample2)
                
                # Store the correlation value in the matrix
                col_name1 = ff_file1.split('.')[0]
                col_name2 = ff_file2.split('.')[0]
                corr_mat.loc[col_name1, col_name2] = corr


        #Save the correlation matrix for a cancer type
        corr_mat.to_csv('correlation_results/cfDNA_only/coad_cfDNA_cor.csv', index= True)

# %%
cfDNA_cor = calc_correlations(ff_path, ff_feature)


