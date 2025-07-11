# %%
import numpy as np
import pandas as pd
import glob
import os
from scipy.stats import pearsonr

# %%
peak_names = pd.read_csv('brca_correlaties/TCGA_data/high_var_peaks.csv', header=None).values
ff_feature = 'cov_x'

# %%
# ATAC
sample_ids = pd.read_csv('brca_correlaties/TCGA_data/coad_sample_ids.csv', header=None).values
atac_data =  pd.read_csv('brca_correlaties/TCGA_data/high_var_coad.csv', header=None).values

# %%
# cfDNA
ff_path = 'TSSClassification/results/04-16_10-42_high_peakvar_brca'

# %%
## Funtcion to check for constan values
def is_constant(array):
    """Check if a NumPy array contains constant values (all elements are equal)."""
    return np.all(array == array[0])

# %%
def calc_correlations(atac_data, atac_samples, ff_path, ff_feature):
        file_names = [f for f in os.listdir(ff_path) if f.endswith(".csv")]
        column_headers = [f.split('.')[0] for f in file_names]

        flat_sample_ids = [x[0] if isinstance(x, (list, tuple, np.ndarray)) else x for x in atac_samples]

        corr_mat = pd.DataFrame(
            np.nan,
            index= flat_sample_ids,
            columns=column_headers
        )
        
        
        for ff_file in file_names:
            cf_file_path = os.path.join(ff_path, ff_file)
            cf_sample = pd.read_csv(cf_file_path, header=0, dtype='object', usecols=[ff_feature])
            # cf_sample = cf_sample.drop(index=0).reset_index(drop=True)
            cf_sample = cf_sample.astype(float).values.flatten()
            
            for i in range(atac_data.shape[1]):
                atac_sample = atac_data[:, i]
                # atac_sample, cf_sample = atac_sample.align(ff_sample, join='inner') 
            
                # Skip correlation calculation if one of the columns has NaN, or is constant
                if np.isnan(atac_sample).all():
                    print(f"Skipping an atac sample due to all NaN values.")
                    continue

                if np.isnan(cf_sample).all():
                    print(f"Skipping a fragment feature file due to all NaN values.")
                    continue

                if is_constant(atac_sample):
                    print(f"Skipping an atac sample due to constant values.")
                    continue

                if is_constant(cf_sample):
                    print(f"Skipping a fragment feature file due to constant values.")
                    continue
                
                # Calculate the Pearson correlation
                corr, _ = pearsonr(atac_sample, cf_sample)
                
                # Store the correlation value in the matrix
                col_name = ff_file.split('.')[0]
                corr_mat.loc[flat_sample_ids[i], col_name] = corr
        
        #Save the correlation matrix for a cancer type
        corr_mat.to_csv('correlation_results/all_samples/high_peak_var/peakvar_coad_atac_brca_cfDNA.csv', index= True) 

# %%
atac_cfDNA = calc_correlations(atac_data, sample_ids, ff_path, ff_feature)


