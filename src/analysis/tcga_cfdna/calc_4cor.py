# %%
import numpy as np
import pandas as pd
import glob
import os
from scipy.stats import pearsonr

# %%
peak_names = pd.read_csv('brca_correlaties/TCGA_data/high_var_covx_peaks.csv', header=None).values
ff_feature = 'cov_x'

# %%
# ATAC
coad_sample_ids = pd.read_csv('brca_correlaties/TCGA_data/coad_sample_ids.csv', header=None).values
coad_atac_data =  pd.read_csv('brca_correlaties/TCGA_data/high_covxvar_coad.csv', header=None).values
brca_sample_ids = pd.read_csv('brca_correlaties/TCGA_data/brca_sample_ids.csv', header=None).values
brca_atac_data =  pd.read_csv('brca_correlaties/TCGA_data/high_covxvar_brca.csv', header=None).values

# %%
# cfDNA
brca_cf_data = 'TSSClassification/results/04-21_12-08_high_covxvar_brca'
coad_cf_data = 'TSSClassification/results/04-21_12-09_high_covxvar_coad'

# %%
bed_df = pd.read_csv("TSSClassification/beds/dts_beds/covx_high_var_peaks.bed", sep='\t', header=None)

# %%
## Funtcion to check for constan values
def is_constant(array):
    """Check if a NumPy array contains constant values (all elements are equal)."""
    return np.all(array == array[0])

# %%

#Flattening the atac sample ids to function as indexes of the corr matrix
coad_flat_sample_ids = [x[0] if isinstance(x, (list, tuple, np.ndarray)) else x for x in coad_sample_ids]
brca_flat_sample_ids = [x[0] if isinstance(x, (list, tuple, np.ndarray)) else x for x in brca_sample_ids]

#Fethcing the file name sto function as column headers for the corr matrix
brca_file_names = [f for f in os.listdir(brca_cf_data) if f.endswith(".csv")]
brca_column_headers = [f.split('.')[0] for f in brca_file_names]

coad_file_names = [f for f in os.listdir(coad_cf_data) if f.endswith(".csv")]
coad_column_headers = [f.split('.')[0] for f in coad_file_names]

#Loading all the cfDNA data
coad_cf_array = np.zeros((len(bed_df), len(coad_file_names)))
for j, ff_file in enumerate(coad_file_names):
    cf_file_path = os.path.join(coad_cf_data, ff_file)
    cf_sample = pd.read_csv(cf_file_path, header=0, dtype='object', usecols=[ff_feature])
    # cf_sample = cf_sample.drop(index=0).reset_index(drop=True)
    cf_sample = cf_sample.astype(float).values.flatten()
    coad_cf_array[:, j] = cf_sample

brca_cf_array = np.zeros((len(bed_df), len(brca_file_names)))
for j, ff_file in enumerate(brca_file_names):
    cf_file_path = os.path.join(brca_cf_data, ff_file)
    cf_sample = pd.read_csv(cf_file_path, header=0, dtype='object', usecols=[ff_feature])
    # cf_sample = cf_sample.drop(index=0).reset_index(drop=True)
    cf_sample = cf_sample.astype(float).values.flatten()
    brca_cf_array[:, j] = cf_sample

# %%
def calc_correlations(coad_atac_data, brca_atac_data, coad_sample_ids, brca_sample_ids, brca_cf_array, coad_cf_array, bed_df):

    corr_mat1 = pd.DataFrame(
        np.nan,
        index= brca_flat_sample_ids,
        columns=brca_column_headers
    )
    
    corr_mat2 = pd.DataFrame(
        np.nan,
        index= brca_flat_sample_ids,
        columns= coad_column_headers
    )
    
    corr_mat3 = pd.DataFrame(
        np.nan,
        index= coad_flat_sample_ids,
        columns= coad_column_headers
    )

    corr_mat4 = pd.DataFrame(
        np.nan,
        index= coad_flat_sample_ids,
        columns= brca_column_headers
    )
    
    #brca brca    
    for i in range(brca_atac_data.shape[1]):
        atac_sample = brca_atac_data[:, i]
        for j in range(brca_cf_array.shape[1]):
            cf_sample = brca_cf_array[:, j]
            
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
            
            corr_mat1.loc[brca_flat_sample_ids[i], brca_column_headers[j]] = corr

    #brca coad    
    for i in range(brca_atac_data.shape[1]):
        atac_sample = brca_atac_data[:, i]
        for j in range(coad_cf_array.shape[1]):
            cf_sample = coad_cf_array[:, j]

            # # Skip correlation calculation if one of the columns has NaN, or is constant
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
            
            corr_mat2.loc[brca_flat_sample_ids[i], coad_column_headers[j]] = corr

    #Coad Coad    
    for i in range(coad_atac_data.shape[1]):
        atac_sample = coad_atac_data[:, i]
        for j in range(coad_cf_array.shape[1]):
            cf_sample = coad_cf_array[:, j]

            # # Skip correlation calculation if one of the columns has NaN, or is constant
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
            
            corr_mat3.loc[coad_flat_sample_ids[i], coad_column_headers[j]] = corr
    
    #Coad brca    
    for i in range(coad_atac_data.shape[1]):
        atac_sample = coad_atac_data[:, i]
        for j in range(brca_cf_array.shape[1]):
            cf_sample = brca_cf_array[:, j]
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
            
            corr_mat4.loc[coad_flat_sample_ids[i], brca_column_headers[j]] = corr
    

        
    #Save the correlation matrix for a cancer type
    corr_mat1.to_csv('correlation_results/all_samples/no_outliers/high_covxvar/covxvar_brca_atac_brca_cfDNA.csv', index= True) 
    corr_mat2.to_csv('correlation_results/all_samples/no_outliers/high_covxvar/covxvar_brca_atac_coad_cfDNA.csv', index= True) 
    corr_mat3.to_csv('correlation_results/all_samples/no_outliers/high_covxvar/covxvar_coad_atac_coad_cfDNA.csv', index= True) 
    corr_mat4.to_csv('correlation_results/all_samples/no_outliers/high_covxvar/covxvar_coad_atac_brca_cfDNA.csv', index= True) 

# %%
coad_atac_coad_cfDNA = calc_correlations(coad_atac_data, brca_atac_data, coad_sample_ids, brca_sample_ids, brca_cf_array, coad_cf_array, bed_df)

# %%


