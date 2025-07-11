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


# %%
scrna_path = "data/scrna/tissues"
# Bed voor matchen gene names
bed_df = pd.read_csv("TSSClassification/beds/dts_beds/tabula_match.bed", sep='\t', header=None)
print("Bed is ingeladen")


# %%
# cfDNA
ff_path = 'TSSClassification/results/05-01_16-42_brca_tabula'

# Cell types van paper ohalen
reference_df = pd.read_excel("rankjes/paper_data.xlsx", usecols=["cell_type", "tissue"])
reference_df.dropna(subset=["cell_type", "tissue"], inplace=True)
reference_df["normalized_cell_type"] = [name.rsplit("_", 1)[0] for name in reference_df["cell_type"]]
reference_df["normalized_tissue"] = reference_df["tissue"].str.lower().str.replace(r'[-\s]+', '_', regex=True)

# Load one cfDNA file to extract all fragment feature column names
first_file = glob.glob(os.path.join(ff_path, '*.csv'))[0]
fragment_features = [
    col for col in pd.read_csv(first_file, nrows=1).columns
    if col != 'Unnamed: 0'
]
print(f"Gevonden fragment features: {fragment_features}")

h5ad_files = [f for f in os.listdir(scrna_path) if f.endswith(".h5ad")]

def strip_suffix(cell_type):
        for suffix in suffixes_to_remove:
            if cell_type.endswith(suffix):
                return cell_type[: -len(suffix)]
        return cell_type


total_matched_cell_types = 0  # initialize counter

#anndata loop
for file_name in h5ad_files:
    print(f"\nStart processing {file_name}")
    start_time = time.time()

    adata = ad.read_h5ad(os.path.join(scrna_path, file_name))

    # Normalize cell type names
    adata.obs['cell_type'] = (adata.obs['cell_type'].str.replace(", ", "_").str.replace(" ", "_").str.replace("-", "_"))

    bed_genes = bed_df.iloc[:, 3]
    ann_genes = adata.var['feature_name']
    matching_genes = ann_genes.isin(bed_genes)
    filtered_adata = adata[:, matching_genes]
    suba = filtered_adata[filtered_adata.obs["assay"] == "10x 3' v3", :]

    tissue_name = file_name.replace(".h5ad", "")

    #filter reference and strip the suffixes
    suffixes_to_remove = ["_Bone", "_Small", "_Large", "_Salivary", "_Lymph"]
    filtered_reference = reference_df[reference_df["normalized_tissue"] == tissue_name].copy()
    filtered_reference["normalized_cell_type"] = filtered_reference["normalized_cell_type"].apply(strip_suffix)

    valid_cell_types = set(filtered_reference["normalized_cell_type"].unique())
    current_cell_types = set(suba.obs['cell_type'].unique())
    matching_types = current_cell_types & valid_cell_types

    suba = suba[suba.obs["cell_type"].isin(matching_types)]
    
    matched_count = suba.obs['cell_type'].nunique()
    total_matched_cell_types += matched_count
    print(f"Filtering completed: {matched_count} unique matched cell types")


    rna_mat = suba.X
    if issparse(rna_mat):
        rna_mat = rna_mat.tocsr()

    # Calculate mean per cell type
    cell_types = suba.obs["cell_type"].values
    unique_cell_types, inverse_idx = np.unique(cell_types, return_inverse=True)
    mean_rna_mat = np.zeros((len(unique_cell_types), rna_mat.shape[1]), dtype=np.float32)

    ##Saven hier 
    for i, cell_type_idx in enumerate(range(len(unique_cell_types))):
        mask = np.where(inverse_idx == cell_type_idx)[0]
        submat = rna_mat[mask]
        mean = submat.mean(axis=0)
        mean_rna_mat[i, :] = np.asarray(mean).ravel()

    rna_mat_z = zscore(mean_rna_mat, axis=1, ddof=0)

    # Feature loop
    for feature in fragment_features:
        print(f"Calculating correlation for feature: {feature}")

        file_names = [f for f in os.listdir(ff_path) if f.endswith(".csv")]
        cf_array = np.zeros((len(file_names), len(bed_df)), dtype=np.float32)

        for j, ff_file in enumerate(file_names):
            cf_file_path = os.path.join(ff_path, ff_file)
            ff_vals = pd.read_csv(cf_file_path, usecols=[feature], dtype=float)[feature].values
            cf_array[j, :] = ff_vals

        cf_array_z = zscore(cf_array, axis=1, ddof=0)
        corr_mat = np.dot(rna_mat_z, cf_array_z.T) / rna_mat_z.shape[1]

        output_dir = f'correlation_results/scrna/brca/{feature}'
        os.makedirs(output_dir, exist_ok=True)
        
        sample_ids = [f.split('.')[0] for f in file_names]
        corr_df = pd.DataFrame(corr_mat, index=unique_cell_types, columns=sample_ids)

     
        corr_df.to_csv(os.path.join(output_dir, f'{tissue_name}.csv'))

    print(f"Finished {tissue_name} in {time.time() - start_time:.1f} sec")


print(f"Total matched cell types across all tissues: {total_matched_cell_types}")

