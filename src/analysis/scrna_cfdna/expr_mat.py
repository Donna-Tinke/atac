# %%
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
scrna_path = "data/scrna/tissues/stomach"
# Bed voor matchen gene names
bed_df = pd.read_csv("TSSClassification/beds/dts_beds/tabula_match.bed", sep='\t', header=None)
print("Bed is ingeladen")

# %%
# # Cell types van paper ohalen
# reference_df = pd.read_excel("rankjes/paper_data.xlsx", usecols=["cell_type", "tissue"])
# reference_df.dropna(subset=["cell_type", "tissue"], inplace=True)
# reference_df["normalized_cell_type"] = [name.rsplit("_", 1)[0] for name in reference_df["cell_type"]]
# reference_df["normalized_tissue"] = reference_df["tissue"].str.lower().str.replace(r'[-\s]+', '_', regex=True)

# %%
h5ad_files = [f for f in os.listdir(scrna_path) if f.endswith(".h5ad")]

def strip_suffix(cell_type):
        for suffix in suffixes_to_remove:
            if cell_type.endswith(suffix):
                return cell_type[: -len(suffix)]
        return cell_type


total_matched_cell_types = 0  # initialize counter

# %%
#anndata loop
for file_name in h5ad_files:
    print(f"\nStart processing {file_name}")
    start_time = time.time()

    adata = ad.read_h5ad(os.path.join(scrna_path, file_name))

    # Normalize cell type names
    adata.obs['cell_type'] = (adata.obs['cell_type'].str.replace(", ", "_").str.replace(" ", "_").str.replace("-", "_"))

    #Filteren op matchende bed locaties
    bed_genes = bed_df.iloc[:, 3]
    ann_genes = adata.var['feature_name']
    matching_genes = ann_genes.isin(bed_genes)
    filtered_adata = adata[:, matching_genes]
    suba = filtered_adata[filtered_adata.obs["assay"] == "10x 3' v3", :]

    tissue_name = file_name.replace(".h5ad", "")

    #Strip the suffixes en filter op aleen matchende cell types
    suffixes_to_remove = ["_Bone", "_Small", "_Large", "_Salivary", "_Lymph"]
    # suba = reference_df[reference_df["normalized_tissue"] == tissue_name].copy()
    # filtered_reference["normalized_cell_type"] = filtered_reference["normalized_cell_type"].apply(strip_suffix)

    # valid_cell_types = set(filtered_reference["normalized_cell_type"].unique())
    # current_cell_types = set(suba.obs['cell_type'].unique())
    # matching_types = current_cell_types & valid_cell_types

    # suba = suba[suba.obs["cell_type"].isin(matching_types)]
    
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

    for i, cell_type_idx in enumerate(range(len(unique_cell_types))):
        mask = np.where(inverse_idx == cell_type_idx)[0]
        submat = rna_mat[mask]
        mean = submat.mean(axis=0)
        mean_rna_mat[i, :] = np.asarray(mean).ravel()

    #Convert to pd dataframe
    mean_df = pd.DataFrame(mean_rna_mat, index=unique_cell_types)

    # Expression matrixje saven
    output_dir = "data/scrna/expr_mats"
    mean_df.to_csv(os.path.join(output_dir, f'{tissue_name}.csv'))

# %%



