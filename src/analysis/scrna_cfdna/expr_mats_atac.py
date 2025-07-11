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
scrna_path = "data/scrna/tissues/ovary"


# %%
h5ad_files = [f for f in os.listdir(scrna_path) if f.endswith(".h5ad")]


# %%
atac_genes = pd.read_csv('atac_expres/atac_genes.csv')
atac_genes = atac_genes['Gene'].astype(str)

# %%
#anndata loop
for file_name in h5ad_files:
    print(f"\nStart processing {file_name}")
    start_time = time.time()

    adata = ad.read_h5ad(os.path.join(scrna_path, file_name))

    #Filteren op matchende genes
    ann_genes = adata.var['feature_name'].astype(str)
    matching_genes_mask = ann_genes.isin(atac_genes)
    filtered_adata = adata[:, matching_genes_mask]
    suba = filtered_adata[filtered_adata.obs["assay"] == "10x 3' v3", :]

    tissue_name = file_name.replace(".h5ad", "")

    rna_mat = suba.X
    print(rna_mat.shape)
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
    output_dir = "data/scrna/expr_mats_atac/atac_match"
    mean_df.to_csv(os.path.join(output_dir, f'{tissue_name}.csv'))


