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

# %%
reference_df = pd.read_excel("rankjes/paper_data.xlsx", usecols=["cell_type", "tissue", "rank", "GC_code", "depth"])
reference_df.dropna(subset=["cell_type", "tissue"], inplace=True)
reference_df["normalized_cell_type"] = [name.rsplit("_", 1)[0] for name in reference_df["cell_type"]]
reference_df["normalized_tissue"] = reference_df["tissue"].str.lower().str.replace(r'[-\s]+', '_', regex=True)

# %%
h5ad_files = [f for f in os.listdir(scrna_path) if f.endswith(".h5ad")]

# %%
def strip_suffix(cell_type):
        for suffix in suffixes_to_remove:
            if cell_type.endswith(suffix):
                return cell_type[: -len(suffix)]
        return cell_type

# %%
filtered_rows = []


for file_name in h5ad_files:

    full_path = os.path.join(scrna_path, file_name)
    adata = ad.read_h5ad(full_path)

    adata.obs['cell_type'] = [name.replace(", ", "_").replace(" ", "_").replace("-", "_") for name in adata.obs['cell_type']]
    adata = adata[adata.obs["assay"] == "10x 3' v3", :]

    tissue_name = os.path.basename(file_name).replace(".h5ad", "")
    
    suffixes_to_remove = ["_Bone", "_Small", "_Large", "_Salivary", "_Lymph"]
    filtered_reference = reference_df[reference_df["normalized_tissue"] == tissue_name].copy()
    filtered_reference["normalized_cell_type"] = filtered_reference["normalized_cell_type"].apply(strip_suffix)

    valid_cell_types = set(filtered_reference["normalized_cell_type"].unique())
    current_cell_types = set(adata.obs['cell_type'].unique())
    matching_types = current_cell_types & valid_cell_types

    filtered_reference = filtered_reference[filtered_reference["normalized_cell_type"].isin(matching_types)]

    filtered_rows.append(filtered_reference)


matched_reference_df = pd.concat(filtered_rows, ignore_index=True)

# Filter for only 10x coverage 
filtered_depth_df = matched_reference_df[matched_reference_df["depth"] == "30X"].copy()

# %%
original_columns = filtered_depth_df.columns.tolist()
sorted_df = filtered_depth_df.sort_values(by=["GC_code", "rank"])

# %%
def assign_compressed_ranks(group):
    group = group.copy()
    group["rank"] = range(1, len(group) + 1)  # assign new ranks 1..n
    return group

# %%
recomputed_df = (
    sorted_df.groupby("GC_code", group_keys=False)
    .apply(assign_compressed_ranks)
)
recomputed_df = recomputed_df[original_columns]


# %%
recomputed_df.to_csv("rankjes/30x_matched_reference_ranks.csv", index=False)



