# inspect_h5ad.py

import anndata as ad

# Load the .h5ad file in backed mode so it doesn’t overload memory
adata = ad.read_h5ad("data/scrna/scrna_ann.h5ad", backed='r')

# Print basic info
print(f"Shape: {adata.shape}")              # Number of cells × genes
print(f"obs columns: {adata.obs.columns}")  # Metadata per cell
print(f"var columns: {adata.var.columns}")  # Metadata per gene
