# %%
import numpy as np
import pandas as pd
import os
import glob
from scipy.stats import pearsonr
import anndata as ad
from scipy.sparse import csr_matrix
from scipy.stats import zscore

# %%
atac_dir = 'data/atac_files/celltype_mats'

# %%
scrna_dir = 'data/scrna/expr_mats_atac/atac_match'

# %%
def load_csv_matrices(folder_path):
    matrices = {}
    for file in os.listdir(folder_path):
        if file.endswith('.csv'):
            name = os.path.splitext(file)[0]
            df = pd.read_csv(os.path.join(folder_path, file), index_col=0)
            matrices[name] = df
            print(name)
    return matrices

# %%
atac_matrices = load_csv_matrices(atac_dir)

# %%
for scrna_file in os.listdir(scrna_dir):
    print(f"Now working with scRNA {scrna_file}")
    if not scrna_file.endswith('.csv'):
        continue

    scrna_name = os.path.splitext(scrna_file)[0]
    scrna_path = os.path.join(scrna_dir, scrna_file)
    scrna_df = pd.read_csv(scrna_path, index_col=0)
    scrna_df = scrna_df.dropna(axis=0, how='any')

    for atac_name, atac_df in atac_matrices.items():
        print(f"Now working with ATAC {atac_name}")

        # Find matching cell types (index labels)
        common_cell_types = atac_df.index.intersection(scrna_df.index)
        if len(common_cell_types) == 0:
            print(f"No common cell types between {atac_name} and {scrna_name}. Skipping.")
            continue

        # Subset both matrices
        atac_sub = atac_df.loc[common_cell_types]
        scrna_sub = scrna_df.loc[common_cell_types]

        # Combine and correlate
        combined = np.concatenate([atac_sub, scrna_sub], axis=0)
        corr = np.corrcoef(combined)

        # Construct correlation matrix (ATAC rows × scRNA columns)
        corr_matrix = pd.DataFrame(
            corr[:atac_sub.shape[0], atac_sub.shape[0]:],
            index=common_cell_types,
            columns=common_cell_types
        )

        diag_corrs = corr_matrix.values.diagonal()

      
        matched_df = pd.DataFrame({
            "Cell_Type": common_cell_types,
            "Pearson_Correlation": diag_corrs
        })

        matched_df.to_csv(f"correlation_results/heat/common/{atac_name}_{scrna_name}.csv", index=False)


print("All correlation matrices and mean summary saved.")

# %%



