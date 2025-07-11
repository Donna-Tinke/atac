# %%
import pandas as pd
import os
import anndata as ad

# %%
scrna_path = "data/scrna/tissues"

# %%
reference_df = pd.read_excel("rankjes/paper_data.xlsx", usecols=["cell_type", "tissue"])
reference_df.dropna(subset=["cell_type", "tissue"], inplace=True)
reference_df["normalized_cell_type"] = [name.rsplit("_", 1)[0] for name in reference_df["cell_type"]]
reference_df["normalized_tissue"] = reference_df["tissue"].str.lower().str.replace(r'[-\s]+', '_', regex=True)

# %%
h5ad_files = [f for f in os.listdir(scrna_path) if f.endswith(".h5ad")]

# %%
cell_name_rows = []

for file_name in h5ad_files:

    full_path = os.path.join(scrna_path, file_name)
    adata = ad.read_h5ad(full_path)
    adata.obs['cell_type'] = [name.replace(", ", "_").replace(" ", "_").replace("-", "_") for name in adata.obs['cell_type']] 
    adata = adata[adata.obs["assay"] == "10x 3' v3", :]

    tissue_name = os.path.basename(file_name).replace(".h5ad", "")

    # Strip suffixes from normalized_cell_type in filtered_reference
    suffixes_to_remove = ["_Bone", "_Small", "_Large", "_Salivary", "_Lymph"]

    def strip_suffix(cell_type):
        for suffix in suffixes_to_remove:
            if cell_type.endswith(suffix):
                return cell_type[: -len(suffix)]
        return cell_type

    filtered_reference = reference_df[reference_df["normalized_tissue"] == tissue_name].copy()
    filtered_reference["normalized_cell_type"] = filtered_reference["normalized_cell_type"].apply(strip_suffix)

    paper_types = sorted(filtered_reference["normalized_cell_type"].unique())
    anndata_types = sorted(adata.obs["cell_type"].unique())
                                         
    # Find matches
    matches = sorted(set(paper_types).intersection(set(anndata_types)))
    paper_only = sorted(set(paper_types) - set(matches))
    anndata_only = sorted(set(anndata_types) - set(matches))

    # Combine: matches first, then mismatches
    combined_paper = matches + paper_only + [None] * (len(anndata_only) - len(paper_only)) if len(anndata_only) > len(paper_only) else matches + paper_only
    combined_anndata = matches + anndata_only + [None] * (len(paper_only) - len(anndata_only)) if len(paper_only) > len(anndata_only) else matches + anndata_only

    # Pad to equal length
    max_len = max(len(combined_paper), len(combined_anndata))
    combined_paper += [None] * (max_len - len(combined_paper))
    combined_anndata += [None] * (max_len - len(combined_anndata))

    # Build tissue-specific DataFrame
    df_tissue = pd.DataFrame({
        f"{tissue_name}_paper_types": combined_paper,
        f"{tissue_name}_anndata_types": combined_anndata
    })

    cell_name_rows.append(df_tissue)


# Combine all tissue DataFrames horizontally
final_df = pd.concat(cell_name_rows, axis=1)

# Save to CSV
final_df.to_csv("rankjes/celltype_comparison2.csv", index=False)

# Optional preview
final_df.head()


# %%



