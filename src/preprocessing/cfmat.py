import pandas as pd
import os

cfdna_matrix = {}

cfdna_folder = 'TSSClassification/results/04-15_10-29_all_coad'

# Loop over cfDNA sample files
for file in os.listdir(cfdna_folder):
    if not file.endswith(".csv"):
        continue

    sample_name = file.replace(".csv", "")
    df = pd.read_csv(os.path.join(cfdna_folder, file), index_col=0)

    # Assuming df has one column with values per genomic bin
    cfdna_matrix[sample_name] = df["cov_x"]


# Combine into a single DataFrame (rows = location, columns = sample)
cfdna_df = pd.DataFrame(cfdna_matrix)

# Optional: sort rows by location if desired
cfdna_df = cfdna_df.sort_index()

# Save to CSV
cfdna_df.to_csv("data/TCGA/coad_cf_mat.csv")