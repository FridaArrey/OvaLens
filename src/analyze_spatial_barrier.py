import scimap as sm
import scanpy as sc
import pandas as pd

# 1. Load the data
adata = sc.read_h5ad("data/external/scimap_processed_single_cell_TMA_dataset.h5ad")

# 2. Gate the cells (Simplified Proteomics Gating)
# We define a "Barrier Cell" as high aSMA and a "Tumor Cell" as high ECadherin
# Adjusting thresholds based on common TMA intensities
sc.pp.log1p(adata)
adata.obs['cell_type'] = 'Stroma'
adata.obs.loc[adata[:, 'ECadherin'].X.flatten() > 1.0, 'cell_type'] = 'Tumor'
adata.obs.loc[adata[:, 'aSMA'].X.flatten() > 1.5, 'cell_type'] = 'Barrier_CAF'

print(f"Cell Composition:\n{adata.obs['cell_type'].value_counts()}")

# 3. Calculate Spatial Interactions (Neighborhood Analysis)
# This finds which cells are consistently next to each other
adata = sm.tl.spatial_interaction(adata, method='radius', radius=30, label='cell_type')

# 4. Calculate the Physical Barrier Ratio
counts = adata.obs['cell_type'].value_counts()
barrier_ratio = (counts['Barrier_CAF'] / counts['Tumor']) * 100
print(f"\n--- 🛰️ Physical Barrier Ratio: {barrier_ratio:.2f}% ---")

# 5. Save the interaction matrix for the report
interaction_df = adata.uns['spatial_interaction']
interaction_df.to_csv("results/spatial_interaction_matrix.csv")

print("✅ Spatial analysis complete. Results saved to results/spatial_interaction_matrix.csv")
