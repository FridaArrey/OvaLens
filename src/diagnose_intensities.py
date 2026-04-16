import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("data/external/scimap_processed_single_cell_TMA_dataset.h5ad")

print("--- 🔬 Protein Intensity Ranges ---")
# Check the raw numbers for your key markers
markers = ['ECadherin', 'aSMA', 'Vimentin']
for m in markers:
    vals = adata[:, m].X.flatten()
    print(f"{m} -> Min: {vals.min():.2f}, Max: {vals.max():.2f}, Mean: {vals.mean():.2f}")

print("\n--- 📝 Column Names in adata.obs ---")
print(adata.obs.columns.tolist())
