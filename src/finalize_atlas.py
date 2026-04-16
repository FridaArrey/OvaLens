import scanpy as sc
import pandas as pd
import numpy as np

# 1. Load the Atlas
print("Loading atlas...")
adata = sc.read_h5ad("data/processed/GSE118828_atlas.h5ad")
adata.obs_names_make_unique()

# 2. Crucial Cleaning: Remove empty cells and unexpressed genes
print("Cleaning data...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Handle potential NaNs from the merge
if pd.isna(adata.X).any():
    print("Fixing NaNs...")
    adata.X = np.nan_to_num(adata.X)

# 3. Normalization
print("Normalizing and Log-transforming...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 4. Highly Variable Genes (The fix for the Bin Edges error)
print("Finding Highly Variable Genes...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# 5. Dimensionality Reduction & Clustering
print("Running PCA, Neighbors, and UMAP...")
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5) 

# 6. Final Report
counts = adata.obs['leiden'].value_counts(normalize=True) * 100
print("\n--- 📊 Final Cluster Distribution ---")
print(counts)

# Save the Final Clustered version
adata.write("data/processed/GSE118828_atlas_clustered.h5ad")
print("\n✅ Success! Finalized Atlas saved.")
