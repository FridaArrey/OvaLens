import scanpy as sc
import matplotlib.pyplot as plt

# Load the newly clustered atlas
adata = sc.read_h5ad("data/processed/GSE118828_atlas_clustered.h5ad")

# Plot UMAP with Leiden clusters and the key Barrier Marker
sc.pl.umap(adata, color=['leiden', 'COL11A1', 'FAP'], 
           title=['Cell Clusters', 'COL11A1 (Barrier Marker)', 'FAP (CAF Marker)'],
           show=False, frameon=False)

plt.savefig("results/plots/barrier_validation_umap.png", bbox_inches='tight')
print("✅ Validation UMAP saved to results/plots/barrier_validation_umap.png")
