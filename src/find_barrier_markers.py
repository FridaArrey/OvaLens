import scanpy as sc
import pandas as pd
import warnings

warnings.filterwarnings("ignore")

# Load the clustered atlas
print("Loading clustered atlas...")
adata = sc.read_h5ad("data/processed/GSE118828_atlas_clustered.h5ad")

# Run Differential Expression (Rank Genes Groups)
print("Ranking genes for Cluster 2...")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Extract the Top 10 genes for Cluster 2 (The Barrier)
result = sc.get.rank_genes_groups_df(adata, group="2")
top_10 = result.head(10)

print("\n--- Top Differential Markers for the Cluster 2 Barrier ---")
print(top_10[['names', 'logfoldchanges', 'pvals_adj']])

# Ensure directory exists and save
import os
os.makedirs("results/markers", exist_ok=True)
top_10.to_csv("results/markers/barrier_cluster_2_top_markers.csv", index=False)
print("\n✅ Markers saved to results/markers/barrier_cluster_2_top_markers.csv")
