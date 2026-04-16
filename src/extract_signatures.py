import scanpy as sc
import os
import pandas as pd
import numpy as np

def extract_signatures():
    atlas_path = "data/processed/GSE118828_atlas.h5ad"
    if not os.path.exists(atlas_path):
        print("Error: Atlas file not found.")
        return

    print("Loading Atlas...")
    adata = sc.read_h5ad(atlas_path)
    
    print("Sanitizing Matrix...")
    X_data = pd.DataFrame(adata.X)
    X_data = X_data.fillna(0).apply(pd.to_numeric, errors='coerce').fillna(0)
    adata.X = X_data.values

    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_counts=1)
    
    print(f"Sanitized Shape: {adata.shape}")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    print("Selecting Highly Variable Genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=1500, flavor='seurat')
    adata = adata[:, adata.var.highly_variable].copy()

    print("Clustering populations (using igraph)...")
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    # Using the recommended 'igraph' flavor for speed
    sc.tl.leiden(adata, resolution=0.5, flavor='igraph', n_iterations=2)
    
    print("Ranking genes...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    os.makedirs("results", exist_ok=True)
    result_df = sc.get.rank_genes_groups_df(adata, group=None)
    result_df.to_csv("results/cluster_markers.csv")
    
    print("\n✅ Success! Markers extracted to results/cluster_markers.csv")
    # Identify the top marker for the first cluster
    top_gene = result_df[result_df['group']=='0']['names'].iloc[0]
    print(f"Top marker for Cluster 0: {top_gene}")

if __name__ == "__main__":
    extract_signatures()
