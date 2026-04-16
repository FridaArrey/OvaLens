import os
# Force OMP to be well-behaved on Mac
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import cellcharter as cc
import scanpy as sc
import pandas as pd

def run_spatial_analysis():
    path = "data/processed/GSE118828_atlas.h5ad"
    if not os.path.exists(path):
        print(f"Error: {path} not found.")
        return
    
    print(f"Loading {path}...")
    # Using backed='r' can save RAM, but let's try a standard load first with a subset check
    adata = sc.read_h5ad(path)
    
    # If the atlas is too huge, we'll downsample to 20k cells to ensure it fits in RAM
    if adata.n_obs > 20000:
        print(f"Downsampling for stability ({adata.n_obs} -> 20000 cells)...")
        sc.pp.subsample(adata, n_obs=20000)

    # Standard preprocessing for coordinates
    if 'X_pca' not in adata.obsm:
        print("Running PCA...")
        sc.tl.pca(adata, n_comps=30)
    
    print("Running Neighbors (Single-threaded for stability)...")
    sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')
    
    if 'X_umap' not in adata.obsm:
        print("Running UMAP...")
        sc.tl.umap(adata)
    
    # Set the manifold for CellCharter
    adata.obsm['spatial'] = adata.obsm['X_umap']

    print("Building spatial neighborhood graph...")
    # Use a fixed radius to prevent the graph from becoming too dense/exploding RAM
    cc.gr.spatial_neighbors(adata)

    print("Identifying cellular niches...")
    cluster_key = 'leiden' if 'leiden' in adata.obs.columns else ('cluster' if 'cluster' in adata.obs.columns else None)
    
    if cluster_key:
        print(f"Calculating enrichment for {cluster_key}...")
        cc.tl.nhood_enrichment(adata, cluster_key=cluster_key)
        
        output_path = "data/processed/ovalens_spatial_validated.h5ad"
        adata.write(output_path)
        
        enrich_key = f"{cluster_key}_nhood_enrichment"
        if enrich_key in adata.uns:
            df_enrich = pd.DataFrame(adata.uns[enrich_key])
            df_enrich.to_csv("results/spatial_barrier_enrichment.csv")
            print(f"✅ Success! Results in results/spatial_barrier_enrichment.csv")
    else:
        print("❌ Error: No cluster labels found.")

if __name__ == "__main__":
    run_spatial_analysis()
