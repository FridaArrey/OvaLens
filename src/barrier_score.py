import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

def calculate_barrier_score(adata, tumor_key="Epithelial/Tumor", shield_key="Shield (Fibroblasts)"):
    """
    Calculates the 'Barrier Score' for OvaLens.
    Measures the density of Fibroblast 'Shields' surrounding Tumor clusters.
    """
    print(f"Calculating Barrier Score for {tumor_key} vs {shield_key}...")
    
    # Extract coordinates (assuming Visium or similar spatial format)
    if 'spatial' not in adata.obsm:
        print("Error: No spatial coordinates found in adata.obsm['spatial']")
        return 0.0
        
    coords = adata.obsm['spatial']
    
    # Identify indices for tumor and shield cells
    # Note: Labels should be calibrated using reference atlas GSE118828
    if 'cell_type' not in adata.obs:
        print("Error: 'cell_type' column missing in adata.obs")
        return 0.0

    tumor_mask = adata.obs['cell_type'] == tumor_key
    shield_mask = adata.obs['cell_type'] == shield_key
    
    if not (any(tumor_mask) and any(shield_mask)):
        print(f"Warning: Missing required cell types ({tumor_key} or {shield_key}) for calculation.")
        return 0.0

    # Fit Nearest Neighbors on the 'Shield' cells
    nn = NearestNeighbors(n_neighbors=5).fit(coords[shield_mask])
    
    # Find distance from each Tumor cell to the nearest 5 Shield cells
    distances, _ = nn.kneighbors(coords[tumor_mask])
    
    # The Barrier Score: Inverse of average distance 
    # (High score = Tumor is tightly surrounded by a fibroblast shield)
    avg_distances = np.mean(distances, axis=1)
    barrier_score = 1 / (np.mean(avg_distances) + 1e-9)
    
    adata.uns['ova_lens_barrier_score'] = barrier_score
    print(f"Final Barrier Score: {barrier_score:.4f}")
    
    return barrier_score

if __name__ == "__main__":
    print("OvaLens Barrier Engine Initialized. Ready for spatial deconvolution.")
