import scanpy as sc
import os
import matplotlib.pyplot as plt
import pandas as pd

def generate_annotated_umap():
    processed_adata_path = "data/processed/GSE118828_atlas.h5ad"
    if not os.path.exists(processed_adata_path):
        print("Error: Processed atlas not found.")
        return

    print("Loading Atlas and re-calculating embeddings...")
    adata = sc.read_h5ad(processed_adata_path)
    
    # Ensure clustering exists
    sc.pp.pca(adata, n_comps=40)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.leiden(adata, resolution=0.5)
    sc.tl.umap(adata)

    # Biological Mapping
    cluster_names = {
        '0': 'Stromal (DCN+)',
        '1': 'Endothelial (VWF+)',
        '2': 'Smooth Muscle',
        '3': 'Pericytes',
        '4': 'Macrophages',
        '5': 'T-cells',
        '6': 'Epithelial',
        '7': 'Fibroblasts (COL1A1+)',
        '8': 'Ciliated/Fallopian',
        '9': 'Tumor A (LCN2+)',
        '10': 'Tumor B (KRT17+)'
    }
    
    adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_names)

    # Calculate Barrier Score (Fibroblasts / Tumor Cells)
    counts = adata.obs['cell_type'].value_counts()
    fibro_count = counts.get('Fibroblasts (COL1A1+)', 0)
    tumor_count = counts.get('Tumor A (LCN2+)', 0) + counts.get('Tumor B (KRT17+)', 0)
    
    barrier_score = fibro_count / tumor_count if tumor_count > 0 else 0

    print(f"\n--- OvaLens Metrics ---")
    print(f"Fibroblast Count: {fibro_count}")
    print(f"Tumor Cell Count: {tumor_count}")
    print(f"Calculated Barrier Score: {barrier_score:.2f}")
    print(f"-----------------------\n")

    # Plotting
    print("Generating Annotated UMAP...")
    sc.set_figure_params(dpi=150, fontsize=10)
    
    # We use the 'cell_type' column for color now
    sc.pl.umap(adata, 
               color='cell_type', 
               title=f'OvaLens TME (Barrier Score: {barrier_score:.2f})',
               legend_loc='right margin',
               show=False)
    
    os.makedirs('figures', exist_ok=True)
    plt.savefig('figures/ova_lens_annotated_umap.png', bbox_inches='tight')
    print("✅ Success! Annotated plot saved to: figures/ova_lens_annotated_umap.png")

if __name__ == "__main__":
    generate_annotated_umap()
