import scanpy as sc
import os
import warnings

# Suppress annoying library warnings
warnings.filterwarnings("ignore")

# Check if the file is where we expect it
FILE_PATH = "data/external/scimap_processed_single_cell_TMA_dataset.h5ad"

def verify():
    if os.path.exists(FILE_PATH):
        print(f"--- 🛰️ OvaLens Spatial Data Verification ---")
        adata = sc.read_h5ad(FILE_PATH)
        print(f"✅ SUCCESS! Loaded TMA Atlas with {adata.n_obs} cells.")
        
        # Display the markers (Proteins) available for the Barrier analysis
        print("\nAvailable Protein Markers:")
        print(adata.var_names.tolist())
        
        # Check if the 'spatial_domain' column exists (the UTAG/Neighborhood result)
        if 'spatial_domain' in adata.obs.columns:
            print(f"\n✅ Spatial Domains found: {adata.obs['spatial_domain'].unique().tolist()}")
        else:
            print("\n⚠️ No 'spatial_domain' found. We will need to calculate these neighborhoods.")
            
    else:
        print(f"❌ File not found at: {os.path.abspath(FILE_PATH)}")
        print("Please ensure the .h5ad files were moved to data/external/")

if __name__ == "__main__":
    verify()
