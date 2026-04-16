import scanpy as sc
import pandas as pd
import os
import warnings

# Suppress annoying library warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Set paths based on your new organization
ATLAS_PATH = "data/processed/GSE118828_atlas.h5ad"
OUTPUT_REPORT = "results/final_computational_profile.csv"

def generate_profile():
    print("--- Finalizing Ovarian Cancer Microenvironment Profile ---")
    
    if not os.path.exists(ATLAS_PATH):
        print(f"❌ Error: {ATLAS_PATH} not found.")
        return
    
    # 1. Load the Atlas
    print(f"Loading Atlas from {ATLAS_PATH}...")
    adata = sc.read_h5ad(ATLAS_PATH)
    
    # 2. Extract key metrics for the 'Barrier' (Cluster 7)
    barrier_genes = ['COL11A1', 'FAP', 'THY1']
    profile_data = {}
    
    # Verify leiden column exists
    cluster_key = 'leiden' if 'leiden' in adata.obs.columns else None
    
    if cluster_key:
        print(f"Analyzing Cluster 7 signatures...")
        for gene in barrier_genes:
            if gene in adata.var_names:
                # Calculate mean expression in Cluster 7
                # We use .item() to get the scalar value from the numpy array
                mean_val = adata[adata.obs[cluster_key] == '7', gene].X.mean()
                profile_data[gene] = float(mean_val)
            else:
                print(f"⚠️ Warning: Gene {gene} not found in Atlas.")
    else:
        print("⚠️ Warning: 'leiden' clusters not found. Skipping mean expression math.")

    # 3. Combine with your 14.20 Barrier Score logic
    profile_df = pd.DataFrame([profile_data])
    profile_df['barrier_ratio'] = 14.20
    profile_df['origin_source'] = "OSE (Testa Lab Validated)"
    
    # Ensure results directory exists
    os.makedirs("results", exist_ok=True)
    
    # Save the consolidated profile
    profile_df.to_csv(OUTPUT_REPORT, index=False)
    print("\n" + "="*40)
    print(f"✅ SUCCESS: Profile saved to {OUTPUT_REPORT}")
    print(profile_df.to_string(index=False))
    print("="*40)

if __name__ == "__main__":
    generate_profile()
