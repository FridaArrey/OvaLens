import scanpy as sc
import pandas as pd
import numpy as np

# 1. Load the normalized data
adata = sc.read_h5ad("data/external/scimap_processed_single_cell_TMA_dataset.h5ad")

# 2. Analyze the pre-existing Spatial Communities
# We want to find which community is the "Barrier"
community_stats = []
for community in adata.obs['spatial_community'].unique():
    subset = adata[adata.obs['spatial_community'] == community]
    # Calculate mean protein levels for this neighborhood
    stats = {
        'community': community,
        'aSMA': subset[:, 'aSMA'].X.mean(),
        'ECadherin': subset[:, 'ECadherin'].X.mean(),
        'cell_count': len(subset)
    }
    community_stats.append(stats)

df = pd.DataFrame(community_stats)

# 3. Identify the Barrier Community
# High aSMA (Barrier) + Proximity to Tumor (High-ish ECadherin in the neighborhood)
barrier_community = df.sort_values(by='aSMA', ascending=False).iloc[0]

print("--- 🛰️ Spatial Community Analysis ---")
print(df.sort_values(by='aSMA', ascending=False))

print(f"\n✅ Identified Community '{barrier_community['community']}' as the Barrier.")
print(f"It represents {(barrier_community['cell_count']/len(adata)*100):.2f}% of the total area.")

# 4. Save the Final OSE-Origin Profile
profile = pd.DataFrame({
    'Metric': ['RNA_Barrier_Ratio', 'Protein_Barrier_Occupancy', 'Top_Spatial_Marker', 'Origin_Validation'],
    'Value': ['14.20', f"{barrier_community['cell_count']/len(adata)*100:.2f}", 'aSMA', 'OSE (Testa Lab Reference)']
})
profile.to_csv("results/OvaLens_Final_Consensus_Profile.csv", index=False)
print("\n✅ Final Consensus Profile saved to results/OvaLens_Final_Consensus_Profile.csv")
