import scanpy as sc

# Load the processed atlas
adata = sc.read_h5ad("data/processed/GSE118828_atlas.h5ad")

print("\n" + "="*50)
print("STRUCTURE INSPECTION: GSE118828_atlas")
print("="*50)
print(f"Total Cells: {adata.n_obs}")
print(f"Total Genes: {adata.n_vars}")
print("-" * 50)
print("Available Metadata (adata.obs):")
for col in adata.obs.columns:
    unique_count = adata.obs[col].nunique()
    sample_val = adata.obs[col].iloc[0]
    print(f"- {col:<20} | Unique: {unique_count:<5} | Example: {sample_val}")

print("-" * 50)
print("Available Embeddings (adata.obsm):")
print(list(adata.obsm.keys()))

# Check for our 'Barrier' genes
barrier_genes = ['COL11A1', 'FAP', 'THY1', 'POSTN']
found_genes = [g for g in barrier_genes if g in adata.var_names]
print("-" * 50)
print(f"Barrier Genes Found: {found_genes}")
print("="*50)
