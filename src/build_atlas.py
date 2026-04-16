import os
import glob
import pandas as pd
import scanpy as sc

def build_hgsoc_atlas():
    data_dir = "data/raw"
    # Matching the specific CSV.GZ pattern from your ls output
    sample_files = glob.glob(os.path.join(data_dir, "GSM*.csv.gz"))
    
    adatas = []
    print(f"Found {len(sample_files)} sample files. Starting harmonization...")

    for file in sample_files:
        # Extracting a cleaner sample ID (e.g., GSM3348308_TUMOR)
        base_name = os.path.basename(file)
        parts = base_name.split('_')
        sample_id = f"{parts[0]}_{parts[2]}" if len(parts) > 2 else parts[0]
        
        print(f"Processing: {sample_id}")
        
        try:
            # These files are CSVs based on your ls output
            # We transpose (.T) because typically genes are rows in these tables
            tmp_adata = sc.read_csv(file).T
            tmp_adata.obs['sample_id'] = sample_id
            
            # Label tissue type based on filename for later analysis
            tmp_adata.obs['tissue_source'] = "Tumor" if "TUMOR" in base_name.upper() else "Environment"
            if "MET" in base_name.upper():
                tmp_adata.obs['status'] = "Metastatic"
            else:
                tmp_adata.obs['status'] = "Primary/Other"

            # Basic QC: filter out low-quality cells
            sc.pp.filter_cells(tmp_adata, min_genes=200)
            adatas.append(tmp_adata)
        except Exception as e:
            print(f"Could not process {file}: {e}")

    if not adatas:
        print("No valid data files found. Check file extensions or delimiters.")
        return None

    # Concatenate all samples into one master OvaLens Atlas
    print("Merging samples...")
    adata = sc.concat(adatas, join='outer')
    adata.obs_names_make_unique()
    
    # Save the processed atlas
    os.makedirs("data/processed", exist_ok=True)
    adata.write("data/processed/GSE118828_atlas.h5ad")
    print("\n✅ Success! Atlas saved to data/processed/GSE118828_atlas.h5ad")
    return adata

if __name__ == "__main__":
    build_hgsoc_atlas()
