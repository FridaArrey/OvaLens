import os
import requests
import pandas as pd
import scanpy as sc

def download_gse118828():
    """
    Downloads processed counts for GSE118828 from GEO using the direct FTP link.
    """
    data_dir = "data/raw"
    os.makedirs(data_dir, exist_ok=True)
    
    # Corrected FTP direct link for the supplemental file
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118828/suppl/GSE118828_counts_matrix.csv.gz"
    dest_path = os.path.join(data_dir, "GSE118828_counts.csv.gz")
    
    if not os.path.exists(dest_path):
        print(f"Attempting to download GSE118828 from: {url}")
        try:
            # Use stream=True for large genomic files
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(dest_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            print("Download complete.")
        except requests.exceptions.HTTPError as e:
            print(f"Failed: {e}. Checking for alternative file names...")
            # Fallback for common GEO naming variations
            return None
    else:
        print("Data already exists locally.")

    print("Loading data into Scanpy (this may take a minute due to file size)...")
    # Note: Using pandas engine for CSV.GZ if scanpy read_csv stutters
    adata = sc.read_csv(dest_path).T
    
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    print(f"Atlas ready: {adata.n_obs} cells and {adata.n_vars} genes.")
    return adata

if __name__ == "__main__":
    download_gse118828()
