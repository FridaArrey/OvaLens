import pyreadr
import os
import pandas as pd

def ingest_data():
    external_dir = "data/external"
    os.makedirs(external_dir, exist_ok=True)
    
    # We will focus on the large dataset first as it contains the labels
    dataset_path = "CellOfOriginDataset.rds"
    
    if os.path.exists(dataset_path):
        print(f"Reading {dataset_path} (137MB)... This may take a moment.")
        try:
            # Read the RDS file
            result = pyreadr.read_r(dataset_path)
            
            # Seurat objects in pyreadr often show up as multiple components
            # We want the metadata which contains the 'Origin' labels
            if None in result:
                df = result[None]
                output_path = f"{external_dir}/testa_metadata.csv"
                df.to_csv(output_path)
                print(f"✅ Success! Metadata saved to {output_path}")
            else:
                print("Checking object keys:", result.keys())
                # Try to find the metadata component
                for key in result.keys():
                    if 'meta' in key.lower():
                        result[key].to_csv(f"{external_dir}/testa_metadata.csv")
                        print(f"✅ Found and saved metadata from key: {key}")

        except Exception as e:
            print(f"❌ Error: {e}")
            print("\nIf this failed, it's because the file is a complex R Class.")
            print("PLAN B: I will provide the OriPrint Gene Signatures directly.")
    else:
        print("File not found.")

if __name__ == "__main__":
    ingest_data()
