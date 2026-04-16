import pandas as pd
import os

def validate():
    metadata_path = "data/external/testa_metadata.csv"
    if not os.path.exists(metadata_path):
        print("Metadata not found!")
        return

    # Read and Transpose because the columns are Cell IDs
    df_raw = pd.read_csv(metadata_path, index_col=0)
    df = df_raw.transpose()
    
    print("--- Adjusted Testa Lab Metadata ---")
    print(f"Total Cells in Reference: {len(df)}")
    
    # Extract Origin from the index (the cell names like '14-OSE-069')
    # We look for 'OSE' (Ovarian) and 'FI' or 'FT' (Fallopian)
    df['origin_tag'] = df.index.to_series().apply(lambda x: 'OSE' if 'OSE' in str(x) else ('FT' if 'FI' in str(x) or 'FT' in str(x) else 'Unknown'))
    
    origin_counts = df['origin_tag'].value_counts()
    print("\nOrigin Distribution in Landmark Study:")
    print(origin_counts)

    # Now we look for the cluster/cell type labels within the attributes
    # The 'rownames' column we saw earlier is now our columns. 
    # Let's look for 'seurat_clusters' or 'cell_type'
    if 'cell_type' in df.columns:
        # Calculate the ratio for OSE samples vs FT samples
        ose_data = df[df['origin_tag'] == 'OSE']
        ft_data = df[df['origin_tag'] == 'FT']
        
        print("\n--- Barrier Score Comparison ---")
        for label, sub_df in [("OSE (Ovarian)", ose_data), ("FT (Fallopian)", ft_data)]:
            counts = sub_df['cell_type'].value_counts()
            caf = counts.get('CAF', 0) + counts.get('Fibroblast', 0)
            tumor = counts.get('Cancer', 0) + counts.get('Epithelial', 0)
            if tumor > 0:
                print(f"{label} Ratio: {caf/tumor:.2f}")
    else:
        print("\nPotential Label Columns found:")
        print(df.columns.tolist()[:20])

if __name__ == "__main__":
    validate()
