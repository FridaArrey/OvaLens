import pandas as pd
import numpy as np

def analyze():
    df = pd.read_csv("data/processed/validation/testa_epigenetic_ref.csv", index_col=0).transpose()
    df['origin_tag'] = df.index.to_series().apply(lambda x: 'OSE' if 'OSE' in str(x) else ('FT' if 'FI' in str(x) or 'FT' in str(x) else 'Unknown'))
    
    # Filter to only known origins
    known = df[df['origin_tag'].isin(['OSE', 'FT'])]
    
    # Calculate the mean methylation for OSE vs FT across all probes
    # (Simplified: looking for the most variable probes)
    ose_means = known[known['origin_tag'] == 'OSE'].iloc[:, :-1].mean()
    ft_means = known[known['origin_tag'] == 'FT'].iloc[:, :-1].mean()
    
    diff = (ose_means - ft_means).abs().sort_values(ascending=False)
    
    print("--- Top Discriminatory Epigenetic Probes ---")
    print(diff.head(10))
    print("\nThese probes represent the 'switches' that distinguish Ovarian from Fallopian origins.")

if __name__ == "__main__":
    analyze()
