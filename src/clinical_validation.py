import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
import numpy as np

# 1. Load the TCGA Matrix
tcga_df = pd.read_csv('data/clinical_validation/OV_TCGA.csv')

# 2. Strict Numeric Conversion & Cleanup
tcga_df['X0'] = pd.to_numeric(tcga_df['X0'], errors='coerce')
tcga_df['X1'] = pd.to_numeric(tcga_df['X1'], errors='coerce')
tcga_df['X75'] = pd.to_numeric(tcga_df['X75'], errors='coerce')
tcga_df = tcga_df.dropna(subset=['X0', 'X1', 'X75'])

# 3. Dynamic Thresholding
# If median doesn't work, we split based on the mean to ensure two groups
mean_risk = tcga_df['X75'].mean()
tcga_df['barrier_group'] = np.where(tcga_df['X75'] > mean_risk, 'High Barrier', 'Low Barrier')

# If one group is still empty, we do a 50/50 percentile split
if tcga_df['barrier_group'].nunique() < 2:
    tcga_df = tcga_df.sort_values('X75')
    mid_point = len(tcga_df) // 2
    tcga_df.iloc[:mid_point, tcga_df.columns.get_loc('barrier_group')] = 'Low Barrier'
    tcga_df.iloc[mid_point:, tcga_df.columns.get_loc('barrier_group')] = 'High Barrier'

print(f"✅ Data processed. Groups created: {tcga_df['barrier_group'].value_counts().to_dict()}")

# 4. Survival Analysis
kmf = KaplanMeierFitter()
plt.figure(figsize=(10, 6))

for group, color in zip(['High Barrier', 'Low Barrier'], ['red', 'blue']):
    mask = tcga_df['barrier_group'] == group
    group_data = tcga_df[mask]
    if not group_data.empty:
        kmf.fit(group_data['X0'], group_data['X1'], label=group)
        kmf.plot_survival_function(ax=plt.gca(), color=color, lw=2)

plt.title('OvaLens Clinical Validation: TCGA Survival by Barrier Risk', fontsize=14)
plt.xlabel('Standardized Survival Time', fontsize=12)
plt.ylabel('Overall Survival Probability', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)

plt.savefig('results/barrier_clinical_validation_KM.png', dpi=300)
print(f"✅ Final plot saved: results/barrier_clinical_validation_KM.png")
