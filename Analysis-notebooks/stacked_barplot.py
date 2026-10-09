'''
Stacked barplot for HQ-MAGs

Usage:
    python3 stacked_barplot.py
'''

import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
# import seaborn as sns
import os

out_dir = 'plots'
os.makedirs(out_dir, exist_ok=True)

# font settings
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['font.size'] = 13
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'

df = pd.read_csv('master_mag_list_final.csv', sep=',')
# keep only high-quality MAGs
df = df[df['HQ-MAG'] == 1]
df = df[~df['Assembly_Name'].isin(['R006bc2010.E06', 'R006bc2010.F06'])]
df = df.reset_index(drop=True)

# metrics to display
metrics = ['HQ-MAG', 'cMAG', 'cir_cMAG']
df['name_short_3'] = df['name_short_3'].str.replace('control', 'PCR-free', regex=False)
df['name_short_3'] = df['name_short_3'].str.replace('Takara PrimeSTAR LongSeq', 'Takara LongSeq', regex=False)
grouped = df.groupby('name_short_3')[metrics].sum()
kit_order = ['PCR-free', 'Takara', 'NEB', 'Kura']

def _prefix_rank(name):
    for i, prefix in enumerate(kit_order):
        if name.startswith(prefix):
            return i
    raise ValueError(f'Unrecognized name_short_3 prefix: {name}')

def _dilution_rank(name):
    m = re.search(r'(\d*\.?\d+)x', name)
    return -float(m.group(1)) if m else 0

def _amount_rank(name):
    m = re.search(r'(\d+)\s*ng', name)
    return -float(m.group(1)) if m else 0

grouped = grouped.loc[sorted(grouped.index, key=lambda n: (_prefix_rank(n), _dilution_rank(n), _amount_rank(n)))]

# plot
fig, axes = plt.subplots(1, 1, figsize=(11, 6))

plot_data = pd.DataFrame({
    'Circular cMAG': grouped['cir_cMAG'],
    'HQ-MAG (non cMAG)': grouped['HQ-MAG'] - grouped['cMAG'],
})

colors = ['#55A868', '#4C72B0']
bottom = pd.Series(0, index=grouped.index)
for label, color in zip(plot_data.columns, colors):
    axes.bar(grouped.index, plot_data[label], bottom=bottom, label=label, color=color, alpha=0.85, width=0.7)
    bottom += plot_data[label]

axes.set_xticks(range(len(grouped.index)))
axes.set_xticklabels(grouped.index, rotation=45, ha='right', fontsize=15)
axes.set_xlim(-0.5, len(grouped.index) - 0.5)
axes.set_ylim(0, 4)
axes.set_yticks(range(5))
axes.tick_params(axis='y', labelsize=15)
axes.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=14)
axes.set_ylabel('Count', fontsize=15)
fig.tight_layout()
fig.savefig(f'{out_dir}/mags_stacked_plot.png', dpi=300, facecolor='white', bbox_inches='tight')

print("Done.")
