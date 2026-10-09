'''
Stacked relative-abundance plot of the ESKAPE mock community

Usage:
    python3 eskape_abundance_plot.py
'''

import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

out_dir = 'plots'
os.makedirs(out_dir, exist_ok=True)

# font settings
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['font.size'] = 12
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'

df = pd.read_csv('eskape_abundance_kit.tsv', sep='\t')

# species name
sp_cols = [c for c in df.columns if c.startswith('d__')]
sp_short = {
    'Enterococcus_B lactis': 'E. faecium',
    'Staphylococcus aureus': 'S. aureus',
    'Klebsiella pneumoniae': 'K. pneumoniae',
    'Acinetobacter baumannii': 'A. baumannii',
    'Pseudomonas aeruginosa': 'P. aeruginosa',
}
rename = {c: sp_short[re.search(r's__([^|]+)', c).group(1)] for c in sp_cols}
df = df.rename(columns=rename)

# stacking order
sp_order = ['E. faecium', 'S. aureus', 'K. pneumoniae', 'A. baumannii', 'P. aeruginosa']
sp_colors = {'E. faecium': '#E39B6B', 'S. aureus': '#E3C567', 'A. baumannii': '#7EA6D8',
             'K. pneumoniae': '#6BAA75', 'P. aeruginosa': '#C9625F'}

df['name_short_3'] = df['name_short_3'].str.replace('control', 'PCR-free', regex=False)
df['name_short_3'] = df['name_short_3'].str.replace('Takara PrimeSTAR LongSeq', 'Takara LongSeq', regex=False)

abund = df.set_index('name_short_3')[sp_order]

# order conditions
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

abund = abund.loc[sorted(abund.index, key=lambda n: (_prefix_rank(n), _dilution_rank(n), _amount_rank(n)))]

# plot
fig, ax = plt.subplots(1, 1, figsize=(11, 6))

x = range(len(abund.index))
bottom = pd.Series(0.0, index=abund.index)
for sp in sp_order:
    ax.bar(x, abund[sp], bottom=bottom, label=sp, color=sp_colors[sp], width=0.7,
           edgecolor='white', linewidth=0.4)
    bottom += abund[sp]

ax.set_xticks(list(x))
ax.set_xticklabels(abund.index, rotation=45, ha='right', fontsize=13)
ax.set_xlim(-0.5, len(abund.index) - 0.5)
ax.set_ylim(0, 100)
ax.tick_params(axis='y', labelsize=13)
ax.set_ylabel('Relative abundance (%)', fontsize=13)
legend = ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=len(sp_order),
                   frameon=False, fontsize=12, handlelength=1.2, columnspacing=1.0)
for t in legend.get_texts():
    t.set_fontstyle('italic')
for side in ('top', 'right'):
    ax.spines[side].set_visible(False)

fig.tight_layout()
fig.savefig(f'{out_dir}/eskape_abundance_plot.png', dpi=300, facecolor='white', bbox_inches='tight')
print("Done.")
