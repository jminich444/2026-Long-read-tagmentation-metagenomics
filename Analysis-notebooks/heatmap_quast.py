'''
Heatmap for HQ-MAGs

Usage:
    python3 heatmap_quast.py
'''

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import seaborn as sns

out_dir = 'plots'
os.makedirs(out_dir, exist_ok=True)

# font settings
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['font.size'] = 10
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'

df = pd.read_csv('master_mag_list_final.csv', sep=',')

# keep only high-quality MAGs
df = df[df['HQ-MAG'] == 1]
df = df[~df['Assembly_Name'].isin(['R006bc2010.E06', 'R006bc2010.F06'])]
df = df.reset_index(drop=True)

_object_metric_cols = ['genome_fraction_pct', 'N50', 'NGA50', 'misassemblies',
                       'mismatches_per_100kb', 'indels_per_100kb',
                       'total_rRNA_genes', 'tRNA_genes', 'rRNA_operons', 'n_viruses',
                       '5S_rRNA', '16S_rRNA', '23S_rRNA', 'GC_Content',
                       'relocations', 'translocations', 'inversions']
df[_object_metric_cols] = df[_object_metric_cols].apply(pd.to_numeric)

# species
df['sp'] = df['species'].str.extract(r's__(.*)')[0]
sp_short = {
    'Enterococcus_B lactis': 'E. faecium',
    'Acinetobacter baumannii': 'A. baumannii',
    'Klebsiella pneumoniae': 'K. pneumoniae',
    'Pseudomonas aeruginosa': 'P. aeruginosa',
}
df['sp_label'] = df['sp'].map(sp_short)
sp_order = ['E. faecium', 'A. baumannii', 'K. pneumoniae', 'P. aeruginosa']
sp_colors = {'E. faecium': '#e66101', 'A. baumannii': '#5e3c99',
             'K. pneumoniae': '#008837', 'P. aeruginosa': '#d01c8b'}

# kit
def kit_from_name(name):
    if name.startswith('control (gtube)'):
        return 'PCR-free'
    elif name.startswith('Takara PrimeSTAR LongSeq'):
        return 'Takara LongSeq'
    elif name.startswith('Kura Decodifi'):
        return 'Kura Decodifi'
    elif name.startswith('NEB Q5 XT'):
        return 'NEB Q5 XT'
    else:
        raise ValueError(f'Unrecognized name_short_3 value: {name}')

df['kit_label'] = df['name_short_3'].apply(kit_from_name)
kit_order = ['PCR-free', 'Takara LongSeq', 'Kura Decodifi', 'NEB Q5 XT']
kit_colors = {'PCR-free': 'grey', 'Takara LongSeq': 'mediumpurple',
              'NEB Q5 XT': 'orange', 'Kura Decodifi': 'skyblue'}
CONTROL = 'PCR-free'

METRICS = {
    'Completeness':         dict(label='Completeness (%)',   fmt='{:.1f}', kind='mean'),
    'Contamination':        dict(label='Contamination (%)',  fmt='{:.2f}', kind='mean'),
    'cMAG':                 dict(label='cMAG',               fmt='{:.0f}', kind='binary'),
    'cir_cMAG':             dict(label='Circular cMAG',      fmt='{:.0f}', kind='binary'),
    'Total_Contigs':        dict(label='Total contigs',      fmt='{:.2f}', kind='mean'),
    'N50':                  dict(label='N50 (Mb)',           fmt='{:.2f}', kind='log', scale=1e-6),
    'NGA50':                dict(label='NGA50 (Mb)',         fmt='{:.2f}', kind='log', scale=1e-6),
    'genome_fraction_pct':  dict(label='Genome fraction (%)', fmt='{:.2f}', kind='mean'),
    'misassemblies':        dict(label='Misassemblies',      fmt='{:.2f}', kind='mean'),
    'mismatches_per_100kb': dict(label='Mismatches /100kb',  fmt='{:.2f}', kind='mean'),
    'indels_per_100kb':     dict(label='Indels /100kb',      fmt='{:.2f}', kind='mean'),
    'rRNA_operons':         dict(label='rRNA operons',       fmt='{:.2f}', kind='mean'),
    'tRNA_genes':           dict(label='tRNA genes',         fmt='{:.0f}', kind='mean'),
    'n_viruses':            dict(label='Viruses detected',   fmt='{:.2f}', kind='mean'),
    'total_rRNA_genes':     dict(label='Total rRNA genes', fmt='{:.2f}', kind='mean'),
    'GC_Content':           dict(label='GC content',         fmt='{:.2f}', kind='mean'),
    # supplementary-only metrics
    '16S_rRNA':             dict(label='16S rRNA',           fmt='{:.1f}', kind='mean'),
    '23S_rRNA':             dict(label='23S rRNA',           fmt='{:.1f}', kind='mean'),
    '5S_rRNA':              dict(label='5S rRNA',            fmt='{:.1f}', kind='mean'),
    'relocations':          dict(label='Relocations',        fmt='{:.1f}', kind='mean'),
    'translocations':       dict(label='Translocations',     fmt='{:.1f}', kind='mean'),
    'inversions':           dict(label='Inversions',         fmt='{:.1f}', kind='mean'),
}

# main figure
MAIN_GROUPS = {
    'Composition':    ['GC_Content'],
    'Circularity':    ['cMAG', 'cir_cMAG', 'Total_Contigs'],
    'Quality':        ['Completeness', 'Contamination'],
    'Gene content':   ['total_rRNA_genes', 'rRNA_operons', 'tRNA_genes', 'n_viruses'],
    'Accuracy':       ['N50', 'NGA50', 'genome_fraction_pct', 'misassemblies',
                       'mismatches_per_100kb', 'indels_per_100kb'],
}

SUPP_GROUPS = {
    'rRNA genes':     ['16S_rRNA', '23S_rRNA', '5S_rRNA'],
    'Misassembly\ntypes': ['relocations', 'translocations', 'inversions'],
}

def build_matrices(groups):
    """Return (colour matrix, text matrix) with rows = metrics, cols = (species, kit)."""
    metrics = [m for ms in groups.values() for m in ms]

    # values used for colouring
    tdf = df[['sp_label', 'kit_label'] + metrics].copy()
    for m in metrics:
        if METRICS[m]['kind'] == 'log':
            tdf[m] = np.log10(tdf[m])

    tmean = tdf.groupby(['sp_label', 'kit_label'])[metrics].mean()
    sd = tdf.groupby('sp_label')[metrics].std(ddof=0)
    cols = pd.MultiIndex.from_product([sp_order, kit_order], names=['sp', 'kit'])
    color = pd.DataFrame(np.nan, index=metrics, columns=cols)
    text = pd.DataFrame('', index=metrics, columns=cols)

    for sp in sp_order:
        if (sp, CONTROL) not in tmean.index:
            raise ValueError(f'No PCR-free control found for {sp}')
        for kit in kit_order:
            if (sp, kit) not in tmean.index:
                continue
            sub = df[(df['sp_label'] == sp) & (df['kit_label'] == kit)]
            for m in metrics:
                spec = METRICS[m]
                # colour
                if kit != CONTROL:
                    d = tmean.loc[(sp, kit), m] - tmean.loc[(sp, CONTROL), m]
                    s = sd.loc[sp, m]
                    color.loc[m, (sp, kit)] = 0.0 if (s == 0 or np.isnan(s)) else d / s
                # text
                if spec['kind'] == 'binary':
                    text.loc[m, (sp, kit)] = f'{int(sub[m].sum())}/{len(sub)}'
                else:
                    if spec['kind'] == 'log':
                        v = 10 ** tmean.loc[(sp, kit), m]
                    else:
                        v = tmean.loc[(sp, kit), m]
                    v = v * spec.get('scale', 1)
                    text.loc[m, (sp, kit)] = spec['fmt'].format(v)
    return color, text


# plot
def plot_heatmap(groups, fname, figsize, vlim=2.5, fs_scale=1.5):
    fs = lambda size: size * fs_scale
    annot_size = fs(8.5)
    color, text = build_matrices(groups)
    n_rows, n_cols = color.shape
    n_kits = len(kit_order)

    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(bottom=0.17)
    cax = fig.add_axes([0.20, 0.11, 0.30, 0.025])

    # colour cells
    sns.heatmap(color.values.astype(float), ax=ax, cbar_ax=cax,
                cmap='RdBu_r', center=0, vmin=-vlim, vmax=vlim,
                linewidths=0.6, linecolor='white', annot=False,
                cbar_kws={'orientation': 'horizontal'}, mask=color.isna().values)
    cax.set_xlabel('Difference from PCR-free control within species\n(red = higher, blue = lower)',
                   fontsize=fs(9.5))
    cax.tick_params(labelsize=fs(9))

    # grey control cells
    for j, (sp, kit) in enumerate(color.columns):
        for i in range(n_rows):
            if kit == CONTROL and text.iloc[i, j] != '':
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor='#e3e3e3',
                                       edgecolor='white', lw=0.6))

    # real values as cell text
    for i in range(n_rows):
        for j in range(n_cols):
            t = text.iloc[i, j]
            if not t:
                continue
            c = color.iloc[i, j]
            dark = (not np.isnan(c)) and abs(c) > 0.65 * vlim
            ax.text(j + 0.5, i + 0.5, t, ha='center', va='center',
                    fontsize=annot_size, color='white' if dark else 'black',
                    fontweight='bold' if kit_order[j % n_kits] == CONTROL else 'normal')

    # species block separators
    for k in range(1, len(sp_order)):
        ax.axvline(k * n_kits, color='black', lw=1.8)
    
    # metric group separators + group labels on the right
    r = 0
    for gname, ms in groups.items():
        if r > 0:
            ax.axhline(r, color='black', lw=1.8)
        r += len(ms)

    # axes labels / ticks
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels([METRICS[m]['label'] for m in color.index], rotation=0, fontsize=fs(10))
    ax.tick_params(axis='y', length=0)
    ax.set_xticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    trans = ax.get_xaxis_transform()
    for j, (_, kit) in enumerate(color.columns):
        ax.add_patch(Rectangle((j + 0.05, 1.01), 0.9, 0.03, transform=trans,
                               facecolor=kit_colors[kit], edgecolor='none', clip_on=False))
    # species headers
    for k, sp in enumerate(sp_order):
        x0, x1 = k * n_kits, (k + 1) * n_kits
        ax.add_patch(Rectangle((x0 + 0.05, 1.05), n_kits - 0.1, 0.045, transform=trans,
                               facecolor=sp_colors[sp], edgecolor='none', clip_on=False))
        ax.text((x0 + x1) / 2, 1.105, sp, transform=trans, ha='center', va='bottom',
                fontsize=fs(11), fontstyle='italic', fontweight='bold')

    # legend
    handles = [Patch(facecolor=kit_colors[k], label=k) for k in kit_order]
    fig.legend(handles=handles, title='Kit', loc='lower left', bbox_to_anchor=(0.58, 0.035),
               ncol=2, frameon=False, fontsize=fs(9.5),
               title_fontproperties={'size': fs(10), 'weight': 'bold'})

    fig.savefig(f'{out_dir}/{fname}', dpi=300, facecolor='white', bbox_inches='tight')
    plt.close(fig)

plot_heatmap(MAIN_GROUPS, 'master_heatmap_new_test.png', figsize=(15, 8.5))

print("Done")
