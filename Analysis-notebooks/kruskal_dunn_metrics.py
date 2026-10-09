'''
Kruskal-Wallis test with Dunn's post-hoc when significant

Usage:
    python3 kruskal_dunn_metrics.py
'''

import os
import re
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from scikit_posthocs import posthoc_dunn

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

out_dir = 'plots'
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv('master_mag_list_final.csv', sep=',')
df = df[df['HQ-MAG'] == 1]
df = df[~df['Assembly_Name'].isin(['R006bc2010.E06', 'R006bc2010.F06'])]

# kit short name
def kit(name):
    if 'control' in name:
        return 'PCR-free'
    if 'Takara PrimeSTAR LongSeq' in name:
        return 'Takara LongSeq'
    for k in ['NEB Q5 XT', 'Kura Decodifi']:
        if k in name:
            return k
    return 'other'

def reaction_volume(name):
    m = re.search(r'(\d+(?:\.\d+)?)x', name)
    return f'{m.group(1)}x' if m else None

df['kit'] = df['name_short_3'].apply(kit)
df['reaction_volume'] = df['name_short_3'].apply(reaction_volume)

metrics = [
    'Completeness', 'Contamination', 'cMAG', 'cir_cMAG', 'Total_Contigs', 
    '16S_rRNA', '23S_rRNA', '5S_rRNA', 'total_rRNA_genes',
    'rRNA_operons', 'tRNA_genes', 'n_viruses',
    'genome_fraction_pct', 'GC_Content', 'N50', 'NGA50', 'misassemblies',
    'relocations', 'translocations', 'inversions',
    'mismatches_per_100kb', 'indels_per_100kb',
]
for col in metrics:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# labels
metric_labels = {
    'Completeness': 'Completeness (%)',
    'Contamination': 'Contamination (%)',
    'Total_Contigs': 'Total contigs',
    'cMAG': 'Complete MAG',
    'cir_cMAG': 'Circular complete MAG',
    '16S_rRNA': '16S rRNA',
    '23S_rRNA': '23S rRNA',
    '5S_rRNA': '5S rRNA',
    'total_rRNA_genes': 'rRNA genes',
    'rRNA_operons': 'rRNA Operons',
    'tRNA_genes': 'tRNA Genes',
    'n_viruses': 'Number of viruses',
    'genome_fraction_pct': 'Genome fraction',
    'GC_Content': 'GC content (%)',
    'N50': 'N50',
    'NGA50': 'NGA50',
    'misassemblies': 'Misassemblies',
    'relocations': 'Relocations',
    'translocations': 'Translocations',
    'inversions': 'Inversions',
    'mismatches_per_100kb': 'Mismatches per 100kb',
    'indels_per_100kb': 'Indels per 100kb',
}

# y-axis labels
metric_ylabels = {
    'Completeness': 'Completeness (%)',
    'Contamination': 'Contamination (%)',
    'Total_Contigs': 'Count',
    'cMAG': 'cMAG',
    'cir_cMAG': 'Circular cMAG',
    '16S_rRNA': '16S rRNA copies',
    '23S_rRNA': '23S rRNA copies',
    '5S_rRNA': '5S rRNA copies',
    'total_rRNA_genes': 'Total rRNA genes',
    'rRNA_operons': 'rRNA operons',
    'tRNA_genes': 'tRNA genes',
    'n_viruses': 'Number of viruses',
    'genome_fraction_pct': 'Fraction (%)',
    'GC_Content': 'GC content (%)',
    'N50': 'Length (bp)',
    'NGA50': 'Length (bp)',
    'misassemblies': 'Misassemblies (count)',
    'relocations': 'Relocations (count)',
    'translocations': 'Translocations (count)',
    'inversions': 'Inversions (count)',
    'mismatches_per_100kb': 'Mismatches per 100kb',
    'indels_per_100kb': 'Indel rate/100kb',
}

conditions = {
    'dna_input_ng': sorted(df['dna_input_ng'].dropna().unique()),
    'kit': ['PCR-free', 'Takara LongSeq', 'NEB Q5 XT', 'Kura Decodifi'],
    'reaction_volume': sorted(df['reaction_volume'].dropna().unique(), key=lambda v: float(v.rstrip('x'))),
}

# plot titles
condition_labels = {
    'dna_input_ng': 'DNA input',
    'kit': 'Kit',
    'reaction_volume': 'Reaction volume',
}

def stars(pval):
    if pval < 0.001:
        return '***'
    if pval < 0.01:
        return '**'
    if pval < 0.05:
        return '*'
    return ''

def pval_number_format(pval):
    base = '0.00E+00' if pval < 1e-4 else '0.0000'
    suffix = stars(pval)
    return f'{base}"{suffix}"' if suffix else base

def kw_stat_line(stat, pval):
    p_str = f'{pval:.2e}' if pval < 1e-4 else f'{pval:.4f}'
    return f'H = {stat:.2f}, p = {p_str}{stars(pval)}'

posthoc_results = {cond: [] for cond in conditions}
table_rows = []
pval_formats = {}
kw_stats = {}

for col in metrics:
    row = {'metric': col}
    for cond, order in conditions.items():
        groups = [df[df[cond] == g][col].dropna().values for g in order]
        groups = [g for g in groups if len(g)]
        if len(groups) < 2:
            row[(cond, 'H_statistic')] = float('nan')
            row[(cond, 'p_value')] = float('nan')
            continue

        stat, pval = kruskal(*groups)
        row[(cond, 'H_statistic')] = stat
        row[(cond, 'p_value')] = pval
        pval_formats[(col, cond)] = pval_number_format(pval)
        kw_stats[(col, cond)] = (stat, pval)

        if pval < 0.05:
            sub = df[[cond, col]].dropna()
            dunn = posthoc_dunn(sub, val_col=col, group_col=cond, p_adjust='holm')
            present = [g for g in order if g in dunn.index]
            for i, a in enumerate(present):
                for b in present[i + 1:]:
                    posthoc_results[cond].append({
                        'metric': col, f'{cond}_a': a, f'{cond}_b': b,
                        'p_value_adj': dunn.loc[a, b],
                        'sig': stars(dunn.loc[a, b]),
                    })
    table_rows.append(row)

metric_names = [r['metric'] for r in table_rows]
col_tuples = [(cond, stat) for cond in conditions for stat in ('H_statistic', 'p_value')]
results_table = pd.DataFrame([{k: v for k, v in r.items() if k != 'metric'} for r in table_rows])
results_table = results_table[col_tuples]
results_table.columns = pd.MultiIndex.from_tuples(results_table.columns)
results_table.index = metric_names
results_table.index.name = 'metric'

# out_path = f'{out_dir}/kruskal_wallis_by_condition.xlsx'
# with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
#     results_table.to_excel(writer, sheet_name='Kruskal-Wallis', index=True)
#     ws = writer.sheets['Kruskal-Wallis']
#     header_rows = 3  # condition row + stat row + index-name row
#     cond_list = list(conditions)
#     for r, col in enumerate(metric_names, start=header_rows + 1):
#         for cond in cond_list:
#             key = (col, cond)
#             if key not in pval_formats:
#                 continue
#             p_col = 2 + cond_list.index(cond) * 2 + 1  # 1-indexed: A=metric, then H,p pairs
#             ws.cell(row=r, column=p_col).number_format = pval_formats[key]

#     for cond, rows in posthoc_results.items():
#         if not rows:
#             continue
#         pd.DataFrame(rows).to_excel(writer, sheet_name=f'Dunn_{cond}'[:31], index=False)

kit_colors = {
    'PCR-free': 'grey',
    'Takara LongSeq': 'mediumpurple',
    'NEB Q5 XT': 'orange',
    'Kura Decodifi': 'skyblue',
}

def group_colors(cond, labels):
    if cond == 'kit':
        return [kit_colors.get(g, '#888888') for g in labels]
    cmap = plt.get_cmap('tab10')
    return [cmap(i % 10) for i in range(len(labels))]


# multi-panel figure
def draw_boxplot_panel(ax, cond, metric, order, stat, pval, panel_label=None):
    data_groups = [(g, df[df[cond] == g][metric].dropna().values) for g in order]
    data_groups = [(g, d) for g, d in data_groups if len(d)]
    labels = [g for g, d in data_groups]
    data = [d for g, d in data_groups]
    colors = group_colors(cond, labels)

    # boxplot
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=True)
    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.75)
    for median in bp['medians']:
        median.set_color('black')

    rng = np.random.default_rng(42)
    for i, d in enumerate(data, start=1):
        jitter = rng.uniform(-0.1, 0.1, size=len(d))
        ax.scatter(np.full(len(d), i) + jitter, d, color='black', s=8, alpha=0.6, zorder=3)

    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
    ax.tick_params(axis='y', labelsize=9)
    ax.set_ylabel(metric_ylabels[metric], fontsize=11)
    title = f'{panel_label}  {metric_labels[metric]}' if panel_label else metric_labels[metric]
    ax.set_title(title, fontsize=12, fontweight='bold', loc='left', pad=18)
    ax.text(0.5, 1.02, kw_stat_line(stat, pval), transform=ax.transAxes,
            ha='center', va='bottom', fontsize=9)

    label_idx = {g: i + 1 for i, g in enumerate(labels)}
    sig_pairs = [
        r for r in posthoc_results[cond]
        if r['metric'] == metric and r['p_value_adj'] < 0.05
        and r[f'{cond}_a'] in label_idx and r[f'{cond}_b'] in label_idx
    ]

    sig_pairs.sort(key=lambda r: abs(label_idx[r[f'{cond}_b']] - label_idx[r[f'{cond}_a']]))
    y_max = max(d.max() for d in data)
    y_min = min(d.min() for d in data)
    step = (y_max - y_min) * 0.08 or (abs(y_max) * 0.08 if y_max else 1)
    y_cursor = y_max + step * 1.5

    for r in sig_pairs:
        x1, x2 = label_idx[r[f'{cond}_a']], label_idx[r[f'{cond}_b']]
        ax.plot([x1, x1, x2, x2], [y_cursor, y_cursor + step * 0.3, y_cursor + step * 0.3, y_cursor],
                lw=1.0, color='black')
        ax.text((x1 + x2) / 2, y_cursor + step * 0.3, r['sig'], ha='center', va='bottom', fontsize=10)
        y_cursor += step * 1.3

    ax.set_ylim(top=y_cursor + step)


def build_composite_figure(cond):
    sig_metric_set = {r['metric'] for r in posthoc_results[cond]}
    sig_metrics = [m for m in metrics if m in sig_metric_set]
    if not sig_metrics:
        return None

    ncols = min(5, int(np.ceil(np.sqrt(len(sig_metrics)))))
    nrows = -(-len(sig_metrics) // ncols)  # ceil division
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4.5 * nrows))
    axes = np.atleast_1d(axes).flatten()

    panel_letters = string.ascii_uppercase
    for i, metric in enumerate(sig_metrics):
        stat, pval = kw_stats[(metric, cond)]
        draw_boxplot_panel(axes[i], cond, metric, conditions[cond], stat, pval,
                            panel_label=f'({panel_letters[i % 26]})')

    for ax in axes[len(sig_metrics):]:
        ax.axis('off')

    fig.text(0.5, 0.955, '* p<0.05    ** p<0.01    *** p<0.001', ha='center', fontsize=11, style='italic')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = f'{out_dir}/kw_{cond}.png'
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out_path

# Build figure
for cond in conditions:
    composite_path = build_composite_figure(cond)
    if composite_path:
        print(f'>> [Dunn] {cond}: composite figure written to {composite_path}')

print("Done")
