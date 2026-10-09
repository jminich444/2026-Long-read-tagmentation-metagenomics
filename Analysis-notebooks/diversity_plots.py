"""
Fecal - Alpha and beta diversity plots

Usage:
    python3 diversity_plots.py
"""

from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.patches import Ellipse

excel_file = "diversity_tables_fecal.xlsx"
kingdoms_in = ["bacteria", "virus"]
outdir = Path("plots")

# alpha boxplot
alpha_metric = "sobs"
alpha_label = "Richness"
alpha_sheet = "alpha"
stats_sheet = "kruskal_wallis"
boxplot_group = "data_type"
boxplot_group_remap = {"EED pool": "Complex-fecal"}
box_width = 0.75
box_slot_in = 1.3
box_height_in = 5.5
box_palette = ["#A69CB5", "#8DA0A8", "#C9B896", "#B0B0B0", "#9BB68C"]
kingdom_colours = {"Bacteria": "#9BB68C", "Virus": "#B8A6CF"}

# aitchison pcoa
pcoa_sheet = "pcoa_aitchison"
variance_sheet = "pcoa_variance"
metadata_sheet = "metadata"
metadata_key = "sample_fastq"
pcoa_colour_by = "host_subject_id"
pcoa_colour_title = "Host subject ID"
host_remap = {"D5": "603", "H6": "798", "EED": "Complex-fecal"}
pcoa_shape_by = "libprep_mm"
pcoa_shape_title = "Library prep"
highlight_col = "sample_name"
highlight_values = ["D5", "H6"]
fmt = "png"
dpi = 300
font_size = 13
pcoa_font_size = 16
pcoa_figsize = (9.5, 7.5)
palette = "muted"
ellipse_scale = 2.4477468306808925

# read excel sheets
sheets = pd.read_excel(excel_file, sheet_name=None)

def read_sheet(table, kingdom):
    return sheets[f"{table}_{kingdom}"].copy()

# alpha-diversity boxplot
def alpha_boxplot_combined(kingdoms_list):
    meta = sheets[metadata_sheet].drop_duplicates(metadata_key)
    tables = []
    for k in kingdoms_list:
        alpha = read_sheet(alpha_sheet, k)
        if highlight_col and highlight_col not in alpha.columns and "sample_id" in alpha.columns:
            alpha = alpha.merge(meta[[metadata_key, highlight_col]], left_on="sample_id",
                                right_on=metadata_key, how="left")
        alpha["kingdom"] = k.capitalize()
        tables.append(alpha)
    alpha = pd.concat(tables, ignore_index=True).dropna(subset=[boxplot_group])
    alpha[boxplot_group] = alpha[boxplot_group].astype(str)
    if boxplot_group_remap:
        alpha[boxplot_group] = alpha[boxplot_group].replace(boxplot_group_remap)
    groups = sorted(alpha[boxplot_group].unique())
    kingdoms = [k.capitalize() for k in kingdoms_list]

    # split D5 / H6 samples
    names = [str(v) for v in highlight_values] if highlight_values else []
    is_hl = (alpha[highlight_col].astype(str).isin(names)
             if names and highlight_col in alpha.columns
             else pd.Series(False, index=alpha.index))

    fig, ax = plt.subplots(figsize=(box_slot_in * len(groups) + 2, box_height_in))
    sns.boxplot(data=alpha, x=boxplot_group, y=alpha_metric, order=groups,
                hue="kingdom", hue_order=kingdoms, palette=kingdom_colours,
                width=box_width, showfliers=False, ax=ax)
    sns.stripplot(data=alpha[~is_hl], x=boxplot_group, y=alpha_metric, order=groups,
                  hue="kingdom", hue_order=kingdoms, palette={k: "black" for k in kingdoms},
                  dodge=True, size=3, legend=False, ax=ax)
    handles, labels = ax.get_legend_handles_labels()

    if is_hl.any():
        hl_pool = ["#D55E00", "#56B4E9", "#009E73", "#CC79A7"]
        hl_colours = {n: hl_pool[i % len(hl_pool)] for i, n in enumerate(names)}
        n_k = len(kingdoms)
        offsets = {k: (i - (n_k - 1) / 2) * box_width / n_k for i, k in enumerate(kingdoms)}
        hl = alpha[is_hl]
        rng = np.random.default_rng(0)
        x = (hl[boxplot_group].map(groups.index) + hl["kingdom"].map(offsets)
             + rng.uniform(-0.05, 0.05, len(hl)))
        for n in names:
            m = (hl[highlight_col].astype(str) == n).to_numpy()
            ax.scatter(x[m], hl.loc[m, alpha_metric], s=36, color=hl_colours[n],
                       edgecolor="black", linewidth=0.6, zorder=5, label=n)
            handles.append(ax.collections[-1])
            labels.append(n)
    ax.legend(handles, labels, title="", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False,
              fontsize=font_size - 2, handletextpad=0.3)

    ax.set_xlabel("")
    ax.set_ylabel(alpha_label)
    ax.set_title(alpha_label, fontweight="bold")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    out = outdir / f"{alpha_metric}_boxplot_combined.{fmt}"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


# aitchison pcoa
def add_confidence_ellipse(ax, x, y, colour):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    if not np.all(np.isfinite(cov)) or cov[0, 0] <= 0 or cov[1, 1] <= 0:
        return
    pearson = np.clip(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]), -1, 1)
    ellipse = Ellipse((0, 0),
                      width=2 * np.sqrt(1 + pearson), height=2 * np.sqrt(1 - pearson),
                      facecolor=colour, alpha=0.15, edgecolor=colour, linewidth=1.8)
    ellipse.set_transform(
        transforms.Affine2D()
        .rotate_deg(45)
        .scale(np.sqrt(cov[0, 0]) * ellipse_scale, np.sqrt(cov[1, 1]) * ellipse_scale)
        .translate(x.mean(), y.mean())
        + ax.transData
    )
    ax.add_patch(ellipse)


def aitchison_pcoa(k):
    pc = read_sheet(pcoa_sheet, k)

    # metadata columns
    meta = sheets[metadata_sheet].drop_duplicates(metadata_key)
    if host_remap and "host_subject_id" in meta.columns:
        meta["host_subject_id"] = meta["host_subject_id"].astype(str).replace(host_remap)
    wanted = list(dict.fromkeys(c for c in (pcoa_colour_by, pcoa_shape_by, highlight_col) if c))
    pc = pc.drop(columns=[c for c in wanted if c in pc.columns])
    pc = pc.merge(meta[[metadata_key, *wanted]], left_on="sample_id", right_on=metadata_key, how="left")

    pc = pc.dropna(subset=[pcoa_colour_by]).copy()
    pc[pcoa_colour_by] = pc[pcoa_colour_by].astype(str)

    # marker shape
    shape_by = pcoa_shape_by or None
    if shape_by:
        pc[shape_by] = pc[shape_by].astype(str)
        if pc[shape_by].nunique() > 6:
            print(f"  [{k}] {shape_by}: {pc[shape_by].nunique()} values, too many for markers -> single shape")
            shape_by = None

    # highlighted samples
    is_hl = pd.Series(False, index=pc.index)
    if highlight_col and highlight_col in pc.columns and highlight_values:
        is_hl = pc[highlight_col].astype(str).isin([str(v) for v in highlight_values])
    bg, hl = pc[~is_hl], pc[is_hl]

    groups = sorted(bg[pcoa_colour_by].unique())
    pal = palette if len(groups) <= 10 else "husl" 
    colours = dict(zip(groups, sns.color_palette(pal, len(groups))))

    # % variance
    xlabel, ylabel = "PC1", "PC2"
    if variance_sheet and f"{variance_sheet}_{k}" in sheets:
        v = read_sheet(variance_sheet, k).set_index("metric").loc["aitchison"]
        xlabel = f"PC1 ({v['PC1_pct']:.1f}%)"
        ylabel = f"PC2 ({v['PC2_pct']:.1f}%)"

    # distinct marker per prep
    style_order = sorted(pc[shape_by].dropna().unique()) if shape_by else None
    marker_pool = ["o", "s", "^", "X", "D", "P", "v", "*"]
    markers = {c: marker_pool[i % len(marker_pool)]
               for i, c in enumerate(style_order)} if style_order else True

    fig, ax = plt.subplots(figsize=pcoa_figsize)
    for group in groups: 
        pts = bg[bg[pcoa_colour_by] == group]
        add_confidence_ellipse(ax, pts["PC1"].to_numpy(), pts["PC2"].to_numpy(), colours[group])
    sns.scatterplot(data=bg, x="PC1", y="PC2",
                    hue=pcoa_colour_by, hue_order=groups, palette=colours,
                    style=shape_by, style_order=style_order, markers=markers, s=20, ax=ax)

    # highlighted samples
    if len(hl):
        sns.scatterplot(data=hl, x="PC1", y="PC2", style=shape_by, style_order=style_order,
                        markers=markers, color="black", s=40, ax=ax, legend=False, zorder=5)
        for name, g in hl.groupby(highlight_col):
            ax.annotate(str(name), (g["PC1"].mean(), g["PC2"].mean()), fontsize=pcoa_font_size - 2,
                        fontweight="bold", zorder=6, xytext=(6, 6), textcoords="offset points")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"Aitchison Distance | {k.capitalize()}", fontweight="bold")
    leg = ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=pcoa_font_size - 2)
    headings = {pcoa_colour_by: pcoa_colour_title, shape_by: pcoa_shape_title}
    for txt in leg.get_texts(): 
        if txt.get_text() in headings:
            txt.set_text(headings[txt.get_text()])
            txt.set_fontweight("bold")
    fig.tight_layout()
    out = outdir / f"aitchison_pcoa_{k}.{fmt}"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    plt.rcParams["font.size"] = font_size
    outdir.mkdir(parents=True, exist_ok=True)
    alpha_boxplot_combined(kingdoms_in)
    with plt.rc_context({"font.size": pcoa_font_size}):
        for k in kingdoms_in:
            aitchison_pcoa(k)
