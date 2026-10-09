"""
Coverage-depth plots

Usage:
    python3 coverage_depth_plots.py
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess

PANEL_DIR = Path("plots")
PANEL_DIR.mkdir(parents=True, exist_ok=True)

DPI = 300
FIG_SIZE = (7, 9)
PANEL_W, PANEL_H = FIG_SIZE[0] / 3, FIG_SIZE[1] / 3
SCALE = 0.45
FONT = {"title": 4.5, "label": 4.5, "tick": 4, "suptitle": 6}
DEPTH_LABEL = "Normalized depth"

# kits from top to bottom
KIT_ORDER = ["Takara PrimeSTAR LongSeq", "NEB Q5 XT", "Kura Decodifi"]
KIT_COLORS = {"Takara PrimeSTAR LongSeq": "mediumpurple", "NEB Q5 XT": "orange",
              "Kura Decodifi": "skyblue"}
REPLICATES = ["", "(Q)"]
KIT_LABELS = {"Takara PrimeSTAR LongSeq": "Takara LongSeq"}
SPECIES_LABELS = {"Enterococcus_B lactis": "Enterococcus faecium"}

windows = pd.read_csv("all_MAGs.25kb.depth_GC.tsv", sep="\t")
summary = pd.read_csv("all_MAGs.coverage_summary.tsv", sep="\t")

# split Kit
for table in (windows, summary):
    table["Kit_Name"] = table["Kit"].map(lambda k: next(n for n in KIT_ORDER if k.startswith(n)))
    table["Replicate"] = table["Kit"].map(lambda k: "(Q)" if k.endswith("(Q)") else "")
    # label shown on plots
    table["Kit_Label"] = (table["Kit_Name"].replace(KIT_LABELS) + " " + table["Replicate"]).str.strip()

def save(fig, name):
    fig.savefig(PANEL_DIR / f"{name}.png", dpi=DPI)
    plt.close(fig)

# Position along each MAG vs normalized depth
def plot_position_depth(ax, df, color):
    # place contigs end to end
    contig_end = df.groupby("Contig", sort=False)["End"].max()
    offset = contig_end.cumsum() - contig_end 
    x_mb = ((df["Start"] + df["End"]) / 2 + df["Contig"].map(offset)) / 1e6

    ax.plot(x_mb, df["Normalized_Depth"], linewidth=0.8 * SCALE, alpha=0.85, color=color)
    ax.scatter(x_mb, df["Normalized_Depth"], s=6 * SCALE**2, alpha=0.8, zorder=3, color=color)
    ax.axhline(1, linestyle="--", linewidth=0.9 * SCALE, color="grey")    # MAG median
    for boundary in contig_end.cumsum().iloc[:-1]:                # contig boundaries
        ax.axvline(boundary / 1e6, linewidth=0.5 * SCALE, alpha=0.3, color="grey")


# GC % vs normalized depth
def plot_gc_depth(ax, df, color):
    gc, depth = df["GC_Percent"], df["Normalized_Depth"]
    rho, pvalue = spearmanr(gc, depth)
    ax.scatter(gc, depth, s=12 * SCALE**2, alpha=0.6, edgecolors="none", color=color)
    smooth = lowess(endog=depth, exog=gc, frac=0.5, it=1, return_sorted=True)
    ax.plot(smooth[:, 0], smooth[:, 1], linewidth=1.8 * SCALE, color="black")
    ax.axhline(1, linestyle="--", linewidth=0.9 * SCALE, color="grey")    # MAG median
    ax.text(0.97, 0.03, f"ρ = {rho:.2f}, P = {pvalue:.1e}\nn = {len(df)}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8 * SCALE,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5 * SCALE))

# species panel
def species_panel(species, kind):
    sp_df = windows[windows["species"] == species]
    fig, axes = plt.subplots(len(KIT_ORDER), len(REPLICATES), figsize=(PANEL_W, PANEL_H),
                             sharey=True, squeeze=False)
    for r, kit_name in enumerate(KIT_ORDER):
        for c, replicate in enumerate(REPLICATES):
            df = sp_df[(sp_df["Kit_Name"] == kit_name) & (sp_df["Replicate"] == replicate)]
            ax = axes[r, c]
            ax.set_title(f"{KIT_LABELS.get(kit_name, kit_name)} {replicate}".strip(),
                         fontsize=FONT["title"], pad=2)
            ax.grid(axis="y", alpha=0.2)
            ax.tick_params(labelsize=FONT["tick"], length=2, pad=1)
            if df.empty: 
                continue
            if kind == "position":
                plot_position_depth(ax, df, KIT_COLORS[kit_name])
            else:
                plot_gc_depth(ax, df, KIT_COLORS[kit_name])

    xlabel = "MAG position (Mb)" if kind == "position" else "Window GC content (%)"
    axes[0, 0].set_ylim(bottom=0)
    for ax in axes[-1, :]:
        ax.set_xlabel(xlabel, fontsize=FONT["label"], labelpad=1)
    # one y-axis label for the whole grid
    fig.supylabel(DEPTH_LABEL, fontsize=FONT["label"], x=0.08 / PANEL_W)
    fig.suptitle(SPECIES_LABELS.get(species, species), fontstyle="italic", fontweight="bold",
                 fontsize=FONT["suptitle"])
    # margins given in inches, converted to fractions of the panel
    fig.subplots_adjust(left=0.32 / PANEL_W, right=1 - 0.05 / PANEL_W,
                        bottom=0.32 / PANEL_H, top=1 - 0.33 / PANEL_H, wspace=0.15, hspace=0.45)

    suffix = ".GC_depth" if kind == "gc" else ".position_depth"
    save(fig, SPECIES_LABELS.get(species, species).replace(" ", "_") + suffix)


# CV of window depth across MAGs
def cv_panel():
    cv = summary.assign(kit_rank=summary["Kit_Name"].map(KIT_ORDER.index))
    cv = cv.sort_values(["species", "kit_rank", "Replicate"], ascending=[False, False, False])
    # abbreviated genus
    species = cv["species"].replace(SPECIES_LABELS).str.replace(r"^(\w)\w* ", r"\1. ", regex=True)

    fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
    ax.barh(species + " | " + cv["Kit_Label"], cv["CV_Depth"], color=cv["Kit_Name"].map(KIT_COLORS))
    ax.set_xlabel("CV of 25-kb window depth", fontsize=FONT["label"], labelpad=1)
    ax.tick_params(axis="y", labelsize=FONT["label"], length=1.5, pad=1)
    ax.tick_params(axis="x", labelsize=FONT["tick"], length=1.5, pad=1)
    ax.grid(axis="x", alpha=0.2)
    # title on the whole panel
    fig.suptitle("Depth variability across MAGs", fontsize=FONT["suptitle"], fontweight="bold")
    fig.subplots_adjust(left=1.2 / PANEL_W, right=1 - 0.08 / PANEL_W,
                        bottom=0.3 / PANEL_H, top=1 - 0.33 / PANEL_H)
    save(fig, "all_MAGs.depth_CV")

for species in sorted(windows["species"].unique()):
    print(f">> {species}")
    species_panel(species, "gc")
    species_panel(species, "position")
cv_panel()

print("Done.")
