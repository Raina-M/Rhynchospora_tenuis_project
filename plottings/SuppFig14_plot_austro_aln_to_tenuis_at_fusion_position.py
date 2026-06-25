#!/usr/bin/env python3
"""
Two stacked genome-mapping plots (TSV1 top, TSV2 bottom) on one PDF page.
- Input TSV format (tab-delim, no header):
    col1: sp1_chr
    col2: sp1_start
    col3: sp1_end
    col4: sp2_chr   (expected 'chr1' for all rows)
    col5: sp2_start
    col6: sp2_end
- X axis displayed in Mbp (position / 1e6)
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties, findSystemFonts

# -----------------------------
# User parameters - change these
# -----------------------------
tsv1 = "austro_to_REC_h1_synteny.tsv"
tsv2 = "austro_to_REC_h2_synteny.tsv"
out_pdf = "two_mappings.pdf"
figure_size = (8, 9)   # width, height in inches
hspace = 2            # vertical spacing between the two plots (subplots_adjust hspace)
rect_height = 0.5        # rectangle height in track units
tick_fontsize = 13
label_fontsize = 13
title_fontsize = 13
# chr1 lengths
max_pos1 = 212800643
max_pos2 = 215061937

# Set font
font_path = "/home/mzhang/.local/share/fonts/arial/ARIAL.TTF"
arial_font = FontProperties(fname=font_path)

# -----------------------------
# Load files function
# -----------------------------
def load_mapping_tsv(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"TSV file not found: {path}")
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["sp1_chr", "sp1_start", "sp1_end",
                            "sp2_chr", "sp2_start", "sp2_end"],
                     dtype={"sp1_chr": str, "sp2_chr": str})
    # sanity: ensure numeric coordinates
    for c in ("sp1_start","sp1_end","sp2_start","sp2_end"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    # drop rows with NA coords
    df = df.dropna(subset=["sp2_start", "sp2_end"])
    return df

# -----------------------------
# Plotting function for one dataframe into a given axis
# -----------------------------
def plot_stacked_tracks(ax, df, max_pos, arial_fontprop=None, title_text=""):
    # Desired order bottom->top
    desired_order = ["chr2_scaffold_7", "chr2_scaffold_8", "chr2_scaffold_9",
                 "chr2_scaffold_10", "chr2_scaffold_11", "chr2_scaffold_12",
                 "chr3_scaffold_13", "chr3_scaffold_14", "chr3_scaffold_15",
                 "chr3_scaffold_16", "chr3_scaffold_17", "chr3_scaffold_18"]

    # Keep only chromosomes present in this df, respecting desired order
    present = list(df["sp1_chr"].unique())
    sp1_chrs = [c for c in desired_order if c in present]

    # map chr -> y index (1..n)
    chr_to_y = {chr_name: i+1 for i, chr_name in enumerate(sp1_chrs)}
    n_chrs = len(sp1_chrs)

    # color map (consistent within this plot)
    color_map = {chr_id: plt.cm.tab20(i % 20) for i, chr_id in enumerate(sp1_chrs)}

    # plot rectangles (convert positions to Mbp)
    for _, row in df.iterrows():
        chr_name = row["sp1_chr"]
        if chr_name not in chr_to_y:
            continue
        y = chr_to_y[chr_name]
        start_m = row["sp2_start"] / 1e6
        end_m = row["sp2_end"] / 1e6
        width = max(end_m - start_m, 1e-6)  # avoid zero width
        rect = patches.Rectangle((start_m, y - rect_height/2),
                                 width=width,
                                 height=rect_height,
                                 facecolor=color_map[chr_name],
                                 linewidth=0,
                                 alpha=0.9)
        ax.add_patch(rect)

    # Formatting
    #ax.set_xlim(0, max_pos / 1e6)
    ax.set_xlim(127, 131.4)
    for label in ax.get_xticklabels():
        label.set_fontproperties(arial_font)
        label.set_fontsize(tick_fontsize)
    
    ax.set_ylim(0.5, n_chrs + 0.5)
    ax.set_yticks(list(range(1, n_chrs + 1)))
    ax.set_yticklabels(["Chr2 scaffold 7", "Chr2 scaffold 8", "Chr2 scaffold 9",
                        "Chr2 scaffold 10", "Chr2 scaffold 11", "Chr2 scaffold 12",
                        "Chr3 scaffold 13", "Chr3 scaffold 14", "Chr3 scaffold 15",
                        "Chr3 scaffold 16", "Chr3 scaffold 17", "Chr3 scaffold 18"],
                        fontsize=tick_fontsize, fontproperties=arial_font)
    ax.set_xlabel("R. tenuis Chr1 position (Mbp)", fontproperties=arial_fontprop, fontsize=label_fontsize)
    ax.set_ylabel("R. austrobrasiliensis chromosomes", fontproperties=arial_fontprop, fontsize=label_fontsize)
    if title_text:
        ax.set_title(title_text, fontproperties=arial_font, fontsize=title_fontsize)

    # light vertical grid lines
    ax.grid(axis='x', linestyle='--', alpha=0.25)
    

# -----------------------------
# Main run
# -----------------------------
def main():

    # Load dataframes
    df1 = load_mapping_tsv(tsv1)
    df2 = load_mapping_tsv(tsv2)

    # create PDF and figure
    with PdfPages(out_pdf) as pdf:
        fig = plt.figure(figsize=figure_size)

        # top plot (TSV1)
        ax1 = fig.add_subplot(2, 1, 1)
        plot_stacked_tracks(ax1, df1, max_pos1, arial_fontprop=arial_font, title_text="Mapping to REC Chr1 h1")

        # bottom plot (TSV2)
        ax2 = fig.add_subplot(2, 1, 2)
        plot_stacked_tracks(ax2, df2, max_pos2, arial_fontprop=arial_font, title_text="Mapping to REC Chr1 h2")

        # Adjust spacing between top & bottom plots
        fig.subplots_adjust(left=0.15, right=0.95, hspace=hspace)

        # Save page
        pdf.savefig(fig)
        plt.close(fig)

    print(f"Saved PDF to: {out_pdf}")

if __name__ == "__main__":
    main()

