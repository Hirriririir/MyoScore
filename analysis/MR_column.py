#!/usr/bin/env python3
"""
Fig 7F: MR direction column — designed to align with Fig 7A right side.
Same gene order, same row height, narrow format for Illustrator splicing.

Outputs:
  Fig7F_MR_muscle_eqtl.pdf  — standalone heatmap (original)
  Fig7F_MR_column.pdf       — narrow column matching 7A row layout
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
mpl.rcParams['font.size'] = 8

from pathlib import Path
PROJECT_ROOT = Path(".")
OUT = str(PROJECT_ROOT / "MyoScore/figures")

# Load MR results
mr = pd.read_csv(PROJECT_ROOT / "GWAS_TWAS/mr_results/MR1_muscle_eqtl_results.csv")

# Gene order matching Fig 7A exactly (top to bottom)
gene_order = ['ACSS2', 'TMEM52', 'GGT7',  # positive (↑H)
              'CEP250', 'UQCC1', 'CPNE1', 'RSRC2', 'RBL2', 'SNRPC', 'GATAD1']  # negative (↑D)

gene_directions = {
    'ACSS2': +1, 'GGT7': +1, 'TMEM52': +1,
    'CEP250': -1, 'UQCC1': -1, 'CPNE1': -1,
    'RSRC2': -1, 'RBL2': -1, 'SNRPC': -1, 'GATAD1': -1
}

outcome_order = ['Grip strength', 'Walking pace', 'Appendicular lean mass', 'Whole body FFM']
outcome_short = ['Grip', 'Pace', 'ALM', 'FFM']

# Build concordance matrix
n_genes = len(gene_order)
n_out = len(outcome_order)
conc_matrix = np.zeros((n_genes, n_out))
pval_matrix = np.ones((n_genes, n_out))
beta_matrix = np.zeros((n_genes, n_out))

for i, gene in enumerate(gene_order):
    for j, outcome in enumerate(outcome_order):
        row = mr[(mr['gene'] == gene) & (mr['outcome'] == outcome)]
        if len(row) > 0:
            row = row.iloc[0]
            beta = row['wald_beta']
            p = row['wald_p']
            d = gene_directions[gene]
            beta_matrix[i, j] = beta
            pval_matrix[i, j] = p
            # Concordant if: positive gene has positive beta, negative gene has negative beta
            if d == 1:
                conc_matrix[i, j] = 1 if beta > 0 else -1
            else:
                conc_matrix[i, j] = 1 if beta < 0 else -1

# ============================================================
# Narrow column version (for splicing with Fig 7A)
# ============================================================
fig, ax = plt.subplots(figsize=(3.2, 5.5))  # narrow to match 7A height

# Use diamond markers like Fig 7A
for i, gene in enumerate(gene_order):
    for j, outcome in enumerate(outcome_order):
        conc = conc_matrix[i, j]
        p = pval_matrix[i, j]
        beta = beta_matrix[i, j]

        # Color: blue=concordant, red=discordant, gray=ns
        if p >= 0.05:
            color = '#cccccc'
            edge = '#999999'
        elif conc > 0:
            color = '#4575b4'  # concordant blue
            edge = '#2b5c8a'
        else:
            color = '#d73027'  # discordant red
            edge = '#a52019'

        # Size proportional to -log10(p)
        size = min(max(-np.log10(p) * 15, 20), 120)

        ax.scatter(j, n_genes - 1 - i, marker='D', s=size, c=color,
                  edgecolors=edge, linewidth=0.5, zorder=3)

        # Significance stars
        if p < 0.001:
            sig = '***'
        elif p < 0.01:
            sig = '**'
        elif p < 0.05:
            sig = '*'
        else:
            sig = ''

        if sig:
            ax.text(j, n_genes - 1 - i - 0.35, sig, ha='center', va='top',
                   fontsize=5, color=color)

# Grid lines
for i in range(n_genes):
    ax.axhline(y=i, color='#eeeeee', linewidth=0.3, zorder=0)
for j in range(n_out):
    ax.axvline(x=j, color='#eeeeee', linewidth=0.3, zorder=0)

# Separator between positive and negative genes
ax.axhline(y=n_genes - 3.5, color='black', linewidth=0.8, linestyle='--')

# Labels
ax.set_yticks(range(n_genes))
ax.set_yticklabels(gene_order[::-1], fontsize=7)
ax.set_xticks(range(n_out))
ax.set_xticklabels(outcome_short, rotation=45, ha='right', fontsize=7)

ax.set_xlim(-0.5, n_out - 0.5)
ax.set_ylim(-0.5, n_genes - 0.5)

# Title
ax.set_title('Muscle eQTL\ncis-MR', fontsize=8, fontweight='bold')

# Remove box
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.tick_params(left=False, bottom=False)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#4575b4',
           markeredgecolor='#2b5c8a', markersize=7, label='Concordant (P<0.05)'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#d73027',
           markeredgecolor='#a52019', markersize=7, label='Discordant (P<0.05)'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#cccccc',
           markeredgecolor='#999999', markersize=7, label='Not significant'),
]
ax.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.22),
         fontsize=5.5, frameon=False, ncol=1)

plt.tight_layout()
fig.savefig(f"{OUT}/Fig7F_MR_column.pdf", dpi=300, bbox_inches='tight')
fig.savefig(f"{OUT}/Fig7F_MR_column.png", dpi=300, bbox_inches='tight')
print("Saved Fig7F_MR_column.pdf/.png")

# Print concordance summary
concordant = int((conc_matrix > 0).sum())
sig_concordant = int(((conc_matrix > 0) & (pval_matrix < 0.05)).sum())
total_sig = int((pval_matrix < 0.05).sum())
print(f"\nDirection concordant: {concordant}/{n_genes * n_out} ({100*concordant/(n_genes*n_out):.0f}%)")
print(f"Significant & concordant: {sig_concordant}/{total_sig} significant tests")

# Per-gene summary
print(f"\nPer-gene concordance:")
for i, gene in enumerate(gene_order):
    n_conc = int((conc_matrix[i] > 0).sum())
    n_sig = int((pval_matrix[i] < 0.05).sum())
    print(f"  {gene:8s}: {n_conc}/4 concordant, {n_sig}/4 significant")
