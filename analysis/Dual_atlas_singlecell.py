#!/usr/bin/env python3
"""
Fig 7G: Single-cell MyoScore validation across TWO independent atlases.

Atlas 1: HLMA (Lai et al. Nature 2024) — 292K cells, 15-99 yr
Atlas 2: Sanger SKM (Eraslan et al.) — 183K cells, 15-75 yr

Multi-panel:
  a) MyoScore by cell type, young vs old (HLMA)
  b) MyoScore by cell type, young vs old (Sanger) — replication
  c) Dimension heatmap Δ(young-old) — combined pattern
  d) MyoScore vs age (donor-level, both atlases)
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
mpl.rcParams['font.size'] = 8

COLOR_YOUNG = "#72c95e"
COLOR_OLD = "#50327b"
COLOR_HEALTHY = "#f4e030"
COLOR_UNHEALTHY = "#50327b"

from pathlib import Path

PROJECT_ROOT = Path(".")
OUT = str(PROJECT_ROOT / "MyoScore/figures")
DATA = str(PROJECT_ROOT / "MyoScore/data")

genes_df = pd.read_csv(
    PROJECT_ROOT / "GWAS_TWAS/gmhs_dimensions/dimension_genes_v33_telomere.csv"
)
DIMENSION_WEIGHTS = {'Strength': 0.252, 'Mass': 0.177, 'LeanMuscle': 0.243,
                     'Youth': 0.242, 'Resilience': 0.087}


def compute_myoscore(adata, genes_df, gene_col='var_names'):
    """Compute 5-dimension MyoScore on an AnnData object."""
    if gene_col == 'SYMBOL':
        # Sanger uses ENSEMBL as var_names, SYMBOL in var
        symbol_to_idx = {s: i for i, s in enumerate(adata.var['SYMBOL'].values)}
        sc_genes = set(adata.var['SYMBOL'].values)
    else:
        sc_genes = set(adata.var_names)

    dim_scores = {}
    for dim in ['Strength', 'Mass', 'LeanMuscle', 'Youth', 'Resilience']:
        dg = genes_df[genes_df['dimension'] == dim].copy()
        dg = dg[dg['ID'].isin(sc_genes)]
        dg = dg.sort_values('weight', ascending=False).drop_duplicates('ID')
        if len(dg) == 0:
            continue
        gene_list = dg['ID'].tolist()
        dirs = dg.set_index('ID')['direction_v3'].to_dict()
        wts = dg.set_index('ID')['weight'].to_dict()

        if gene_col == 'SYMBOL':
            idx = [symbol_to_idx[g] for g in gene_list]
            expr = adata.X[:, idx]
        else:
            expr = adata[:, gene_list].X

        if hasattr(expr, 'toarray'):
            expr = expr.toarray()

        mu = np.nanmean(expr, axis=0)
        sd = np.nanstd(expr, axis=0)
        sd[sd == 0] = 1
        expr_z = (expr - mu) / sd
        dir_arr = np.array([dirs[g] for g in gene_list])
        wt_arr = np.array([wts[g] for g in gene_list])
        wt_arr = wt_arr / wt_arr.sum()
        score = (expr_z * dir_arr) @ wt_arr
        dim_scores[dim] = score

    composite = np.zeros(adata.shape[0])
    for dim, w in DIMENSION_WEIGHTS.items():
        if dim in dim_scores:
            composite += w * dim_scores[dim]

    # Store in obs
    for dim in dim_scores:
        adata.obs[f'MyoScore_{dim}'] = dim_scores[dim]
    adata.obs['MyoScore'] = composite

    # Scale to 0-100
    for col in [c for c in adata.obs.columns if c.startswith('MyoScore')]:
        vals = adata.obs[col].values
        adata.obs[col] = 100 * (vals - vals.min()) / (vals.max() - vals.min())

    return adata


# ============================================================
# Load / compute Atlas 1: HLMA
# ============================================================
hlma_cache = f"{DATA}/sc_myoscore_obs_cache.parquet"
if os.path.exists(hlma_cache):
    print("Loading HLMA cached obs...")
    obs_hlma = pd.read_parquet(hlma_cache)
else:
    print("Computing HLMA MyoScore...")
    # HLMA atlas h5ad file: Lai et al. Nature 2024 (obtain from original publication)
    adata_hlma = sc.read_h5ad(f"{DATA}/HLMA_All_CellType_scsn_RNA.h5ad")
    if adata_hlma.X.max() > 50:
        sc.pp.normalize_total(adata_hlma, target_sum=1e4)
        sc.pp.log1p(adata_hlma)
    adata_hlma = compute_myoscore(adata_hlma, genes_df)
    obs_hlma = adata_hlma.obs.copy()
    obs_hlma.to_parquet(hlma_cache)

print(f"  HLMA: {len(obs_hlma)} cells, age range {obs_hlma['age'].min()}-{obs_hlma['age'].max()}")

# ============================================================
# Load / compute Atlas 2: Sanger SKM
# ============================================================
sanger_cache = f"{DATA}/sc_myoscore_sanger_obs_cache.parquet"
if os.path.exists(sanger_cache):
    print("Loading Sanger cached obs...")
    obs_sanger = pd.read_parquet(sanger_cache)
else:
    print("Computing Sanger MyoScore...")
    # Sanger SKM atlas h5ad file: Eraslan et al. (obtain from original publication)
    adata_sanger = sc.read_h5ad(str(PROJECT_ROOT / "single_cell_data/SKM_human_pp_cells2nuclei.h5ad"))
    if adata_sanger.X.max() > 50:
        sc.pp.normalize_total(adata_sanger, target_sum=1e4)
        sc.pp.log1p(adata_sanger)
    adata_sanger = compute_myoscore(adata_sanger, genes_df, gene_col='SYMBOL')
    obs_sanger = adata_sanger.obs.copy()
    obs_sanger.to_parquet(sanger_cache)

print(f"  Sanger: {len(obs_sanger)} cells")

# Map Sanger cell types to comparable categories
sanger_ct_map = {
    'MF-I': 'Type I', 'MF-II': 'Type II', 'Specialised MF': 'Specialized MF',
    'MuSC': 'MuSC', 'FB': 'FAP', 'CapEC': 'Endo', 'VenEC': 'Endo', 'ArtEC': 'Endo',
    'Pericyte': 'Pericyte', 'SMC': 'SMC', 'Macrophage': 'Myeloid cell',
    'Monocyte': 'Myeloid cell', 'T-cell': 'Lymphocyte', 'CD4+T': 'Lymphocyte',
    'CD8+T': 'Lymphocyte', 'NK-cell': 'Lymphocyte',
    'MF-IIsc(fg)': 'Type II', 'MF-Isc(fg)': 'Type I',
}
obs_sanger['CellType_mapped'] = obs_sanger['annotation_level0'].map(sanger_ct_map)

# ============================================================
# Statistics
# ============================================================
print("\n=== HLMA Young vs Old ===")
y_hlma = obs_hlma[obs_hlma['age_pop'] == 'young_pop']['MyoScore']
o_hlma = obs_hlma[obs_hlma['age_pop'] == 'old_pop']['MyoScore']
stat_h, p_h = stats.mannwhitneyu(y_hlma, o_hlma, alternative='greater')
d_h = (y_hlma.mean() - o_hlma.mean()) / np.sqrt((y_hlma.std()**2 + o_hlma.std()**2) / 2)
print(f"  Young={y_hlma.mean():.1f}, Old={o_hlma.mean():.1f}, P={p_h:.2e}, d={d_h:.3f}")

print("\n=== Sanger Young vs Old ===")
y_san = obs_sanger[obs_sanger['Age_bin'] == 'young']['MyoScore']
o_san = obs_sanger[obs_sanger['Age_bin'] == 'old']['MyoScore']
stat_s, p_s = stats.mannwhitneyu(y_san, o_san, alternative='greater')
d_s = (y_san.mean() - o_san.mean()) / np.sqrt((y_san.std()**2 + o_san.std()**2) / 2)
print(f"  Young={y_san.mean():.1f}, Old={o_san.mean():.1f}, P={p_s:.2e}, d={d_s:.3f}")

# ============================================================
# FIGURE
# ============================================================
print("\nCreating Fig 7G...")

fig = plt.figure(figsize=(14, 9))
gs = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.3)

common_cts = ['Type I', 'Type II', 'Specialized MF', 'MuSC', 'FAP',
              'Endo', 'Pericyte', 'SMC', 'Myeloid cell']


def plot_violin_panel(ax, obs_df, age_col, young_val, old_val, ct_col, title, panel_label):
    """Plot violin comparing young vs old by cell type."""
    labels = []
    for i, ct in enumerate(common_cts):
        y_data = obs_df[(obs_df[ct_col] == ct) & (obs_df[age_col] == young_val)]['MyoScore'].values
        o_data = obs_df[(obs_df[ct_col] == ct) & (obs_df[age_col] == old_val)]['MyoScore'].values
        if len(y_data) < 10 or len(o_data) < 10:
            labels.append(ct)
            continue
        # Subsample
        if len(y_data) > 3000:
            y_data = np.random.choice(y_data, 3000, replace=False)
        if len(o_data) > 3000:
            o_data = np.random.choice(o_data, 3000, replace=False)

        pos_y = i * 2.2
        pos_o = i * 2.2 + 0.9

        vp1 = ax.violinplot([y_data], positions=[pos_y], showmeans=True,
                             showextrema=False, widths=0.8)
        for body in vp1['bodies']:
            body.set_facecolor(COLOR_YOUNG)
            body.set_alpha(0.7)
            body.set_edgecolor('none')
        vp1['cmeans'].set_color('black')
        vp1['cmeans'].set_linewidth(1)

        vp2 = ax.violinplot([o_data], positions=[pos_o], showmeans=True,
                             showextrema=False, widths=0.8)
        for body in vp2['bodies']:
            body.set_facecolor(COLOR_OLD)
            body.set_alpha(0.7)
            body.set_edgecolor('none')
        vp2['cmeans'].set_color('black')
        vp2['cmeans'].set_linewidth(1)

        # Significance
        _, p = stats.mannwhitneyu(y_data, o_data, alternative='greater')
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
        y_top = max(np.percentile(y_data, 95), np.percentile(o_data, 95))
        color_sig = 'black' if sig != 'ns' else 'gray'
        ax.text(i * 2.2 + 0.45, y_top + 1.5, sig, ha='center', fontsize=6, color=color_sig)

        labels.append(ct)

    ax.set_xticks([i * 2.2 + 0.45 for i in range(len(labels))])
    ax.set_xticklabels(labels, rotation=40, ha='right', fontsize=6.5)
    ax.set_ylabel('MyoScore', fontsize=8)
    ax.set_ylim(10, 55)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(title, fontsize=8)
    ax.text(-0.05, 1.05, panel_label, transform=ax.transAxes, fontsize=11, fontweight='bold')


# --- Panel a: HLMA violin ---
ax_a = fig.add_subplot(gs[0, 0])
plot_violin_panel(ax_a, obs_hlma, 'age_pop', 'young_pop', 'old_pop', 'Annotation',
                  'HLMA (Lai et al., n=292,423)', 'a')

legend_elements = [Patch(facecolor=COLOR_YOUNG, alpha=0.7, label='Young'),
                   Patch(facecolor=COLOR_OLD, alpha=0.7, label='Old')]
ax_a.legend(handles=legend_elements, loc='upper right', fontsize=6, frameon=False)

# --- Panel b: Sanger violin ---
ax_b = fig.add_subplot(gs[0, 1])
plot_violin_panel(ax_b, obs_sanger, 'Age_bin', 'young', 'old', 'CellType_mapped',
                  'Sanger Atlas (n=183,161)', 'b')

# --- Panel c: Dimension heatmap (HLMA + Sanger combined) ---
ax_c = fig.add_subplot(gs[1, 0])

dims = ['Strength', 'Mass', 'LeanMuscle', 'Youth', 'Resilience']
muscle_cts = ['Type I', 'Type II', 'Specialized MF', 'MuSC', 'FAP', 'Endo', 'Myeloid cell']

# Average Δ from both atlases
delta_combined = np.zeros((len(muscle_cts), len(dims)))

for i, ct in enumerate(muscle_cts):
    for j, dim in enumerate(dims):
        col = f'MyoScore_{dim}'
        deltas = []

        # HLMA
        if col in obs_hlma.columns:
            y = obs_hlma[(obs_hlma['Annotation'] == ct) & (obs_hlma['age_pop'] == 'young_pop')][col]
            o = obs_hlma[(obs_hlma['Annotation'] == ct) & (obs_hlma['age_pop'] == 'old_pop')][col]
            if len(y) > 10 and len(o) > 10:
                deltas.append(y.mean() - o.mean())

        # Sanger
        if col in obs_sanger.columns:
            y = obs_sanger[(obs_sanger['CellType_mapped'] == ct) & (obs_sanger['Age_bin'] == 'young')][col]
            o = obs_sanger[(obs_sanger['CellType_mapped'] == ct) & (obs_sanger['Age_bin'] == 'old')][col]
            if len(y) > 10 and len(o) > 10:
                deltas.append(y.mean() - o.mean())

        if deltas:
            delta_combined[i, j] = np.mean(deltas)

cmap_delta = LinearSegmentedColormap.from_list('delta', [COLOR_UNHEALTHY, '#e8e8e8', COLOR_HEALTHY])
im = ax_c.imshow(delta_combined, cmap=cmap_delta, aspect='auto', vmin=-4, vmax=4)

for i in range(len(muscle_cts)):
    for j in range(len(dims)):
        val = delta_combined[i, j]
        color = 'white' if abs(val) > 2.5 else 'black'
        ax_c.text(j, i, f'{val:+.1f}', ha='center', va='center', fontsize=6.5, color=color)

ax_c.set_xticks(range(len(dims)))
ax_c.set_xticklabels(dims, rotation=30, ha='right', fontsize=7)
ax_c.set_yticks(range(len(muscle_cts)))
ax_c.set_yticklabels(muscle_cts, fontsize=7)
ax_c.set_title('Δ MyoScore (Young − Old)\nCombined across two atlases', fontsize=8)
ax_c.text(-0.08, 1.05, 'c', transform=ax_c.transAxes, fontsize=11, fontweight='bold')

divider = make_axes_locatable(ax_c)
cax = divider.append_axes("right", size="5%", pad=0.08)
cb = plt.colorbar(im, cax=cax)
cb.set_label('Δ Score', fontsize=7)
cb.ax.tick_params(labelsize=6)

# --- Panel d: MyoScore vs Age (myofiber-only, donor-level, colored by atlas) ---
ax_d = fig.add_subplot(gs[1, 1])

COLOR_HLMA = "#e07b39"    # warm orange for HLMA
COLOR_SANGER = "#2878b5"  # cool blue for Sanger

# HLMA myofiber-only donors
hlma_mf = obs_hlma[obs_hlma['Annotation'].isin(['Type I', 'Type II'])]
donor_hlma = hlma_mf.groupby('sample').agg({
    'MyoScore': 'mean', 'MyoScore_Youth': 'mean',
    'age': 'first', 'age_pop': 'first'
}).reset_index()
donor_hlma['atlas'] = 'HLMA'

# Sanger myofiber-only donors
sanger_ct_map_local = {'MF-I': 'Type I', 'MF-II': 'Type II',
                       'MF-IIsc(fg)': 'Type II', 'MF-Isc(fg)': 'Type I'}
obs_sanger['CellType_mf'] = obs_sanger['annotation_level0'].map(sanger_ct_map_local)
sanger_mf = obs_sanger[obs_sanger['CellType_mf'].isin(['Type I', 'Type II'])].copy()
age_map = {'15-20': 17.5, '20-25': 22.5, '25-30': 27.5, '35-40': 37.5,
           '50-55': 52.5, '55-60': 57.5, '60-65': 62.5, '70-75': 72.5}
sanger_mf['age_numeric'] = sanger_mf['Age_group'].map(age_map)
donor_sanger = sanger_mf.groupby('DonorID').agg({
    'MyoScore': 'mean', 'MyoScore_Youth': 'mean',
    'age_numeric': 'first', 'Age_bin': 'first'
}).reset_index()
donor_sanger.rename(columns={'age_numeric': 'age', 'Age_bin': 'age_pop'}, inplace=True)
donor_sanger['atlas'] = 'Sanger'
donor_sanger['age_pop'] = donor_sanger['age_pop'].map({'young': 'young_pop', 'old': 'old_pop'})

donors_all = pd.concat([donor_hlma[['MyoScore', 'MyoScore_Youth', 'age', 'age_pop', 'atlas']],
                        donor_sanger[['MyoScore', 'MyoScore_Youth', 'age', 'age_pop', 'atlas']]])

# Z-score normalize MyoScore within each atlas
for col in ['MyoScore', 'MyoScore_Youth']:
    for atlas in ['HLMA', 'Sanger']:
        mask = donors_all['atlas'] == atlas
        vals = donors_all.loc[mask, col]
        donors_all.loc[mask, col + '_z'] = (vals - vals.mean()) / vals.std()

# Plot colored by atlas source
for _, row in donors_all.iterrows():
    color = COLOR_HLMA if row['atlas'] == 'HLMA' else COLOR_SANGER
    marker = 'o' if row['atlas'] == 'HLMA' else 's'
    ax_d.scatter(row['age'], row['MyoScore_z'], c=color, s=55, alpha=0.85,
                edgecolors='white', linewidth=0.5, marker=marker, zorder=3)

# Per-atlas regression lines
for atlas, color, ls in [('HLMA', COLOR_HLMA, '-'), ('Sanger', COLOR_SANGER, '--')]:
    sub = donors_all[donors_all['atlas'] == atlas]
    r_a, p_a = stats.spearmanr(sub['age'], sub['MyoScore_z'])
    slope_a, intercept_a = np.polyfit(sub['age'], sub['MyoScore_z'], 1)
    x_a = np.linspace(sub['age'].min() - 2, sub['age'].max() + 2, 50)
    ax_d.plot(x_a, slope_a * x_a + intercept_a, color=color, linewidth=1.2,
             linestyle=ls, alpha=0.6)

# Overall stats
r, p = stats.spearmanr(donors_all['age'], donors_all['MyoScore_z'])

# Per-atlas stats for annotation
r_h, p_h_corr = stats.spearmanr(donor_hlma['age'], donor_hlma['MyoScore'])
r_s, p_s_corr = stats.spearmanr(donor_sanger['age'], donor_sanger['MyoScore'])

ax_d.set_xlabel('Age (years)', fontsize=8)
ax_d.set_ylabel('MyoScore (myofiber)\n(within-atlas Z-score)', fontsize=8)
ax_d.axhline(y=0, color='gray', linewidth=0.5, linestyle=':')

# Stats text
stats_text = (f'Combined: ρ={r:.2f}, P={p:.3f}\n'
              f'HLMA: ρ={r_h:.2f}, P={p_h_corr:.3f}\n'
              f'Sanger: ρ={r_s:.2f}, P={p_s_corr:.3f}\n'
              f'n={len(donors_all)} donors (myofibers)')
ax_d.text(0.03, 0.97, stats_text, transform=ax_d.transAxes, fontsize=6, va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='gray', linewidth=0.5))

ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)
ax_d.text(-0.1, 1.05, 'd', transform=ax_d.transAxes, fontsize=11, fontweight='bold')

# Legend
from matplotlib.lines import Line2D
leg_d = [Line2D([0],[0], marker='o', color='w', markerfacecolor=COLOR_HLMA, markersize=6, linestyle='', label='HLMA (15–99 yr)'),
         Line2D([0],[0], marker='s', color='w', markerfacecolor=COLOR_SANGER, markersize=6, linestyle='', label='Sanger (15–75 yr)')]
ax_d.legend(handles=leg_d, loc='upper right', fontsize=6, frameon=False)

# Save
fig.savefig(f"{OUT}/Fig7G_singlecell_validation.pdf", dpi=300, bbox_inches='tight')
fig.savefig(f"{OUT}/Fig7G_singlecell_validation.png", dpi=300, bbox_inches='tight')
print("  Saved Fig 7G")

# Save donor scores
donors_all.to_csv(f"{DATA}/sc_myoscore_donor_scores_both_atlases.csv", index=False)

# ============================================================
# Summary
# ============================================================
print("\n=== MANUSCRIPT TEXT ===")
print(f"Fig 7G: Single-cell MyoScore validation across two independent atlases")
print(f"  Total cells: {len(obs_hlma) + len(obs_sanger):,}")
print(f"  Total donors: {len(donors_all)}")
print(f"  HLMA: Young vs Old P = {p_h:.2e}, d = {d_h:.3f}")
print(f"  Sanger: Young vs Old P = {p_s:.2e}, d = {d_s:.3f}")
print(f"  MyoScore vs age (combined): ρ = {r:.2f}, P = {p:.1e}")
print(f"\nDone!")
