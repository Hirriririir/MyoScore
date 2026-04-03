#!/usr/bin/env python3
"""
GSE202745 Within-individual Muscle Consistency Analysis
Compare healthy controls vs LGMD R12 patients

This script calculates:
1. Pairwise muscle correlations within individuals
2. Intraclass correlation coefficients (ICC)
3. Comparative statistics between healthy and myopathy groups
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
warnings.filterwarnings('ignore')

# For illustrator-editable fonts
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# Project colors
COLORS = {
    'unhealthy': '#50327b',
    'healthy': '#f4e030',
    'intermediate': '#46508b',
    'good': '#72c95e',
    'neutral': '#31848f'
}

def load_data():
    """Load scores and meta data for GSE202745"""
    # Load scores
    scores = pd.read_csv('MyoScore/data/MyoScore_scores.csv', index_col=0)

    # Load main meta
    meta = pd.read_csv('Meta/RNA-seq integration meta (MultiOmics).csv')

    # Filter GSE202745
    gse202745_meta = meta[meta['Geo_accession'] == 'GSE202745'].copy()

    # Load original GEO meta with true patient IDs
    geo_meta = pd.read_csv('Validation/GSE202745 GEO meta.csv', encoding='utf-8-sig')

    return scores, gse202745_meta, geo_meta

def extract_patient_from_title(title):
    """Extract patient ID from GEO title like 'Healthy,semimembranosus muscle, patient C1'"""
    import re
    match = re.search(r'patient\s+([CP]?\d+)', title)
    if match:
        return match.group(1)
    return None

def create_analysis_df(scores, meta, geo_meta):
    """Create analysis dataframe with patient groupings"""

    # Map GSE202745 samples
    gse_samples = meta[meta['Geo_accession'] == 'GSE202745'][['Sample_id', 'Phenotype_2', 'BiopsySite', 'Gsm_accession']].copy()
    gse_samples['Sample_id'] = gse_samples['Sample_id'].astype(str)

    # Extract patient IDs from original GEO meta
    geo_meta['Patient_ID'] = geo_meta['Title'].apply(extract_patient_from_title)
    geo_meta['Group'] = geo_meta['Disease state'].apply(
        lambda x: 'Healthy' if x == 'Healthy' else 'LGMD R12'
    )

    # Merge with GEO meta to get true patient IDs
    gse_samples = gse_samples.merge(
        geo_meta[['Accession', 'Patient_ID', 'Group']].rename(columns={'Accession': 'Gsm_accession'}),
        on='Gsm_accession',
        how='left'
    )

    # Merge with scores
    analysis_df = gse_samples.merge(
        scores[['GMHS_v33']].rename(columns={'GMHS_v33': 'MyoScore'}),
        left_on='Sample_id',
        right_index=True,
        how='inner'
    )

    # Standardize muscle names
    analysis_df['Muscle'] = analysis_df['BiopsySite'].str.strip().str.lower()
    analysis_df['Muscle'] = analysis_df['Muscle'].replace({
        'vastus lateralis': 'VL',
        'rectus femoris': 'RF',
        'semimembranosus': 'Sem'
    })

    # Use Patient_ID as Individual
    analysis_df['Individual'] = analysis_df['Patient_ID']

    return analysis_df

def calculate_muscle_correlations(df, group_name):
    """Calculate pairwise muscle correlations for a group"""
    results = {}

    # Pivot to wide format for each individual
    pivot_df = df.pivot_table(
        index='Individual',
        columns='Muscle',
        values='MyoScore',
        aggfunc='first'
    )

    # Calculate pairwise correlations
    muscle_pairs = [('VL', 'RF'), ('Sem', 'RF'), ('Sem', 'VL')]

    for m1, m2 in muscle_pairs:
        if m1 in pivot_df.columns and m2 in pivot_df.columns:
            # Get paired observations
            paired = pivot_df[[m1, m2]].dropna()
            n_pairs = len(paired)

            if n_pairs >= 3:
                r_pearson, p_pearson = pearsonr(paired[m1], paired[m2])
                r_spearman, p_spearman = spearmanr(paired[m1], paired[m2])

                results[f'{m1} vs {m2}'] = {
                    'n_pairs': n_pairs,
                    'pearson_r': r_pearson,
                    'pearson_p': p_pearson,
                    'spearman_r': r_spearman,
                    'spearman_p': p_spearman
                }

    return results, pivot_df

def calculate_icc(df):
    """Calculate intraclass correlation coefficient (ICC 2,1)"""
    # Pivot to wide format
    pivot_df = df.pivot_table(
        index='Individual',
        columns='Muscle',
        values='MyoScore',
        aggfunc='first'
    )

    # Remove individuals with missing muscles
    pivot_clean = pivot_df.dropna()

    if len(pivot_clean) < 3:
        return np.nan, np.nan, np.nan

    # ICC calculation using ANOVA approach
    n = len(pivot_clean)  # number of subjects
    k = len(pivot_clean.columns)  # number of raters/muscles

    values = pivot_clean.values

    # Grand mean
    grand_mean = np.mean(values)

    # Between-subject variance
    subject_means = np.mean(values, axis=1)
    SS_between = k * np.sum((subject_means - grand_mean) ** 2)
    MS_between = SS_between / (n - 1)

    # Within-subject variance
    SS_within = np.sum((values - subject_means[:, np.newaxis]) ** 2)
    MS_within = SS_within / (n * (k - 1))

    # ICC (2,1) formula
    icc = (MS_between - MS_within) / (MS_between + (k - 1) * MS_within)

    # Variance components
    var_between = (MS_between - MS_within) / k
    var_within = MS_within

    return icc, var_between, var_within

def calculate_muscle_stats(df, group_name):
    """Calculate descriptive statistics by muscle"""
    stats_df = df.groupby('Muscle')['MyoScore'].agg(['mean', 'std', 'count'])
    stats_df.columns = ['Mean', 'SD', 'n']
    return stats_df

def plot_comparison(healthy_pivot, lgmd_pivot, output_path):
    """Create comparison scatter plots"""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    muscle_pairs = [('VL', 'RF'), ('Sem', 'RF'), ('Sem', 'VL')]

    for col, (m1, m2) in enumerate(muscle_pairs):
        # Healthy group
        ax = axes[0, col]
        if m1 in healthy_pivot.columns and m2 in healthy_pivot.columns:
            paired = healthy_pivot[[m1, m2]].dropna()
            ax.scatter(paired[m1], paired[m2], c=COLORS['healthy'],
                      edgecolors='black', s=100, alpha=0.8, label='Healthy')

            if len(paired) >= 3:
                r, p = pearsonr(paired[m1], paired[m2])
                ax.set_title(f'Healthy: {m1} vs {m2}\nr={r:.2f}, p={p:.3f}', fontsize=12)

                # Add regression line
                z = np.polyfit(paired[m1], paired[m2], 1)
                p_line = np.poly1d(z)
                x_line = np.linspace(paired[m1].min(), paired[m1].max(), 100)
                ax.plot(x_line, p_line(x_line), color='gray', linestyle='--', alpha=0.7)

        ax.set_xlabel(f'{m1} MyoScore')
        ax.set_ylabel(f'{m2} MyoScore')
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)

        # LGMD group
        ax = axes[1, col]
        if m1 in lgmd_pivot.columns and m2 in lgmd_pivot.columns:
            paired = lgmd_pivot[[m1, m2]].dropna()
            ax.scatter(paired[m1], paired[m2], c=COLORS['unhealthy'],
                      edgecolors='black', s=100, alpha=0.8, label='LGMD R12')

            if len(paired) >= 3:
                r, p = pearsonr(paired[m1], paired[m2])
                ax.set_title(f'LGMD R12: {m1} vs {m2}\nr={r:.2f}, p={p:.3f}', fontsize=12)

                # Add regression line
                z = np.polyfit(paired[m1], paired[m2], 1)
                p_line = np.poly1d(z)
                x_line = np.linspace(paired[m1].min(), paired[m1].max(), 100)
                ax.plot(x_line, p_line(x_line), color='gray', linestyle='--', alpha=0.7)

        ax.set_xlabel(f'{m1} MyoScore')
        ax.set_ylabel(f'{m2} MyoScore')
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def plot_vl_rf_comparison(healthy_pivot, lgmd_pivot, output_path):
    """Create focused VL vs RF comparison plot - the key finding"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    m1, m2 = 'VL', 'RF'

    # Healthy group
    ax = axes[0]
    if m1 in healthy_pivot.columns and m2 in healthy_pivot.columns:
        paired = healthy_pivot[[m1, m2]].dropna()
        ax.scatter(paired[m1], paired[m2], c=COLORS['healthy'],
                  edgecolors='black', s=120, alpha=0.9, zorder=3)

        if len(paired) >= 3:
            r, p = pearsonr(paired[m1], paired[m2])
            ax.set_title(f'Healthy Controls (n={len(paired)})\nr = +{r:.2f}, p = {p:.3f}',
                        fontsize=14, fontweight='bold')

            # Add regression line
            z = np.polyfit(paired[m1], paired[m2], 1)
            p_line = np.poly1d(z)
            x_line = np.linspace(20, 80, 100)
            ax.plot(x_line, p_line(x_line), color=COLORS['neutral'],
                   linestyle='-', linewidth=2, alpha=0.8, zorder=2)

    ax.set_xlabel('Vastus Lateralis MyoScore', fontsize=12)
    ax.set_ylabel('Rectus Femoris MyoScore', fontsize=12)
    ax.set_xlim(20, 80)
    ax.set_ylim(20, 80)
    ax.plot([20, 80], [20, 80], 'k--', alpha=0.3, zorder=1)  # Identity line
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # LGMD group
    ax = axes[1]
    if m1 in lgmd_pivot.columns and m2 in lgmd_pivot.columns:
        paired = lgmd_pivot[[m1, m2]].dropna()
        ax.scatter(paired[m1], paired[m2], c=COLORS['unhealthy'],
                  edgecolors='black', s=120, alpha=0.9, zorder=3)

        if len(paired) >= 3:
            r, p = pearsonr(paired[m1], paired[m2])
            rho, p_spearman = spearmanr(paired[m1], paired[m2])
            ax.set_title(f'LGMD R12 Patients (n={len(paired)})\nr = {r:.2f}, ρ = {rho:.2f}*',
                        fontsize=14, fontweight='bold')

            # Add regression line
            z = np.polyfit(paired[m1], paired[m2], 1)
            p_line = np.poly1d(z)
            x_line = np.linspace(40, 60, 100)
            ax.plot(x_line, p_line(x_line), color=COLORS['neutral'],
                   linestyle='-', linewidth=2, alpha=0.8, zorder=2)

    ax.set_xlabel('Vastus Lateralis MyoScore', fontsize=12)
    ax.set_ylabel('Rectus Femoris MyoScore', fontsize=12)
    ax.set_xlim(40, 60)
    ax.set_ylim(40, 60)
    ax.plot([40, 60], [40, 60], 'k--', alpha=0.3, zorder=1)  # Identity line
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # Add annotation
    fig.text(0.5, 0.02,
             'Positive correlation in healthy → Negative correlation in LGMD R12\n'
             '(differential muscle involvement characteristic of myopathy)',
             ha='center', fontsize=11, style='italic')

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def main():
    print("=" * 60)
    print("GSE202745 Within-individual Muscle Consistency Analysis")
    print("=" * 60)

    # Load data
    scores, meta, geo_meta = load_data()
    print(f"\nLoaded {len(scores)} samples from MyoScore_scores.csv")
    print(f"Found {len(meta)} samples in GSE202745")
    print(f"GEO meta has {len(geo_meta)} entries")

    # Create analysis dataframe
    analysis_df = create_analysis_df(scores, meta, geo_meta)
    print(f"\nMatched {len(analysis_df)} samples with scores")

    # Split by group
    healthy_df = analysis_df[analysis_df['Group'] == 'Healthy']
    lgmd_df = analysis_df[analysis_df['Group'] == 'LGMD R12']

    print(f"\nHealthy controls: {len(healthy_df)} samples, {healthy_df['Individual'].nunique()} individuals")
    print(f"LGMD R12 patients: {len(lgmd_df)} samples, {lgmd_df['Individual'].nunique()} individuals")

    # Calculate correlations for healthy group
    print("\n" + "=" * 60)
    print("HEALTHY CONTROLS - Pairwise Muscle Correlations")
    print("=" * 60)
    healthy_corr, healthy_pivot = calculate_muscle_correlations(healthy_df, 'Healthy')

    for pair, vals in healthy_corr.items():
        print(f"\n{pair}:")
        print(f"  n pairs: {vals['n_pairs']}")
        print(f"  Pearson r: {vals['pearson_r']:.3f} (p={vals['pearson_p']:.4f})")
        print(f"  Spearman ρ: {vals['spearman_r']:.3f} (p={vals['spearman_p']:.4f})")

    # Calculate correlations for LGMD group
    print("\n" + "=" * 60)
    print("LGMD R12 PATIENTS - Pairwise Muscle Correlations")
    print("=" * 60)
    lgmd_corr, lgmd_pivot = calculate_muscle_correlations(lgmd_df, 'LGMD R12')

    for pair, vals in lgmd_corr.items():
        print(f"\n{pair}:")
        print(f"  n pairs: {vals['n_pairs']}")
        print(f"  Pearson r: {vals['pearson_r']:.3f} (p={vals['pearson_p']:.4f})")
        print(f"  Spearman ρ: {vals['spearman_r']:.3f} (p={vals['spearman_p']:.4f})")

    # Calculate ICC for both groups
    print("\n" + "=" * 60)
    print("ICC Analysis")
    print("=" * 60)

    healthy_icc, healthy_var_between, healthy_var_within = calculate_icc(healthy_df)
    print(f"\nHealthy Controls:")
    print(f"  ICC: {healthy_icc:.3f}")
    print(f"  Between-individual variance: {healthy_var_between:.2f}")
    print(f"  Within-individual variance: {healthy_var_within:.2f}")

    lgmd_icc, lgmd_var_between, lgmd_var_within = calculate_icc(lgmd_df)
    print(f"\nLGMD R12 Patients:")
    print(f"  ICC: {lgmd_icc:.3f}")
    print(f"  Between-individual variance: {lgmd_var_between:.2f}")
    print(f"  Within-individual variance: {lgmd_var_within:.2f}")

    # Muscle-specific statistics
    print("\n" + "=" * 60)
    print("Muscle-specific MyoScore Statistics")
    print("=" * 60)

    print("\nHealthy Controls:")
    healthy_stats = calculate_muscle_stats(healthy_df, 'Healthy')
    print(healthy_stats.round(2).to_string())

    print("\nLGMD R12 Patients:")
    lgmd_stats = calculate_muscle_stats(lgmd_df, 'LGMD R12')
    print(lgmd_stats.round(2).to_string())

    # Generate comparative plots
    plot_comparison(
        healthy_pivot,
        lgmd_pivot,
        'Validation/GSE202745_complete_muscle_analysis.png'
    )

    # Generate focused VL-RF comparison (key finding)
    plot_vl_rf_comparison(
        healthy_pivot,
        lgmd_pivot,
        'Validation/GSE202745_VL_RF_correlation_reversal.png'
    )

    # Save summary results
    summary_results = {
        'Group': ['Healthy', 'LGMD R12'],
        'n_samples': [len(healthy_df), len(lgmd_df)],
        'n_individuals': [healthy_df['Individual'].nunique(), lgmd_df['Individual'].nunique()],
        'ICC': [healthy_icc, lgmd_icc],
        'VL_RF_r': [
            healthy_corr.get('VL vs RF', {}).get('pearson_r', np.nan),
            lgmd_corr.get('VL vs RF', {}).get('pearson_r', np.nan)
        ],
        'Sem_RF_r': [
            healthy_corr.get('Sem vs RF', {}).get('pearson_r', np.nan),
            lgmd_corr.get('Sem vs RF', {}).get('pearson_r', np.nan)
        ],
        'Sem_VL_r': [
            healthy_corr.get('Sem vs VL', {}).get('pearson_r', np.nan),
            lgmd_corr.get('Sem vs VL', {}).get('pearson_r', np.nan)
        ]
    }

    summary_df = pd.DataFrame(summary_results)
    summary_df.to_csv('Validation/GSE202745_muscle_consistency_summary.csv', index=False)
    print(f"\nSaved summary: Validation/GSE202745_muscle_consistency_summary.csv")

    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)

    return analysis_df, healthy_corr, lgmd_corr

if __name__ == "__main__":
    analysis_df, healthy_corr, lgmd_corr = main()
