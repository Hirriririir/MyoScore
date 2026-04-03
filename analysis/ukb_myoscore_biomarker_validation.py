# Note: UK Biobank data requires approved access (application 19542)
"""
UKB Validation of MyoScore Novel Genes via Biomarker-Muscle Associations.

Tests whether biomarkers corresponding to MyoScore novel genes (GGT7→GGT,
SLC44A2→Phosphatidylcholine, CYBRD1→indirect iron) associate with muscle
phenotypes in the expected direction.

Analysis 1: Biomarker ~ Muscle phenotype association (adjusted regression)
Analysis 2: Stratified analysis by biomarker quartiles
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# === Paths ===
UKB_TSV = Path("data/ukb_biomarker_data.tsv")  # UK Biobank data (application 19542)
OUT_DIR = Path("output/ukb_myoscore_validation")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# === Column mapping ===
COLS = {
    "eid": "eid",
    # Biomarkers (Tier 1)
    "p30730_i0": "GGT",           # GGT7 → GGT (U/L)
    "p23437_i0": "Phosphatidylcholine",  # SLC44A2 → PC (mmol/l)
    "p23438_i0": "Sphingomyelin",        # SLC44A2 → SM (mmol/l)
    "p23436_i0": "Total_Cholines",       # SLC44A2 → Total cholines (mmol/l)
    # Biomarkers (Tier 2)
    "p30880_i0": "Urate",         # AK9 → Urate
    "p30680_i0": "Calcium",       # CPNE1 → Calcium
    "p30080_i0": "Platelet_count",  # CBFA2T2 → Platelets
    "p30120_i0": "Lymphocyte_count",  # CBFA2T2 → Lymphocytes
    # Muscle phenotypes
    "p46_i0": "Grip_L",           # Grip strength left (Kg)
    "p47_i0": "Grip_R",           # Grip strength right (Kg)
    "p23101_i0": "FFM",           # Whole body fat-free mass (Kg)
    "p23113_i0": "Leg_FFM_R",     # Leg FFM right (Kg)
    "p23114_i0": "Leg_FFM_L",     # Leg FFM left (Kg) -- check
    "p924_i0": "Walk_pace",       # Usual walking pace (categorical)
    # Covariates
    "p21003_i0": "Age",
    "p31": "Sex",
    "p21001_i0": "BMI",
    "p21000_i0": "Ethnicity",
    "p54_i0": "Centre",
}

# === Step 1: Extract target columns from 33GB TSV ===
cache_path = OUT_DIR / "ukb_myoscore_subset.parquet"

if cache_path.exists():
    print(f"Loading cached subset from {cache_path}")
    df = pd.read_parquet(cache_path)
else:
    print(f"Extracting {len(COLS)} columns from UKB TSV (33GB)...")
    print("This may take a few minutes...")

    # Read header to get column indices
    with open(UKB_TSV, "r") as f:
        header = f.readline().strip().split("\t")

    col_indices = {}
    for ukb_col, name in COLS.items():
        if ukb_col in header:
            col_indices[header.index(ukb_col)] = name
        else:
            print(f"  WARNING: {ukb_col} ({name}) not found in header")

    print(f"  Found {len(col_indices)}/{len(COLS)} columns")

    # Use pandas usecols with column names for correct mapping
    ukb_col_names = list(COLS.keys())
    ukb_col_names_found = [c for c in ukb_col_names if c in header]
    print(f"  Reading {len(ukb_col_names_found)} columns via chunked reader...")
    chunks = []
    for chunk in pd.read_csv(
        UKB_TSV, sep="\t", usecols=ukb_col_names_found,
        low_memory=False, na_values=["", "NA", "nan"],
        chunksize=50000,
    ):
        chunks.append(chunk)
    df = pd.concat(chunks, ignore_index=True)
    print(f"  Read {len(df)} rows")
    # Rename columns
    df = df.rename(columns=COLS)

    # Convert numeric columns
    numeric_cols = [c for c in df.columns if c not in ("eid", "Walk_pace", "Sex", "Ethnicity", "Centre")]
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Sex: encode Female=0, Male=1
    df["Sex"] = df["Sex"].map({"Female": 0, "Male": 1, 0: 0, 1: 1})

    # Walk pace: 1=Slow, 2=Steady, 3=Brisk (treat as ordinal)
    df["Walk_pace"] = pd.to_numeric(df["Walk_pace"], errors="coerce")

    # Max grip strength
    df["Grip_max"] = df[["Grip_L", "Grip_R"]].max(axis=1)

    # Cache
    df.to_parquet(cache_path)
    print(f"  Cached to {cache_path}")

print(f"\nDataset: {len(df)} participants, {len(df.columns)} columns")
print(f"Non-null counts:")
for c in df.columns:
    n = df[c].notna().sum()
    if n > 0:
        print(f"  {c:25s}: {n:>7,d}")

# === Step 2: Compute max grip if not cached ===
if "Grip_max" not in df.columns:
    df["Grip_max"] = df[["Grip_L", "Grip_R"]].max(axis=1)

# === Step 3: Analysis 1 - Biomarker-Muscle Association (OLS regression) ===
print("\n" + "=" * 70)
print("Analysis 1: Biomarker-Muscle Phenotype Associations (adjusted OLS)")
print("=" * 70)

from sklearn.preprocessing import StandardScaler

biomarkers = {
    "GGT": {"gene": "GGT7", "direction": "↑H", "expected_sign": "negative",
             "rationale": "High GGT7 expr = healthy → high GGT = unhealthy → neg assoc with grip"},
    "Phosphatidylcholine": {"gene": "SLC44A2", "direction": "↑H", "expected_sign": "positive",
             "rationale": "High SLC44A2 expr = healthy → high PC = healthy → pos assoc with grip"},
    "Sphingomyelin": {"gene": "SLC44A2", "direction": "↑H", "expected_sign": "positive",
             "rationale": "High SLC44A2 expr = healthy → high SM = healthy → pos assoc with grip"},
    "Urate": {"gene": "AK9", "direction": "↑H", "expected_sign": "positive",
             "rationale": "High AK9 expr = healthy → moderate urate = antioxidant benefit"},
    "Calcium": {"gene": "CPNE1", "direction": "↑D", "expected_sign": "negative",
             "rationale": "High CPNE1 expr = disease → high Ca perturbation"},
    "Platelet_count": {"gene": "CBFA2T2", "direction": "↑H", "expected_sign": "positive",
             "rationale": "High CBFA2T2 expr = healthy → normal platelet function"},
}

muscle_outcomes = ["Grip_max", "FFM", "Grip_L", "Grip_R"]
covariates = ["Age", "Sex", "BMI"]

results = []

for biomarker, info in biomarkers.items():
    if biomarker not in df.columns:
        print(f"  SKIP: {biomarker} not in data")
        continue

    for outcome in muscle_outcomes:
        # Complete cases
        cols_needed = [biomarker, outcome] + covariates
        sub = df[cols_needed].dropna()

        if len(sub) < 100:
            continue

        # Standardize biomarker for comparable betas
        sub = sub.copy()
        sub[biomarker] = (sub[biomarker] - sub[biomarker].mean()) / sub[biomarker].std()

        # OLS: outcome ~ biomarker + age + sex + BMI
        import statsmodels.api as sm

        X = sub[[biomarker] + covariates]
        X = sm.add_constant(X)
        y = sub[outcome]

        model = sm.OLS(y, X).fit()

        beta = model.params[biomarker]
        se = model.bse[biomarker]
        p = model.pvalues[biomarker]
        ci_lo, ci_hi = model.conf_int().loc[biomarker]
        n = len(sub)

        # Check direction match
        expected = info["expected_sign"]
        observed = "positive" if beta > 0 else "negative"
        match = "YES" if expected == observed else "NO"

        results.append({
            "Gene": info["gene"],
            "Biomarker": biomarker,
            "Outcome": outcome,
            "N": n,
            "Beta": beta,
            "SE": se,
            "CI_lo": ci_lo,
            "CI_hi": ci_hi,
            "P": p,
            "Expected_direction": expected,
            "Observed_direction": observed,
            "Direction_match": match,
        })

        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"  {info['gene']:10s} | {biomarker:22s} → {outcome:10s} | "
              f"β={beta:+.3f} (SE={se:.3f}) | p={p:.2e} {sig} | "
              f"expect={expected:8s} got={observed:8s} [{match}]")

results_df = pd.DataFrame(results)
results_df.to_csv(OUT_DIR / "UKB_biomarker_muscle_associations.csv", index=False)
print(f"\nSaved: {OUT_DIR / 'UKB_biomarker_muscle_associations.csv'}")

# === Step 4: Analysis 2 - Stratified by biomarker quartiles ===
print("\n" + "=" * 70)
print("Analysis 2: Muscle phenotype by biomarker quartiles")
print("=" * 70)

tier1_biomarkers = ["GGT", "Phosphatidylcholine", "Sphingomyelin"]
tier1_available = [b for b in tier1_biomarkers if b in df.columns]

fig, axes = plt.subplots(len(tier1_available), 2, figsize=(14, 5 * len(tier1_available)))
if len(tier1_available) == 1:
    axes = axes.reshape(1, -1)

colors = ["#72c95e", "#31848f", "#46508b", "#50327b"]  # Q1(healthiest) to Q4

for row_i, biomarker in enumerate(tier1_available):
    info = biomarkers[biomarker]
    sub = df[[biomarker, "Grip_max", "FFM", "Age", "Sex", "BMI"]].dropna()
    sub = sub.copy()

    # Create quartiles
    sub["Quartile"] = pd.qcut(sub[biomarker].rank(method="first"), 4,
                               labels=["Q1", "Q2", "Q3", "Q4"])

    # Adjust for age, sex, BMI (residualize)
    import statsmodels.api as sm
    for outcome in ["Grip_max", "FFM"]:
        X_cov = sm.add_constant(sub[["Age", "Sex", "BMI"]])
        resid = sm.OLS(sub[outcome], X_cov).fit().resid
        sub[f"{outcome}_adj"] = resid + sub[outcome].mean()

    # Plot 1: Grip strength by quartile
    ax = axes[row_i, 0]
    quartile_data = []
    for qi, q in enumerate(["Q1", "Q2", "Q3", "Q4"]):
        vals = sub.loc[sub["Quartile"] == q, "Grip_max_adj"]
        quartile_data.append(vals)
        ax.bar(qi, vals.mean(), color=colors[qi], alpha=0.8,
               yerr=vals.sem(), capsize=5, edgecolor="white", linewidth=1)
        ax.text(qi, vals.mean() + vals.sem() + 0.2, f"n={len(vals):,}",
                ha="center", va="bottom", fontsize=8)

    # Trend test
    trend_r, trend_p = stats.spearmanr(
        sub[biomarker], sub["Grip_max_adj"]
    )

    ax.set_xticks(range(4))
    ax.set_xticklabels(["Q1\n(lowest)", "Q2", "Q3", "Q4\n(highest)"], fontsize=9)
    ax.set_ylabel("Adjusted Grip Strength (Kg)", fontsize=10)
    ax.set_title(f"{info['gene']} → {biomarker}\nvs Grip Strength (ρ={trend_r:.3f}, p={trend_p:.1e})",
                 fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Plot 2: FFM by quartile
    ax = axes[row_i, 1]
    for qi, q in enumerate(["Q1", "Q2", "Q3", "Q4"]):
        vals = sub.loc[sub["Quartile"] == q, "FFM_adj"]
        ax.bar(qi, vals.mean(), color=colors[qi], alpha=0.8,
               yerr=vals.sem(), capsize=5, edgecolor="white", linewidth=1)
        ax.text(qi, vals.mean() + vals.sem() + 0.1, f"n={len(vals):,}",
                ha="center", va="bottom", fontsize=8)

    trend_r2, trend_p2 = stats.spearmanr(sub[biomarker], sub["FFM_adj"])

    ax.set_xticks(range(4))
    ax.set_xticklabels(["Q1\n(lowest)", "Q2", "Q3", "Q4\n(highest)"], fontsize=9)
    ax.set_ylabel("Adjusted Fat-Free Mass (Kg)", fontsize=10)
    ax.set_title(f"{info['gene']} → {biomarker}\nvs FFM (ρ={trend_r2:.3f}, p={trend_p2:.1e})",
                 fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.suptitle("UKB Validation: MyoScore Novel Gene Biomarkers vs Muscle Phenotypes\n"
             "(adjusted for age, sex, BMI; bars = quartile means ± SEM)",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT_DIR / "UKB_biomarker_quartile_plot.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT_DIR / "UKB_biomarker_quartile_plot.pdf", bbox_inches="tight")
print(f"Saved: {OUT_DIR / 'UKB_biomarker_quartile_plot.png'}")

# === Step 5: Summary figure - Forest plot of all associations ===
print("\n" + "=" * 70)
print("Generating forest plot summary")
print("=" * 70)

grip_results = results_df[results_df["Outcome"] == "Grip_max"].copy()
grip_results = grip_results.sort_values("Beta")

fig, ax = plt.subplots(figsize=(10, max(5, len(grip_results) * 0.6)))

y_pos = range(len(grip_results))
for i, (_, row) in enumerate(grip_results.iterrows()):
    color = "#72c95e" if row["Direction_match"] == "YES" else "#e74c3c"
    ax.errorbar(row["Beta"], i, xerr=[[row["Beta"] - row["CI_lo"]],
                [row["CI_hi"] - row["Beta"]]], fmt="o",
                color=color, markersize=8, capsize=4, linewidth=2)

    sig = ""
    if row["P"] < 0.001:
        sig = "***"
    elif row["P"] < 0.01:
        sig = "**"
    elif row["P"] < 0.05:
        sig = "*"
    ax.text(row["CI_hi"] + 0.02, i, f"p={row['P']:.1e} {sig}",
            va="center", fontsize=8, color="gray")

ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8)
ax.set_yticks(y_pos)
labels = [f"{row['Gene']} → {row['Biomarker']}" for _, row in grip_results.iterrows()]
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("Standardized β (per SD biomarker → Grip Strength, Kg)", fontsize=11)
ax.set_title("UKB Validation: MyoScore Novel Gene Biomarker → Grip Strength\n"
             "(green = direction matches MyoScore prediction, red = mismatch)",
             fontsize=12, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(OUT_DIR / "UKB_forest_plot.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT_DIR / "UKB_forest_plot.pdf", bbox_inches="tight")
print(f"Saved: {OUT_DIR / 'UKB_forest_plot.png'}")

# === Print summary ===
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
grip_only = results_df[results_df["Outcome"] == "Grip_max"]
n_match = (grip_only["Direction_match"] == "YES").sum()
n_total = len(grip_only)
n_sig = (grip_only["P"] < 0.05).sum()
print(f"Biomarkers tested: {n_total}")
print(f"Significant (p<0.05): {n_sig}/{n_total}")
print(f"Direction concordant: {n_match}/{n_total}")
print(f"Direction + significant: {((grip_only['Direction_match']=='YES') & (grip_only['P']<0.05)).sum()}/{n_total}")
