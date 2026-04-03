# MyoScore <img src="man/figures/logo.svg" align="right" height="139" />

**A genetically-informed transcriptomic scoring system for quantifying human skeletal muscle health**

MyoScore quantifies skeletal muscle health across five genetically-driven
dimensions based on GWAS-TWAS integration of 27 muscle-related phenotypes.

## Installation

```r
# Install from CRAN
install.packages("MyoScore")

# or Install from GitHub
devtools::install_github("Hirriririir/MyoScore")
```

## Quick Start

```r
library(MyoScore)

# Calculate MyoScore from raw count matrix
scores <- myoscore_score("path/to/raw_counts.csv")

# Or from an R matrix
scores <- myoscore_score(count_matrix)

# View results
head(scores)
#>     Strength_score Mass_score LeanMuscle_score Youth_score Resilience_score MyoScore
#> S1          72.3       65.1             80.2        55.8             68.4     69.2
#> S2          45.1       38.7             42.3        61.2             35.6     44.1
```

## Five Dimensions

| Dimension | Weight | GWAS Basis |
|-----------|--------|------------|
| **Strength** | 25.2% | Grip strength, walking pace |
| **Mass** | 17.7% | Fat-free mass (whole body, limbs) |
| **LeanMuscle** | 24.3% | Thigh fat infiltration MRI |
| **Youth** | 24.2% | Telomere length |
| **Resilience** | 8.7% | Myopathy diagnosis, CK levels |

**Higher score = healthier muscle** (0-100 scale)

## Visualization

```r
# Radar chart (requires fmsb)
myoscore_plot_radar(scores, groups = metadata$condition)

# Grouped boxplot (requires ggplot2)
myoscore_plot_boxplot(scores, groups = metadata$condition)
```

## Input Requirements

- Raw count matrix (genes x samples), **not** TPM/FPKM
- Gene Symbols as row names
- At least 2 samples (recommend >= 20)

## Citation

> Revealing myopathy spectrum: integrating transcriptional and clinical
> features of human skeletal muscles with varying health conditions.
> *Communications Biology*, 2024.
> DOI: [10.1038/s42003-024-06096-7](https://doi.org/10.1038/s42003-024-06096-7)

## License

MIT
