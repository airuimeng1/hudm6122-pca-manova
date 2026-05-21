# HUDM 6122 Final Project — PCA + Factorial MANOVA on China's Provincial Digital Economy

**Measuring Regional Disparities in China's Digital Economy Development:
A Multivariate Analysis Using PCA and Factorial MANOVA, 2013–2022**

Airui Meng · HUDM 6122 Multivariate Analysis I · Teachers College, Columbia University · Spring 2026

This repository contains the R analysis pipeline, panel data, intermediate outputs, and
typeset report for a course project that constructs a Digital Economy Development Index
(DEDI) for 31 mainland Chinese provinces, 2013–2022, and tests for regional and temporal
differences via a 4 × 2 factorial MANOVA on retained principal-component scores.

---

## Headline results

- **PCA**: 4 components retained (Kaiser + Cumulative ≥ 80% rule), explaining 85.8% of
  variance after log-transforming heavy-skew indicators and Varimax-rotating the retained
  loadings. KMO = 0.904 ("meritorious"); Bartlett χ²(231) = 11553.82, *p* < .001.
- **Factorial MANOVA** (DVs = PC1–PC4, *n* = 310 province-years):
  - Region main effect: Pillai *V* = 0.754, *F*(12, 903) = 25.25, *p* < .001, partial η² = 0.25.
  - Period main effect: Pillai *V* = 0.719, *F*(4, 299) = 191.12, *p* < .001, partial η² = 0.72.
  - Region × Period interaction: **non-significant** (*V* = 0.021, *p* = .892, partial η² = 0.007).
- **Substantive interpretation**: every region advanced markedly between 2013–2017 and
  2018–2022, but the relative rank-ordering of regions on a four-dimensional digital-economy
  construct has not changed — *σ*-divergence in absolute spread alongside
  *β*-convergence in relative spread.
- Sensitivity analyses (single-year MANOVA on 2018/2019/2020, year-by-year PCA
  loading stability, and outlier-removed refit) corroborate the main conclusions.

---

## Repository structure

```
hudm6122-pca-manova/
├── README.md                       This file.
├── data/
│   ├── panel_data.xlsx              Cleaned province × year panel, 22 indicators.
│   ├── original_data.xlsx           Raw indicators before interpolation.
│   ├── interpolation_processing.xlsx
│   └── indicator_system.xlsx        3-tier indicator hierarchy.
├── R/
│   ├── 00_setup.R                   Loaders, region map, indicator labels.
│   ├── 01_descriptive.R             Descriptive statistics + missing-value audit.
│   ├── 02_pca.R                     log1p + z-score + KMO/Bartlett + PCA + Varimax.
│   ├── 03_dedi.R                    Variance-weighted DEDI composite (0–100 scale).
│   ├── 04_assumptions.R             Mardia, Box's M, Levene, Mahalanobis diagnostics.
│   ├── 05_manova.R                  Factorial MANOVA + univariate ANOVAs + Games–Howell.
│   └── 06_sensitivity.R             Single-year MANOVA, PCA stability, outlier-drop refit.
├── outputs/
│   ├── figures/                     8 PNG figures at 300 dpi.
│   ├── tables/                      16 CSV tables.
│   ├── pc_scores.csv                PC1–PC4 per province-year.
│   ├── panel_dedi.csv               PC1–PC4 + DEDI per province-year.
│   ├── pca_objects.rds              prcomp(), Varimax rotation, cached objects.
│   ├── manova_results.rds           lm() fit, Manova(), Games-Howell tables.
│   └── assumption_checks.rds        Mardia, Box's M, Levene, Mahalanobis.
└── report/
    ├── final_report.pdf             Typeset write-up (17 pages).
    ├── final_report.tex             LaTeX source.
    └── *.png                        Figures referenced by the report.
```

---

## How to reproduce

The pipeline is staged so that each step writes intermediate artefacts that the next
step reads back. Run from the repository root:

```r
setwd("path/to/hudm6122-pca-manova")

## Stage 0: shared setup (loaded automatically by each downstream script)
source("R/00_setup.R")

## Stage 1: descriptive statistics + imputation audit
source("R/01_descriptive.R")

## Stage 2: PCA — writes pc_scores.csv, pca_objects.rds, scree/cumvar/loadings figures
source("R/02_pca.R")

## Stage 3: DEDI composite + regional/temporal heatmap and trend figures
source("R/03_dedi.R")

## Stage 4: MANOVA assumption diagnostics
source("R/04_assumptions.R")

## Stage 5: factorial MANOVA + Games–Howell post hoc + interaction plot
source("R/05_manova.R")

## Stage 6: sensitivity analyses
source("R/06_sensitivity.R")
```

### Dependencies

R ≥ 4.3 with the following CRAN packages:

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "ggplot2",
  "psych",       # KMO, Bartlett, skewness
  "car",         # Manova, leveneTest, Anova
  "MVN",         # Mardia tests
  "rstatix",     # games_howell_test
  "heplots",     # boxM, etasq
  "knitr", "scales"
))
```

To re-typeset the report PDF, `cd report/` and run `pdflatex final_report.tex` twice
(any TeX Live installation with `amsmath`, `booktabs`, `tabularx`, `multirow`,
`ragged2e`, `subcaption`, `natbib`, and `hyperref` will do).

---

## Method at a glance

The 22 third-level indicators are nested under three first-level dimensions (digital
infrastructure, digital industry, digital environment), following the framework of
He et al.\ (2023). Indicators with |skew| > 1.5 are log-transformed before *z*-score
standardization. Four PCs are retained by the conservative
*k* = max(*k*\_Kaiser, *k*\_80%) rule and Varimax-rotated for interpretability;
the four PCs are labelled **Industry/Innovation Scale**, **Per-capita Penetration**,
**Infrastructure Connectivity**, and **Enterprise Digitalization**.

A composite DEDI index is constructed as the variance-share-weighted, sign-aligned sum
of the retained scores, rescaled to [0, 100], purely for descriptive ranking.
**Inferential** analysis operates directly on the four orthogonal PC scores — the
multicollinearity assumption of MANOVA is therefore automatically satisfied.

The factorial MANOVA fits

$$
\mathbf{Y}_{ij\ell} = \boldsymbol\mu + \boldsymbol\alpha_i + \boldsymbol\beta_j +
(\boldsymbol{\alpha\beta})_{ij} + \boldsymbol\varepsilon_{ij\ell}
$$

with *i* ∈ {East, Central, West, Northeast}, *j* ∈ {Early 2013–17, Late 2018–22},
and *ℓ* indexing province-years within each cell. Because Box's *M* rejects covariance
homogeneity, **Pillai's trace** is reported as the primary multivariate criterion;
Wilks' Λ, Hotelling–Lawley *T*², and Roy's largest root are reported for transparency.
Univariate Type-II ANOVAs decompose the multivariate effect onto each PC, and pairwise
region contrasts use the **Games–Howell** procedure (which does not assume equal
variances or equal cell sizes).

---

## Data sources

- *China Statistical Yearbook* (National Bureau of Statistics)
- *Statistical Yearbook of China's Tertiary Industry*
- *China Science and Technology Statistical Yearbook*
- *Digital Inclusive Finance Index* (Digital Finance Research Center, Peking University)

All for fiscal years 2013–2022; 31 mainland provinces (HK / Macao / Taiwan excluded).
Missing cells were filled by linear interpolation in the time dimension; the imputation
footprint is small (only `x5` IPv4 addresses at 20.3% of cells, `x13` mobile internet
users at 8.4%, and `x21` R&D institutions at 0.3% required any interpolation).

---

## Licence and citation

This is coursework submitted in partial fulfilment of HUDM 6122 (Spring 2026) at
Teachers College, Columbia University. The code and analysis are released for
review and reproduction. The underlying yearbook data are public statistical
releases of the Chinese government; the Digital Inclusive Finance Index belongs to
the Digital Finance Research Center of Peking University.

If you build on this work, please cite it as:

> Meng, Airui (2026). *Measuring Regional Disparities in China's Digital Economy
> Development: A Multivariate Analysis Using PCA and Factorial MANOVA,
> 2013–2022.* HUDM 6122 Final Project, Teachers College, Columbia University.
