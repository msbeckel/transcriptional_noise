# Transcriptional Noise Analysis with Mixed-Effects Models

## Overview
This repository contains an R pipeline to analyze age-associated changes in transcriptional noise across cell types using single-cell RNA sequencing (scRNA-seq) data. The analysis employs mixed-effects models to account for individual variability and provides robust statistical and visual outputs.

Dataset used from:

Buckley, M.T., Sun, E.D., George, B.M. et al. Cell-type-specific aging clocks to quantify aging and rejuvenation in neurogenic regions of the brain. Nat Aging 3, 121–137 (2023). [https://doi.org/10.1038/s43587-022-00335-4](https://doi.org/10.1038/s43587-022-00335-4)

---

## Objectives
1. **Model Transcriptional Noise**:  
   - Quantify the relationship between transcriptional noise and age for each cell type.
   - Compare log-transformed and rank-based modeling approaches.
2. **Generate Insights**:  
   - Identify cell types with significant age-dependent noise changes.
   - Provide interactive and publication-ready visualizations.

---

## Input
### Required Data
- **Seurat Object** (`multi_intergrated_seurat_Dec2020.rds`):  
  - **Metadata Columns**:  
    - `transcriptional_noise`: Precomputed transcriptional noise values (continuous, >0).  
    - `Age`: Numeric age values.  
    - `Celltype.LowRes`: Cell type annotations.  
    - `hash.id`: Subject identifiers (for random effects).  

### Optional
- **Cell Type Order**: A predefined order for visualization consistency (e.g., `["Microglia", "Oligodendro", ...]`).

---

## Output
### Results
- **CSV Files**:  
  - `transcriptional_noise_analysis_results.csv`: Model coefficients, confidence intervals, and FDR-adjusted p-values.  
  - `model_comparison.csv`: AIC comparisons between modeling approaches.  

### Diagnostics
- **Residual Plots** (PNG):  
  - QQ-plots and residual-vs-fitted plots for each cell type.  

### Visualizations
- **Static Plots** (PNG):  
  - Forest plots (effect sizes) and scatter plots (age vs. noise trends).  
- **Interactive Plots** (HTML):  
  - Hover-enabled versions of key plots for exploratory analysis.  

### Updated Data
- **Seurat Object** (`seu_updated.rds`):  
  - Metadata now includes model results and quality metrics.  

---

## Methods  
### Statistical Modeling  
1. **Linear Mixed-Effects Model (LMM)**:  
   - **Formula**: `log(transcriptional_noise) ~ Age + (1 | hash.id)`  
   - **Random Effects**: Subject-level intercepts (`hash.id`)  

2. **Model Validation**:  
   - **DHARMa Residual Diagnostics**:  
     - **Simulated Residuals**: Quantile residuals via parametric bootstrap (1000 simulations)  
     - **Key Checks**:  
       - Uniformity of residuals (QQ-plots)  
       - Residual dispersion  
       - Outlier detection  
       - Temporal/homoscedasticity patterns  
   - **Diagnostic Plots**:  
     - `plotQQunif()`: Uniformity assessment  
     - `plotResiduals()`: Residual vs. predicted values  
