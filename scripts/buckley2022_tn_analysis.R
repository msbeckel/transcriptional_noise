# ========================================================  
# Script: TRANSCRIPTIONAL NOISE ANALYSIS - Buckley et al. 2022 dataset  
# Author: MS Beckel  
# Date: 2025-04-11
# Last Updated: 2025-04-11
# ========================================================  
#.libPaths("/ngs/mllorens/mbeckel/R_libs")

# REQUIRED PACKAGES
# ----------------------------------------------------------
pacman::p_load(
  Seurat,
  dplyr,
  matrixStats,
  proxy,
  progress,
  purrr,
  RcppAnnoy,
  lme4,
  ggplot2,
  plotly,
  MuMIn,
  broom.mixed,
  DHARMa)
library(Seurat)
# Set working directories
# ----------------------------------------------------------
SCRIPT_DIR <- "/home/mllorens/mbeckel/projects/transcriptional_noise/scripts/"
INPUT_DIR  <- "/ngs/mllorens/mbeckel/transcriptional_noise/data/01_processed/"
OUTPUT_DIR <- "/ngs/mllorens/mbeckel/transcriptional_noise/results/"

# Load custom functions
# ----------------------------------------------------------
source(paste0(SCRIPT_DIR, "transcriptional_noise_measures.R"))
source(paste0(SCRIPT_DIR, "glmm_models.R"))

# Load seu_subrat object
seu <- readRDS(file = paste0(INPUT_DIR, "multi_intergrated_seurat_Dec2020.rds"))
seu <- subset(x = seu, idents = names(sort(table(seu$Celltype.LowRes), decreasing = T)[1:5]))
seu$Age_biological <- (35 - seu$Prolif_Lineage_Fraction_of_SVZ*100)

# 2. Run Enge analyses
seu <- distance_to_celltype_mean(seu, batch = "hash.ID", cell_type_col = "Celltype.LowRes", assay = "SCT")
seu <- enge_euclidean_dist(seurat = seu, assay = "SCT")
seu <- gcl_per_cell_type_and_batch(seurat = seu, num_divisions = 10, batch = "hash.ID", cell_type_col = "Celltype.LowRes", assay = "SCT")
seu <- hernando_herraez(seurat = seu, batch = "hash.ID", cell_type_col = "Celltype.LowRes", assay = "SCT")
seu <- knn_distance(seurat = seu, dims = 1:30, k = 10, assay = "SCT")
#
# 3. Access results
# head(seurat@meta.data)
# FeaturePlot(seurat, features = "noise")
# VlnPlot(seurat, features = "noise", group.by = "Celltype.LowRes")


# 1. Análisis con distribución Gamma
gamma_results <- run_glmm_noise_analysis(
  seurat = seu,
  response_var = "knn_distance",
  output_dir = OUTPUT_DIR,
  fixed_effect = "Age",
  random_effect = "hash.ID",
  family = Gamma(link = "log")
)

# 1. Análisis con distribución Gamma
gamma_results_invar <- run_glmm_noise_analysis(
  seurat = seu,
  response_var = "euc_dist_tissue_invar",
  output_dir = OUTPUT_DIR,
  fixed_effect = "Age",
  random_effect = "hash.ID",
  family = Gamma(link = "log")
)

# 1. Análisis con distribución Gamma
gamma_results_cor_dist_median <- run_glmm_noise_analysis(
  seurat = seu,
  response_var = "cor_dist_median",
  output_dir = OUTPUT_DIR,
  fixed_effect = "Age",
  random_effect = "hash.ID",
  family = Gamma(link = "log")
)

# 1. Análisis con distribución Gamma
gamma_results_cor_dist <- run_glmm_noise_analysis(
  seurat = seu,
  response_var = "cor_dist",
  output_dir = OUTPUT_DIR,
  fixed_effect = "Age",
  random_effect = "hash.ID",
  family = Gamma(link = "log")
)

# 1. Análisis con distribución Gamma
gamma_results_euc_dist <- run_glmm_noise_analysis(
  seurat = seu,
  response_var = "euc_dist",
  output_dir = OUTPUT_DIR,
  fixed_effect = "Age",
  random_effect = "hash.ID",
  family = Gamma(link = "log")
)

# 1. Análisis con distribución Gamma
gamma_results_man_dist <- run_glmm_noise_analysis(
  seurat = seu,
  response_var = "man_dist",
  output_dir = OUTPUT_DIR,
  fixed_effect = "Age",
  random_effect = "hash.ID",
  family = Gamma(link = "log")
)


# 1. Análisis con distribución Gamma
seu[["GCL_norm"]] <- scale(seu[["GCL"]][, 1], scale = TRUE, center = FALSE)
gamma_results_GCL_norm <- run_glmm_noise_analysis(
  seurat = seu,
  response_var = "GCL_norm",
  output_dir = OUTPUT_DIR,
  fixed_effect = "Age",
  random_effect = "hash.ID",
  family = Gamma(link = "log")
)

# 1. Análisis con distribución Gamma
gamma_results_Age_bio <- run_glmm_noise_analysis(
  seurat = seu,
  response_var = "Age_biological",
  output_dir = OUTPUT_DIR,
  fixed_effect = "Age",
  random_effect = "hash.ID",
  family = Gamma(link = "log")
)

# 1. Análisis con distribución Gamma
gamma_results_transcriptional_noise <- run_glmm_noise_analysis(
  seurat = seu,
  response_var = "transcriptional_noise",
  output_dir = OUTPUT_DIR,
  fixed_effect = "Age",
  random_effect = "hash.ID",
  family = Gamma(link = "log")
)

