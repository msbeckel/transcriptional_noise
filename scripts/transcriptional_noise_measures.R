# ========================================================  
# Script: TRANSCRIPTIONAL NOISE ANALYSIS - ENGE et al. METHODS  
# Author: MS Beckel  
# Date: 2025-03-20
# Last Updated: 2025-03-20
# ========================================================  
#
# This script implements methods from:
# Enge et al. (2017) Single-Cell Analysis of Human Pancreas 
# Reveals Transcriptional Signatures of Aging and Somatic Mutation Patterns
# ----------------------------------------------------------
#.libPaths("/ngs/mllorens/mbeckel/R_libs")

# REQUIRED PACKAGES
# ----------------------------------------------------------
library(Seurat)     # Single-cell data handling
library(dplyr)      # Metadata manipulation
library(matrixStats)# Efficient statistical computations
library(proxy)      # Distance calculations
library(progress)   # Progress bars
library(purrr)      # Functional programming
library(RcppAnnoy)

# 1. GENERAL HELPER FUNCTIONS
# ----------------------------------------------------------

#' @title Create Expression Bins
#' @description Divides genes into expression bins, excluding extremes
#' @param seurat Seurat object with normalized data and variable features
#' @return List of gene bins ordered by mean expression
#' @details
#' - Requires prior RunFindVariableFeatures()
#' - Creates 10 bins, excluding top/bottom 10% most/least expressed
#' - Output used for invariant gene detection in downstream analyses
get_bins <- function(seurat, assay = "RNA") {
  # 1. Get variable features and their expression means
  genes <- VariableFeatures(seurat)
  expr_data <- GetAssayData(seurat, assay = assay, slot = "data")[genes, ]
  means <- rowMeans(expr_data)
  
  # 2. Order genes by descending expression
  genes_ordered <- names(sort(means, decreasing = TRUE))
  
  # 3. Create 10 uniform bins
  n_bins <- 10
  bin_size <- floor(length(genes_ordered)/n_bins)
  
  # 4. Generate bins excluding extremes
  bins <- map(1:(n_bins-2), ~{
    start <- (.x-1)*bin_size + 1 + bin_size  # Skip first bin
    end <- start + bin_size - 1
    genes_ordered[start:end]
  })
  
  return(compact(bins))  # Remove empty bins if any
}

#' @title Identify Low-Variability Genes
#' @description Selects genes with lowest coefficient of variation in each bin
#' @param seurat Seurat object
#' @param bins Gene bins from get_bins()
#' @param assay Name of assay to use (default: "RNA")
#' @return Character vector of stable genes
#' @details
#' - Calculates coefficient of variation (CV = SD/mean) per gene
#' - Selects bottom 10% lowest CV genes in each bin
#' - Useful for finding expression-stable invariant genes
get_least_variable <- function(seurat, bins, assay = "RNA") {
  # Calculate CV per bin
  lvg <- map(bins, ~{
    # 1. Extract expression for current bin
    expr <- as.matrix(GetAssayData(seurat, assay = assay, slot = "data")[.x, ])
    
    # 2. Compute CV using optimized functions
    cv <- rowSds(expr)/rowMeans(expr)
    
    # 3. Select 10% most stable genes
    n_select <- floor(length(.x)*0.1)
    names(sort(cv))[1:n_select]
  }) %>% unlist()
  
  return(lvg)
}

# 2. CORE ANALYSIS METHODS
# ----------------------------------------------------------

#' @title Transcriptional Noise Calculation (Enge Method 1)
#' @description Computes biological vs technical variation ratio using ERCC controls
#' @param seurat Seurat object with expression data
#' @param batch Metadata column containing batch information
#' @param cell_type_col Metadata column with cell type annotations
#' @return Modified Seurat object with new metadata columns:
#' - cordist_bio: Correlation distance to biological mean
#' - cordist_tech: Correlation distance to technical mean (ERCC)
#' - noise: cordist_bio/cordist_tech ratio
#' @details
#' - Requires ERCC controls named with "ERCC-" prefix
#' - Computes group means by batch-cell_type combinations
#' - Uses correlation distance (1 - Pearson correlation)
enge_transcriptional_noise <- function(seurat, batch, cell_type_col = "cell_type") {
  # Data validation
  if(!cell_type_col %in% colnames(seurat@meta.data)) {
    stop("Cell type column missing: ", cell_type_col)
  }
  
  # 1. Identify ERCC controls
  ercc_genes <- grep("^ERCC-", rownames(seurat), value = TRUE)
  if(length(ercc_genes) == 0) stop("No ERCC controls detected")
  
  # 2. Prepare metadata grouping
  meta <- seurat@meta.data
  meta$group_id <- paste(meta[[batch]], meta[[cell_type_col]], sep = "_")
  
  # 3. Calculate group means
  bio_means <- AggregateExpression(seurat, 
                                   assays = "RNA",
                                   group.by = c(batch, cell_type_col))$RNA
  tech_means <- AggregateExpression(seurat, 
                                    assays = "RNA",
                                    features = ercc_genes,
                                    group.by = c(batch, cell_type_col))$RNA
  
  # 4. Compute distances for each cell
  pb <- progress_bar$new(total = ncol(seurat), format = "Calculating distances [:bar] :percent")
  
  seurat@meta.data <- seurat@meta.data %>% mutate(
    cordist_bio = map_dbl(colnames(seurat), ~{
      dist(rbind(
        GetAssayData(seurat, assay = "RNA")[, .x],
        bio_means[, meta[.x, "group_id"]]
      ), method = "correlation")
    }),
    cordist_tech = map_dbl(colnames(seurat), ~{
      dist(rbind(
        GetAssayData(seurat, assay = "RNA")[ercc_genes, .x],
        tech_means[, meta[.x, "group_id"]]
      ), method = "correlation")
    }),
    noise = cordist_bio / cordist_tech
  )
  
  return(seurat)
}

#' @title Distance to Cell Type Mean (Enge Method 2)
#' @description Computes distances to cell type-batch specific mean expression
#' @param seurat Seurat object
#' @param batch Metadata column with batch IDs
#' @param cell_type_col Metadata column with cell type IDs
#' @param assay Assay to use (default: "RNA")
#' @return Modified Seurat object with metadata columns:
#' - cor_dist: Correlation distance
#' - euc_dist: Euclidean distance
#' - man_dist: Manhattan distance
#' @details
#' - Calculates mean expression per batch-cell_type group
#' - Computes three distance metrics for each cell
#' - TODO: Check NA handling for small groups
distance_to_celltype_mean <- function(seurat, batch, cell_type_col = "cell_type", assay = "RNA") {
  # Input validation
  required_cols <- c(batch, cell_type_col)
  missing <- setdiff(required_cols, colnames(seurat@meta.data))
  if(length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))
  
  # 1. Calculate group means
  mean_expr <- AggregateExpression(seurat, 
                                   assays = assay,
                                   group.by = required_cols)[[assay]]
  colnames(mean_expr) <- gsub("-", "_", colnames(mean_expr))
  
  # 2. Initialize metadata columns
  seurat$cor_dist <- NA
  seurat$euc_dist <- NA
  seurat$man_dist <- NA
  
  # 3. Process each batch-cell_type combination
  groups <- unique(seurat@meta.data[, required_cols])
  pb <- progress_bar$new(total = nrow(groups))
  
  for(i in 1:nrow(groups)) {
    current_batch <- groups[i, batch]
    current_cell_type <- groups[i, cell_type_col]
    group_id <- paste(current_batch, current_cell_type, sep = "_")
    
    # 4. Subset cells
    cells <- rownames(seurat@meta.data)[
      seurat@meta.data[[batch]] == current_batch & 
        seurat@meta.data[[cell_type_col]] == current_cell_type
    ]
    
    # 5. Calculate distances per cell
    if(length(cells) > 0) {
      cell_data <- GetAssayData(seurat, assay = assay)[, cells, drop = FALSE]
      ref_vector <- mean_expr[, gsub("-", "_", group_id)]
      
      distances <- apply(cell_data, 2, function(x) {
        c(
          cor = dist(rbind(x, ref_vector), method = "correlation"),
          euc = dist(rbind(x, ref_vector), method = "euclidean"),
          man = dist(rbind(x, ref_vector), method = "manhattan")
        )
      })
      
      seurat@meta.data[cells, c("cor_dist", "euc_dist", "man_dist")] <- t(distances)
    }
    
    pb$tick()
  }
  
  return(seurat)
}

#' @title Invariant Gene Distance (Enge Method 3)
#' @description Computes Euclidean distance to global mean using invariant genes
#' @param seurat Seurat object
#' @param assay Assay to use (default: "RNA")
#' @return Modified Seurat object with metadata column:
#' - euc_dist_tissue_invar: Euclidean distance to invariant gene mean
#' @details
#' - Identifies stable genes using get_least_variable()
#' - Computes global mean of invariant genes
#' - Calculates per-cell Euclidean distance to this mean
enge_euclidean_dist <- function(seurat, assay = "RNA") {
  # 1. Identify invariant genes
  bins <- get_bins(seurat)
  invariant_genes <- get_least_variable(seurat, bins, assay = assay)
  
  # 2. Calculate global mean expression
  mean_expr <- rowMeans(GetAssayData(seurat, assay = assay)[invariant_genes, ])
  
  # 3. Compute distances
  seurat@meta.data$euc_dist_tissue_invar <- apply(
    GetAssayData(seurat, assay = assay)[invariant_genes, ], 
    2, 
    function(x) dist(rbind(x, mean_expr))
  )
  
  return(seurat)
}

# 3. GCL ANALYSIS FRAMEWORK
# ----------------------------------------------------------

#' @title GCL Calculation Wrapper
#' @description Computes Genuine Correlated Loading (GCL) across cell type-batch groups
#' @param seurat Seurat object
#' @param num_divisions Number of random gene splits
#' @param batch Metadata column with batch IDs
#' @param cell_type_col Metadata column with cell type IDs
#' @return Modified Seurat object with metadata column:
#' - GCL: Mean GCL score per cell group
#' @details
#' - Requires >10 cells per group
#' - Computes GCL for each batch-cell_type combination
#' - Returns mean GCL across iterations per group
gcl_per_cell_type_and_batch <- function(seurat, num_divisions, batch, cell_type_col = "cell_type", assay = "RNA") {
  # Validar columnas requeridas
  required_cols <- c(batch, cell_type_col)
  missing_cols <- setdiff(required_cols, colnames(seurat@meta.data))
  if(length(missing_cols) > 0) {
    stop("Columnas faltantes en metadata: ", paste(missing_cols, collapse = ", "))
  }
  # Inicializar columna GCL
  seurat@meta.data$GCL <- NA_real_
  
  # Obtener combinaciones únicas de cell_type y batch
  combinations <- seurat@meta.data %>%
    distinct(!!sym(cell_type_col), !!sym(batch)) %>%
    mutate(combo = paste(!!sym(cell_type_col), !!sym(batch), sep = "_"))
  pb <- progress_bar$new(
    format = "Procesando :current/:total [:bar] :percent (~:eta)",
    total = nrow(combinations)
  )
  # Iterar por cada combinación
  for(i in 1:nrow(combinations)) {
    ct <- combinations[i, cell_type_col]
    b <- combinations[i, batch]
    current_combo <- combinations$combo[i]
    
    # Filtrar células
    cells <- rownames(seurat@meta.data)[
      seurat@meta.data[[cell_type_col]] == ct & 
        seurat@meta.data[[batch]] == b
    ]
    
    if(length(cells) > 10) {
      subset <- subset(seurat, cells = cells)
      
      # Calcular GCL y asignar a metadata
      seurat@meta.data[cells, "GCL"] <- tryCatch({
        mean(gcl(subset, num_divisions, assay = "SCT"), na.rm = TRUE)  # Media de las divisiones
      }, error = function(e) {
        message("Error en ", current_combo, ": ", e$message)
        NA_real_
      })
    }
    
    pb$tick()
  }
  
  return(seurat)
}

#' @title Base GCL Calculation
#' @description Core implementation of Genuine Correlated Loading metric
#' @param seurat Seurat object subset
#' @param num_divisions Number of gene splits
#' @param assay Assay to use (default: "RNA")
#' @return Numeric vector of GCL values
#' @details
#' - Splits genes randomly into two halves
#' - Computes similarity between gene covariance patterns
#' - Repeat for num_divisions iterations
gcl <- function(seurat, num_divisions, assay = "RNA") {
  expr <- t(GetAssayData(seurat, assay = assay, slot = "data"))
  n_cells <- nrow(expr)
  gcl_values <- numeric(num_divisions)
  
  # Función interna para cálculo de matrices
  compute_Aij <- function(X) {
    d <- as.matrix(dist(X, method = "euclidean"))
    m <- rowMeans(d)
    M <- mean(d)
    Aij <- d - outer(m, m, "+") + M
    Aij - d/n_cells
  }
  
  for(i in 1:num_divisions) {
    # 1. Dividir genes aleatoriamente
    genes <- sample(ncol(expr))
    split <- floor(ncol(expr)/2)
    
    # 2. Calcular matrices para cada mitad
    A1 <- compute_Aij(as.matrix(expr[, genes[1:split]]))
    A2 <- compute_Aij(as.matrix(expr[, genes[(split+1):ncol(expr)]]))
    
    # 3. Calcular correlación entre matrices
    gcl_values[i] <- cor(c(A1), c(A2))
  }
  
  return(gcl_values)
}


# 4. ADDITIONAL METHODS
# ----------------------------------------------------------

#' @title Hernando-Herraez Correlation Analysis
#' @description Computes correlation to cell type-batch medians using HVGs
#' @param seurat Seurat object
#' @param batch Metadata column with batch IDs
#' @param cell_type_col Metadata column with cell type IDs
#' @param assay Assay to use (default: "RNA")
#' @return Modified Seurat object with metadata column:
#' - cor_dist_median: 1 - Spearman correlation to median expression
#' @details
#' - Uses top 500 variable genes
#' - Computes median expression per batch-cell_type group
#' - Calculates correlation distance for each cell
hernando_herraez <- function(seurat, batch, cell_type_col = "cell_type", assay = "RNA") {
  # Validación de columnas
  required_cols <- c(batch, cell_type_col)
  if (!all(required_cols %in% colnames(seurat@meta.data))) {
    stop("Columnas requeridas faltantes: ", paste(setdiff(required_cols, colnames(seurat@meta.data)), collapse = ", "))
  }
  
  # 1. Obtener genes variables
  hvgs <- VariableFeatures(seurat)
  if (length(hvgs) < 500) stop("Ejecutar FindVariableFeatures con nfeatures >= 500 primero")
  hvgs <- hvgs[1:500]
  
  # 2. Calcular medianas por grupo
  median_matrix <- AggregateExpression(
    seurat,
    assays = assay,
    features = hvgs,
    group.by = c(batch, cell_type_col),
    slot = "scale.data",
    fun = "median"
  )[[1]]
  colnames(median_matrix) <- gsub("-", "_", colnames(median_matrix))
  
  # 3. Mapeo de grupos
  group_ids <- gsub("-", "_", paste(seurat@meta.data[[batch]], seurat@meta.data[[cell_type_col]], sep = "_"))
  unique_groups <- colnames(median_matrix)
  
  # 4. Calcular distancias optimizado
  expr_data <- GetAssayData(seurat[[assay]], assay = "scale.data")[hvgs, ]
  
  # Función corregida
  calculate_distances <- function(cell_idx) {
    cell_expr <- expr_data[, cell_idx]
    group <- group_ids[cell_idx]
    group_median <- median_matrix[, group]
    if(cell_idx%%10000 == 0) cat(cell_idx, "\n")
    cor(
      cell_expr,
      group_median,
      method = "spearman", 
      use = "pairwise.complete.obs"
    )
  }
  
  # 5. Aplicación vectorizada
  seurat@meta.data$cor_dist_median <- 1 - sapply(
    seq_len(ncol(seurat)), 
    calculate_distances
  )
  
  return(seurat)
}

# 5. Own METHODS
# ----------------------------------------------------------

#' @title Own Correlation Analysis
#' @description Computes median euclidean distance to KNN neighbors for each cell.
#' @param seurat Seurat object
#' @param assay Assay to use (default: "RNA")
#' @return Modified Seurat object with metadata column:
#' - knn_noise: Median distance to KNN neighbors
#' @details
#' - TODO: for each batch?
#' 
knn_distance <- function(seurat, dims = 1:30, k = 10, assay = "SCT") {
  # 1. Find k-NN in PCA space
  seurat <- FindNeighbors(seurat, reduction = "pca", dims = dims, k.param = k, annoy.metric = "euclidean", assay = assay, return.neighbor = TRUE)
  
  # 2. Extract distances to k nearest neighbors
  knn_distances <- Distances(seurat@neighbors[[1]])  # Sparse matrix of distances
  
  # 3. Calculate median distance per cell
  median_noise <- apply(knn_distances, 1, function(x) median(x[x > 0]))  # Exclude self-distance (0)
  
  # 4. Add to metadata
  seurat$knn_distance <- median_noise
  
  return(seurat)
}
