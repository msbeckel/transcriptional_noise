#!/usr/bin/env python3
"""
preprocess_scrna.py

Reads 10× Genomics MEX files for GEO datasets, filters, normalizes, log-transforms,
calculates mitochondrial gene metrics, and saves QC plots.
"""

import os
import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData
import matplotlib.pyplot as plt

def add_barcode_metadata(adata: AnnData) -> AnnData:
    """
    Splits barcodes in adata.obs_names by '_' and adds parts and sample_id to .obs.

    Expected format: barcode_mouse_conditionRep_lane (e.g., AACGTGAT_M6_APP1_L1)

    - Adds: barcode, mouse, condition, replicate, lane
    - Adds: sample_id = mouse_condition_replicate_lane

    Args:
        adata: AnnData object

    Returns:
        Modified AnnData with new metadata columns in .obs
    """
    parts = [bc.split('_') for bc in adata.obs_names]

    # Ensure all barcodes have exactly 4 parts
    for i, p in enumerate(parts):
        if len(p) != 4:
            raise ValueError(f"Barcode at index {i} does not have 4 parts: {p}")

    # Create initial DataFrame
    df = pd.DataFrame(parts, index=adata.obs_names, columns=["barcode", "age", "condition_full", "lane"])

    # Split condition_full into genotype and replicate using regex
    df[["genotype", "replicate"]] = df["condition_full"].str.extract(r"([A-Za-z]+)(\d+)", expand=True)

    # Build sample_id from age, genotype, replicate, and lane
    df["sample_id"] = df[["age", "genotype", "replicate"]].agg('_'.join, axis=1)

    # Drop the original condition_full column
    #df.drop(columns=["condition_full"], inplace=True)

    # Merge with AnnData.obs
    adata.obs = pd.concat([adata.obs, df.drop(columns="condition_full")], axis=1)

    # Convert all columns in adata.obs to categorical dtype
    for col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype('category')

    return adata

def save_qc_plots(adata, results_dir):
    """
    Saves specific QC plots to the results_dir/qc_plots directory.
    """
    qc_plot_dir = os.path.join(results_dir, "qc_plots")
    os.makedirs(qc_plot_dir, exist_ok=True)
    sc.settings.figdir = qc_plot_dir

    # Violin plots (all in one)
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save="_violin_qc_metrics.png",
        show=False
    )

    # Scatter: total counts vs pct MT
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        save="_scatter_total_vs_mt.png",
        show=False
    )

    # Scatter: total counts vs gene counts
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        save="_scatter_total_vs_genes.png",
        show=False
    )

    # Plot doublet score distribution
    sc.pl.scrublet_score_distribution(adata,   
               save="_doublet_score_hist.png", 
               show=False
    )
    
    # Highly variable genes
    sc.pl.highly_variable_genes(
        adata,
        save="_HVG.png",
        show=False
    )

    # Top 20 highly expressed genes
    sc.pl.highest_expr_genes(adata, 
                             n_top=20, 
                             save="_HVG_top20.png", 
                             show=False
    )

    # PCA plot
    sc.pl.pca_variance_ratio(adata, 
                             n_pcs=50,
                             save="_pca_variance_ratio.png",
                             log=True, 
                             show=False)

    # PCA plot PC1 vs PC2 and PC2vs PC3
    sc.pl.pca( adata,
            color=["leiden","leiden","pct_counts_mt", "pct_counts_mt"],
            dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
            ncols=2,
            size=2,
            save="_pca.png",
            show=False
    )

    # UMAP plot
    sc.pl.umap(
        adata,
        color=["leiden", "predicted_doublet", "log1p_total_counts", "pct_counts_mt"],#, "log1p_n_genes_by_counts"],
        wspace=0.5,
        ncols=2, 
        save="_umap.png",
        show=False
    )

    print(f"✅ All QC plots saved to {qc_plot_dir}")

def preprocess_10x_mtx(input_dir: str, output_h5ad: str, results_dir: str, HVG: int = 2000):
    """
    Preprocesses 10x MEX files. Attempts standard load, falls back to matrix search.
    Calculates mitochondrial QC, filters, normalizes, logs, saves AnnData and plots.
    """
    folder_name = os.path.basename(os.path.normpath(input_dir))
    prefix = f"{folder_name}_"

    print(f"Detected prefix '{prefix}' from folder '{folder_name}'")
    print(f"Reading 10× data with prefix '{prefix}' from {input_dir}…")

    try:
        adata = sc.read_10x_mtx(
            input_dir,
            prefix=prefix,
            cache=False
        )
        print("✅ Successfully read using sc.read_10x_mtx.")

    except Exception as e:
        print(f"⚠️  sc.read_10x_mtx failed: {e}")
        print("🔍 Searching for a file with 'matrix' in its name…")

        matrix_file = next((f for f in os.listdir(input_dir) if "matrix" in f.lower()), None)
        if not matrix_file:
            raise FileNotFoundError("No file with 'matrix' in its name found in input_dir for fallback.")

        mtx_path = os.path.join(input_dir, matrix_file)
        print(f"📄 Found matrix file: {matrix_file}. Reading with sc.read_csv()")
        adata = sc.read_csv(mtx_path, delimiter=" ")
        print("✅ Successfully created AnnData from raw matrix.")


    #Check if the matrix is transposed
    sample_like = adata.var_names.str.match(r"^[ACGTN]{8,}").sum()
    if sample_like > len(adata.var_names) / 2:
        print("Matrix appears transposed. Transposing to match expected format…")
        adata = adata.T

    print(f"  Loaded {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Add metadata from barcodes
    print("🔍 Adding metadata from barcodes…")
    adata = add_barcode_metadata(adata)

    # Annotate mitochondrial genes
    adata.var_names_make_unique()
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("Mt")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    print("🧬 Calculating QC metrics including mitochondrial percentage…")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    # QC filters
    print("📉 Filtering cells with <200 genes…")
    sc.pp.filter_cells(adata, min_genes=200)
    print(f"  {adata.n_obs} cells remain")

    print("📉 Filtering genes expressed in <3 cells…")
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  {adata.n_vars} genes remain")

    # Predict cell doublets 
    print("🔎 Predicting and Filtering Doublets")
    sc.pp.scrublet(adata)
    # Filter out predicted doublets
    #adata = adata[adata.obs['predicted_doublet'] == False].copy()

    # Saving count data
    adata.layers["counts"] = adata.X.copy()

    # Normalize and log-transform
    print("📊 Normalizing to 10 000 counts per cell…")
    sc.pp.normalize_total(adata, target_sum=1e4)

    print("🔎 Applying log1p transform…")
    sc.pp.log1p(adata)

    # Highly variable genes
    print("🔍 Identifying highly variable genes…")
    sc.pp.highly_variable_genes(adata, n_top_genes=HVG, subset=False)
    
    #Dimensionality reduction
    print("🔎 Performing PCA…")
    sc.tl.pca(adata)

    print("🔍 Performing UMAP…")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    print("🔎 Performing tSNE…" )
    sc.tl.tsne(adata, n_pcs=20)
    print("🔍 Performing Leiden clustering…" )
    sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2)

    # Save
    os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
    adata.write(output_h5ad)
    print(f"💾 Preprocessed data saved to {output_h5ad}")

    # Save plots
    print(f"📈 Saving QC plots to {results_dir}/qc_plots")
    save_qc_plots(adata, results_dir)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Preprocess 10× MEX files into an AnnData."
    )
    parser.add_argument(
        "--input_dir", required=True,
        help="Directory containing the 10× files (matrix, genes, barcodes)"
    )
    parser.add_argument(
        "--output_h5ad", required=True,
        help="Path to write the preprocessed AnnData (.h5ad)"
    )
    parser.add_argument(
        "--results_dir", required=True,
        help="Directory to save quality control plots"
    )
    args = parser.parse_args()

    preprocess_10x_mtx(args.input_dir, args.output_h5ad, args.results_dir)
