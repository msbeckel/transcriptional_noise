#!/usr/bin/env python3
"""
run_decibel_noise.py

Loads a preprocessed AnnData (.h5ad), runs Decibel noise metrics including:
 - Enge transcriptional noise (requires ERCC spike-ins)
 - Distance-to-celltype-mean (euclidean, correlation, manhattan)
 - Scallop pipeline (membership stability)
Saves updated AnnData with noise annotations.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import argparse
import scanpy as sc
import scallop as sl
import sys
sys.path.append('/home/mllorens/mbeckel/repos_ext/decibel/module/')
import decibel as dcb
import scipy.sparse as sp


def run_decibel_noise(adata_path: str, results_dir: str, plots_dir: str, batch: str = 'batch', cell_type: str = 'cell_type', condition: str = 'condition', num_divisions: int = 10):
    # 1) Load data
    adata = sc.read_h5ad(adata_path, as_sparse="X")  
    print(f"Loaded AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

    # 2) Preprocess data
    adata.obs["batch"] = adata.obs[batch]
    adata.obs["cell_type"] = adata.obs[cell_type]
    adata.obs["condition"] = adata.obs[condition]

    # 2.1) Group clusters by marker genes
    print("Regrouping clusters by marker genes…")
    cluster_to_cell_type = {
        "0": "DG",        # Cluster 0 → DG
        "7": "DG",        # Cluster 7 → DG
        "4": "CA3",       # Cluster 4 → CA3
        "2": "inhibitory",# Cluster 2 → inhibitory
        "1": "CA1",       # Cluster 1 → CA1
        "3": "CA1",       # Cluster 3 → CA1
        "5": "CA1",       # Cluster 5 → CA1
        "6": "CA1"        # Cluster 6 → CA1
    }

    # Ensure Leiden clusters are strings (if stored as integers)
    adata.obs['leiden'] = adata.obs['leiden'].astype(str)

    # Map clusters to regions and create a new column
    adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_cell_type)  

    print("Marker genes added to obs.")

    # 3) Distance to cell-type mean (three metrics)
    print("Computing distance to cell-type mean…")
    dcb.distance_to_celltype_mean(adata,batch="batch")

    # 3.1)
    print("Computing distance to cell-type mean with invariant genes…")
    dcb.distance_to_celltype_mean_invariant(adata,batch="batch")

    # 4) Scallop pipeline: membership stability → noise = 1 - mean membership
    print("Running Scallop pipeline for stability noise…")
    dcb.scallop_pipeline(adata) 

    # 5) Optionally, compute GCL (global coordination level)
    print("Computing global coordination level (GCL)…")
    gcl_results = dcb.gcl_per_cell_type_and_batch(adata, batch="batch", num_divisions=10)

    # 5.1) Save GCL results to CSV
    gcl_out_path = os.path.join(results_dir, "gcl_results.csv")
    gcl_results.to_csv(gcl_out_path)

    # Aggregate results (e.g., mean GCL per cell type and batch)
    gcl_agg = (
        gcl_results
        .groupby(["cell_type", "condition", "batch"])
        .agg(mean_gcl=("GCL", "mean"))
        .reset_index()
    )
    print(gcl_agg)

    # Merge with adata.obs to map GCL values to cells
    adata.obs = adata.obs.merge(
        gcl_agg,
        how="left",
        on=["cell_type", "condition", "batch"]
    )
    adata.obs["mean_gcl"] = adata.obs["mean_gcl"].astype(float)

    #Plots
    print("📊 Ploting…")
    tn_plot_dir = os.path.join(plots_dir, "noise_measures")
    os.makedirs(tn_plot_dir, exist_ok=True)
    sc.settings.figdir = tn_plot_dir
    # UMAP plot
    sc.pl.umap(
        adata,
        color=['cell_type', 'condition', 'genotype', 'scallop_noise', 'mean_gcl', 'cor_dist', 'euc_dist', 'cor_dist_invar', 'euc_dist_invar'],
        wspace=0.5,
        ncols=3, 
        save="_transcriptional_noise.png",
        show=False
    )

    # UMAP genes plot
    sc.pl.umap(
        adata,
        color=['cell_type', 'condition', 'genotype', 'scallop_noise', 'mean_gcl', 'Prox1', 'Mpped1', 'Mndal', 'Gad1'],
        wspace=0.5,
        ncols=3, 
        save="_marker_genes.png",
        show=False
    )

    # Violin plots (all in one)
    sc.pl.violin(
        adata,
        ['scallop_noise', 'mean_gcl', 'cor_dist', 'euc_dist', 'man_dist', 'cor_dist_invar', 'euc_dist_invar', 'man_dist_invar'],
        groupby='cell_type',
        jitter=0.4,
        multi_panel=True,
        save="_violin_noise_metrics.png",
        show=False
    )

    print("Na values in noise metrics:")
    print(f"Scallop: {adata.obs['scallop_noise'].isna().sum()}")
    print(f"Mean GCL: {adata.obs['mean_gcl'].isna().sum()}")
    print(f"Correlation distance: {adata.obs['cor_dist'].isna().sum()}")


    # 6) Save results
    os.makedirs(results_dir, exist_ok=True)
    out_file = os.path.join(results_dir, "adata_with_decibel_noise.h5ad")
    adata.write(out_file)
    print(f"Saved annotated AnnData to {out_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Decibel transcriptional noise estimation on an AnnData"
    )
    parser.add_argument("--adata", required=True, help="Input preprocessed .h5ad file")
    parser.add_argument("--results_dir", default="decibel_results", help="Directory for output .h5ad")
    parser.add_argument("--plots_dir", default="noise_plots", help="Directory for output plots")
    parser.add_argument("--batch", default="batch", help="obs key for batch/donor")
    parser.add_argument("--cell_type", default="cell_type", help="obs key for cell type")
    parser.add_argument("--condition", default="condition", help="obs key for condition")
    parser.add_argument("--num_divisions", default=10, help="num divisions for GCL")
    args = parser.parse_args()

    run_decibel_noise(args.adata, args.results_dir, args.plots_dir, args.batch, args.cell_type, args.condition, args.num_divisions)
    print("Finished running Decibel noise estimation.")
    print("Results saved to:", args.results_dir)
    print("Done.")
