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


def run_decibel_noise(adata_path: str, results_dir: str, batch: str = 'batch', cell_type: str = 'cell_type', condition: str = 'condition'):
    # 1) Load data
    adata = sc.read_h5ad(adata_path)  
    print(f"Loaded AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

    # 2) Preprocess data
    adata["batch"] = adata.obs[batch]
    adata["cell_type"] = adata.obs[cell_type]
    adata["condition"] = adata.obs[condition]

    # 3) Distance to cell-type mean (three metrics)
    print("Computing distance to cell-type mean…")
    dcb.distance_to_celltype_mean(adata,batch=batch)

    # 3.1)
    print("Computing distance to cell-type mean with invariant genes…")
    dcb.distance_to_celltype_mean_invariant(adata,batch=batch)

    # 4) Scallop pipeline: membership stability → noise = 1 - mean membership
    print("Running Scallop pipeline for stability noise…")
    dcb.scallop_pipeline(adata, batch=batch) 

    # 5) Optionally, compute GCL (global coordination level)
    print("Computing global coordination level (GCL)…")
    dcb.gcl_per_cell_type_and_batch(adata, batch=batch)

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
    parser.add_argument("--batch", default="batch", help="obs key for batch/donor")
    parser.add_argument("--cell_type", default="cell_type", help="obs key for cell type")
    parser.add_argument("--condition", default="condition", help="obs key for condition")
    args = parser.parse_args()

    run_decibel_noise(args.adata, args.results_dir, args.batch_key)
