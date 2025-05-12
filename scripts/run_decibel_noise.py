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
import matplotlib.pyplot as plt
import matplotlib
import argparse
import scanpy as sc
import scallop as sl
import sys
sys.path.append('/home/mllorens/mbeckel/repos_ext/decibel/module/')
import decibel as dcb
import scipy.sparse as sp

project = "PRJNA795276"

if project == "GSE141044" :
    def change_celltype_names(adata):
        # 2.1) Group clusters by marker genes
        print("Regrouping clusters by marker genes…")
        cluster_to_cell_type = {
            "0": "DG1",        # Cluster 0 → DG
            "7": "DG2",        # Cluster 7 → DG
            "4": "CA3",       # Cluster 4 → CA3
            "2": "Inhibitory",# Cluster 2 → inhibitory
            "1": "CA1",       # Cluster 1 → CA1
            "3": "CA1",       # Cluster 3 → CA1
            "5": "CA1",       # Cluster 5 → CA1
            "6": "CA1"        # Cluster 6 → CA1
        }

        # Ensure Leiden clusters are strings (if stored as integers)
        adata.obs['leiden'] = adata.obs['leiden'].astype(str)

        # Map clusters to regions and create a new column
        adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_cell_type)

elif project == "PRJNA795276":
    def change_celltype_names(adata):
        # 2.1) Group clusters by marker genes
        print("Regrouping clusters by marker genes…")
        cluster_to_cell_type = {
            "0": "Microglia",          # Cluster 0 → Microglia
            "1": "Oligodendro",        # Cluster 1 → Oligodendro
            "2": "Neuroblast",         # Cluster 2 → Neuroblast
            "3": "Astrocyte_qNSC",     # Cluster 3 → Astrocyte_qNSC
            "4": "aNSC_NPC",           # Cluster 4 → aNSC_NPC
        }

        # Ensure Leiden clusters are strings (if stored as integers)
        adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)

        # Map clusters to regions and create a new column
        adata.obs['cell_type'] = adata.obs['cell_type'].map(cluster_to_cell_type)

        # add condition variable
        median_age = adata.obs['Age'].median()
        adata.obs['condition'] = np.where(adata.obs['Age'] < median_age, 'young', 'old')

        #add neighbors 
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

if project == "GSE141044" :
    def noise_plots(project, adata, plots_dir):
        
        plot_dir = project + "_noise_plots"
        tn_plot_dir = os.path.join(plots_dir, plot_dir)
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

        # Stacked violin plots
        markers = ['Prox1', 'Mpped1', 'Mndal', 'Gad1']
        sc.pl.stacked_violin(adata, 
                            markers, 
                            groupby='cell_type', 
                            dendrogram=True, 
                            save="_marker_genes.png", 
                            show=False)

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

        # Violin plots (separate)
        # Extract the required data from adata.obs
        df = sc.get.obs_df(
            adata,
            keys=['scallop_noise', 'cor_dist', 'cell_type', 'age', 'genotype']
        )

        # Convert 'age' and 'genotype' to ordered categoricals
        df['age'] = pd.Categorical(
            df['age'],
            categories=['M6', 'M24'],  # Enforce order: M6 first, M24 second
            ordered=True
        )

        df['genotype'] = pd.Categorical(
            df['genotype'],
            categories=['WT', 'APP'],  # Enforce order: M6 first, M24 second
            ordered=True
        )

        # Create a FacetGrid of violin plots, one per cell_type
        g = sns.FacetGrid(
            df,
            col='cell_type',  # One plot per cell type
            col_wrap=3,       # Adjust based on number of cell types
            height=4,         # Height of each subplot
            aspect=1.5,       # Width-to-height ratio
            sharey=True       # Share y-axis across plots
        )

        # Map violin plots to each subplot
        g.map_dataframe(
            sns.violinplot,
            x='age',               # X-axis: age groups
            y='scallop_noise',     # Y-axis: scallop noise
            hue='genotype',        # Color by genotype
            split=True,            # Split violins by genotype
            inner='quart',         # Show quartiles inside violins
            linewidth=1,           # Border line width
            palette='Set2',        # Color palette
            density_norm='width'   # Scale violins to same width
        )

        # Adjust plot aesthetics
        g.set_titles(col_template="{col_name}")  # Subplot titles = cell_type
        g.set_axis_labels("Age", "Scallop")  # Axis labels
        g.add_legend(title='Genotype')            # Add a shared legend

        # Save the figure
        plt.savefig(os.path.join(tn_plot_dir, "scallop_noise_by_cell_type_age_genotype.png"), bbox_inches='tight')
        plt.show()

elif project == "PRJNA795276":
    def noise_plots(project, adata, plots_dir):
        
        plot_dir = project + "_noise_plots"
        tn_plot_dir = os.path.join(plots_dir, plot_dir)
        os.makedirs(tn_plot_dir, exist_ok=True)
        sc.settings.figdir = tn_plot_dir
        # UMAP plot
        sc.pl.umap(
            adata,
            color=['cell_type', 'condition', 'scallop_noise', 'mean_gcl', 'cor_dist', 'euc_dist'],
            wspace=0.5,
            ncols=3, 
            save="_transcriptional_noise.png",
            show=False
        )

        # Stacked violin plots
        markers = ['Prox1', 'Mpped1', 'Mndal', 'Gad1']
        sc.pl.stacked_violin(adata, 
                            markers, 
                            groupby='cell_type', 
                            dendrogram=True, 
                            save="_marker_genes.png", 
                            show=False)

        # Violin plots (all in one)
        sc.pl.violin(
            adata,
            ['scallop_noise', 'mean_gcl', 'cor_dist', 'euc_dist', 'man_dist'],
            groupby='cell_type',
            jitter=0.4,
            multi_panel=True,
            save="_violin_noise_metrics.png",
            show=False
        )

        # Extract the required data from adata.obs
        df = sc.get.obs_df(
            adata,
            keys=['scallop_noise', 'cell_type', 'age']
        )

        # Convert 'age' to an ordered categorical variable
        df['age'] = pd.Categorical(
            df['age'],
            categories=['M6', 'M24'],  # Enforce order: M6 first, M24 second
            ordered=True
        )

        # Set the theme for the plot
        sns.set_theme(style="whitegrid")

        # Create the violin plot
        plt.figure(figsize=(12, 6))
        sns.violinplot(
            data=df,
            x='cell_type',
            y='scallop_noise',
            hue='age',
            split=True,
            inner='quart',
            linewidth=1,
            palette='Set2'
        )

        # Customize the plot
        plt.xlabel("Cell Type")
        plt.ylabel("Scallop Noise")
        plt.legend(title='Age Group', loc='upper right')
        plt.tight_layout()

        # Save the figure
        plt.savefig(os.path.join(tn_plot_dir, "scallop_noise_by_cell_type_age.png"), bbox_inches='tight')
        plt.show()

def run_decibel_noise(project: str, adata_path: str, results_dir: str, plots_dir: str, batch: str = 'batch', cell_type: str = 'cell_type', condition: str = 'condition', num_divisions: int = 10):
    # 1) Load data
    adata = sc.read_h5ad(adata_path, as_sparse="X")  
    print(f"Loaded AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

    # 2) Preprocess data
    adata.obs["batch"] = adata.obs[batch]
    adata.obs["cell_type"] = adata.obs[cell_type]
    adata.obs["condition"] = adata.obs[condition]

    # 2.1) Regroup clusters by marker genes
    change_celltype_names(adata)
    print(f"Regrouped clusters by marker genes: {adata.obs['cell_type'].unique()}")


    # 3) Distance to cell-type mean (three metrics)
    print("Computing distance to cell-type mean…")
    dcb.distance_to_celltype_mean(adata,batch="batch")

    # 3.1)
    print("Computing distance to cell-type mean with invariant genes…")
    print("Skipping this method…")
    #dcb.distance_to_celltype_mean_invariant(adata,batch="batch")

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
    noise_plots(project, adata, plots_dir)


    # 5.2) Check for NaN values in noise metrics
    print("Na values in noise metrics:")
    print(f"Scallop: {adata.obs['scallop_noise'].isna().sum()}")
    print(f"Mean GCL: {adata.obs['mean_gcl'].isna().sum()}")
    print(f"Correlation distance: {adata.obs['cor_dist'].isna().sum()}")


    # 6) Save results
    print("✅Saving results…")
    os.makedirs(results_dir, exist_ok=True)
    adata_file = project + "_adata_with_decibel_noise.h5ad"
    out_file = os.path.join(results_dir, adata_file)
    adata.write(out_file)
    print(f"Saved annotated AnnData to {out_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Decibel transcriptional noise estimation on an AnnData"
    )
    parser.add_argument("--project", default="", help="Project name")
    parser.add_argument("--adata", required=True, help="Input preprocessed .h5ad file")
    parser.add_argument("--results_dir", default="decibel_results", help="Directory for output .h5ad")
    parser.add_argument("--plots_dir", default="noise_plots", help="Directory for output plots")
    parser.add_argument("--batch", default="batch", help="obs key for batch/donor")
    parser.add_argument("--cell_type", default="cell_type", help="obs key for cell type")
    parser.add_argument("--condition", default="condition", help="obs key for condition")
    parser.add_argument("--num_divisions", default=10, help="num divisions for GCL")
    args = parser.parse_args()

    run_decibel_noise(args.project, args.adata, args.results_dir, args.plots_dir, args.batch, args.cell_type, args.condition, args.num_divisions)
    print("Finished running Decibel noise estimation.")
    print("Results saved to:", args.results_dir)
    print("Done.")
