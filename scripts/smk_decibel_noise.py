import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import scanpy as sc
import scallop as sl
import sys
sys.path.append('/home/mllorens/mbeckel/repos_ext/decibel/module/')
import decibel as dcb
import scipy.sparse as sp

(
    infile, out_h5ad, umap_png, violin_png, splitted_violin_png, cell_type, condition, batch, classes
) = sys.argv[1:]

classes = classes if classes else None

# 1) Load data
# Ensure output directory exists
os.makedirs(os.path.dirname(out_h5ad), exist_ok=True)

# Load data
adata = sc.read_h5ad(infile, as_sparse="X")  
print(f"Loaded AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

# 2) Preprocess data
adata.obs["batch"] = adata.obs[batch]
adata.obs["cell_type"] = adata.obs[cell_type]
adata.obs["condition"] = adata.obs[condition]

# 3) Distance to cell-type mean (three metrics)
print("🔎 Computing distance to cell-type mean…")
dcb.distance_to_celltype_mean(adata,batch="batch")

# 4) Scallop pipeline: membership stability → noise = 1 - mean membership
print("🔍Running Scallop pipeline for stability noise…")
dcb.scallop_pipeline(adata) 

# 5) Optionally, compute GCL (global coordination level)
print("🔎Computing global coordination level (GCL)…")
gcl_results = dcb.gcl_per_cell_type_and_batch(adata, batch="batch", num_divisions=10)

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

# 5.2) Check for NaN values in noise metrics
print("Na values in noise metrics:")
print(f"Scallop: {adata.obs['scallop_noise'].isna().sum()}")
print(f"Mean GCL: {adata.obs['mean_gcl'].isna().sum()}")
print(f"Correlation distance: {adata.obs['cor_dist'].isna().sum()}")

#Plots
print("📊 Ploting…")

### UMAP plot
sc.pl.umap(
    adata,
    color=['cell_type', 'condition', 'batch', 'scallop_noise', 'mean_gcl', 'cor_dist'],
    wspace=0.5,
    ncols=3, 
    show=False
)
# manually save to provided path
os.makedirs(os.path.dirname(umap_png), exist_ok=True)
plt.savefig(umap_png, bbox_inches='tight')
plt.close()

### Violin plots (all in one)
sc.pl.violin(
    adata,
    ['scallop_noise', 'mean_gcl', 'cor_dist'],
    groupby='cell_type',
    jitter=0.4,
    multi_panel=True,
    show=False
)

# manually save to provided path
os.makedirs(os.path.dirname(violin_png), exist_ok=True)
plt.savefig(violin_png, bbox_inches='tight')
plt.close()

if classes:
    adata.obs["classes"] = adata.obs[classes]
    # Violin plots (separate)
    # Extract the required data from adata.obs
    df = sc.get.obs_df(
        adata,
        keys=['scallop_noise', 'cor_dist', 'cell_type', 'condition', 'classes']
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
        x='condition',               # X-axis: age groups
        y='scallop_noise',     # Y-axis: scallop noise
        hue='classes',        # Color by genotype
        split=True,            # Split violins by genotype
        inner='quart',         # Show quartiles inside violins
        linewidth=1,           # Border line width
        palette='Set2',        # Color palette
        density_norm='width'   # Scale violins to same width
    )

    # Adjust plot aesthetics
    g.set_titles(col_template="{col_name}")  # Subplot titles = cell_type
    g.set_axis_labels(condition, "Scallop")  # Axis labels
    g.add_legend(title=classes)            # Add a shared legend

    # manually save to provided path
    os.makedirs(os.path.dirname(splitted_violin_png), exist_ok=True)
    plt.savefig(splitted_violin_png, bbox_inches='tight')
    plt.close()

# Save adata
print("✅ Saving adata with decibel analysis.")
adata.write_h5ad(out_h5ad)