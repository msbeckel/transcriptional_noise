import numpy as np
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import sys, os
from scipy.stats import median_abs_deviation

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

# Arguments
# argv[1]: input h5ad
# argv[2]: output processed h5ad
# argv[3]: output loom
# argv[4]: qc violin plot path
# argv[5]: umap plot path
# argv[6]: optional batch label (pass empty string if none)
(
    infile, out_h5ad, out_loom, qc_png, scatter_png, hvg_png, umap_png,
    min_cells, hvg, batch_label
) = sys.argv[1:]

# Convert types
hvg = int(hvg)
min_cells = int(min_cells)
batch_label = batch_label if batch_label else None

# Ensure output directory exists
os.makedirs(os.path.dirname(out_h5ad), exist_ok=True)

# Load data
adata = sc.read_h5ad(infile)

# Annotate batch if provided
if batch_label:
    adata.obs['batch'] = adata.obs[batch_label]

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
    adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20], inplace=True, log1p=True
)

# Violin plots for QC
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    multi_panel=True,
    show=False
)
# manually save to provided path
os.makedirs(os.path.dirname(qc_png), exist_ok=True)
plt.savefig(qc_png, bbox_inches='tight')
plt.close()


# Filter cells and genes based on config thresholds
print(f"📉 Filtering cells by counts…")
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)

# Filter cells and genes based on config thresholds
print(f"📉 Filtering cells by MT counts…")
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)

print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

# Violin plots for QC
sc.pl.scatter(
    adata,
    "total_counts", 
    "n_genes_by_counts", 
    color="pct_counts_mt",
    show=False
)
# manually save to provided path
os.makedirs(os.path.dirname(scatter_png), exist_ok=True)
plt.savefig(scatter_png, bbox_inches='tight')
plt.close()

# Doublet scoring using built-in Scanpy
print("🔎 Predicting and Filtering Doublets")
sc.pp.scrublet(
    adata,
    expected_doublet_rate=0.06,
    batch_key='batch' if 'batch' in adata.obs else None
)
# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalization and log transform
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

# Identify highly variable genes
print("🔎 Calculating HVG...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=hvg
)

# Violin plots for QC
sc.pl.highly_variable_genes(
    adata,
    show=False
)
# manually save to provided path
os.makedirs(os.path.dirname(hvg_png), exist_ok=True)
plt.savefig(hvg_png, bbox_inches='tight')
plt.close()

ad = adata[:, adata.var.highly_variable]
adata = ad

#Dimensionality reduction
print("🔎 Performing PCA…")
sc.tl.pca(adata)

# Batch correction via Harmony if annotated
if 'batch' in adata.obs:
    sce.pp.harmony_integrate(adata, "batch")
    sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=10)
else:
    sc.pp.neighbors(adata, n_neighbors=10)

# UMAP and clustering
print("🔍 Performing UMAP…")
sc.tl.umap(adata)
print("🔎 Performing tSNE…" )
sc.tl.tsne(adata, n_pcs=20)
print("🔍 Performing Leiden clustering…" )
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2)

# UMAP plot
sc.pl.umap(
    adata,
    color=['leiden', 'doublet_score', 'batch'],
    show=False
)

# manually save to provided path
os.makedirs(os.path.dirname(umap_png), exist_ok=True)
plt.savefig(umap_png, bbox_inches='tight')
plt.close()

# Save processed data and Loom
adata.write_h5ad(out_h5ad)
adata.write_loom(out_loom)