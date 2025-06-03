import celltypist
from celltypist import models
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import sys
import os

(
    infile, out_h5ad, model, umap_png
) = sys.argv[1:]

# 1) Load data
# Ensure output directory exists
os.makedirs(os.path.dirname(out_h5ad), exist_ok=True)

# Load data
adata = sc.read_h5ad(infile, as_sparse="X")  
print(f"Loaded AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

# Celt type annotation
predictions = celltypist.annotate(adata, model = model, majority_voting = True)

# Get an `AnnData` with predicted labels embedded into the cell metadata columns.
adata = predictions.to_adata()

# Umap visualization
print("📊 Ploting…")
sc.pl.umap(adata, color = ['predicted_labels', 'majority_voting', 'conf_score'], legend_loc = 'on data')
# manually save to provided path
os.makedirs(os.path.dirname(umap_png), exist_ok=True)
plt.savefig(umap_png, bbox_inches='tight')
plt.close()

# Save adata
print("✅ Saving adata with cell annotation.")
adata.write_h5ad(out_h5ad)