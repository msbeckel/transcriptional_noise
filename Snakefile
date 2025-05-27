configfile: "config.yaml"
samples = config['samples']
qc_settings = config['qc']
batches = config.get('batches', {})
cell_types = config.get('cell_types', {})
conditions = config.get('conditions', {})
classes = config.get('classes', {})

rule all:
    input:
        expand(
            [
                "results/{sample}/adata_decibel.h5ad",
                "results/{sample}/decibel_umap.png",
                "results/{sample}/decibel_violin.png",
                "results/{sample}/decibel_splitted_violin.png",
            ],
            sample=samples
        )

rule preprocess:
    input:
        infile = lambda wc: f"data/00_raw/{wc.sample}/adata.h5ad"
    output:
        h5ad = "results/{sample}/processed.h5ad",
        loom = "results/{sample}/processed.loom",
        qc_png = "results/{sample}/qc_violin.png",
        scatter_png = "results/{sample}/qc_scatter_total_vs_genes.png",
        hvg_png = "results/{sample}/qc_filter_hvg.png",
        umap_png = "results/{sample}/umap.png",
    params:
        min_cells = qc_settings['min_cells'],
        hvg = qc_settings['hvg'],
        batch = lambda wc: batches.get(wc.sample, ""),
    conda:
        "envs/scanpy.yaml"
    shell:
        (
            "python scripts/smk_preprocess_scrnaseq.py "
            "{input.infile} {output.h5ad} {output.loom} {output.qc_png} {output.scatter_png} {output.hvg_png} {output.umap_png} "
            "{params.min_cells} {params.hvg} {params.batch}"
        )

rule decibel:
    input:
        infile = lambda wc: f"results/{wc.sample}/processed.h5ad"
    output:
        h5ad = "results/{sample}/adata_decibel.h5ad",
        umap_png = "results/{sample}/decibel_umap.png",
        violin_png = "results/{sample}/decibel_violin.png",
        splitted_violin_png = "results/{sample}/decibel_splitted_violin.png",
    params:
        batch     = lambda wc: batches.get(wc.sample, ""),
        cell_type = lambda wc: cell_types.get(wc.sample, ""),
        condition = lambda wc: conditions.get(wc.sample, ""),
        classes   = lambda wc: classes.get(wc.sample, ""),
    conda:
        "envs/scanpy.yaml"
    shell:
        (
            "python scripts/smk_decibel_noise.py "
            "{input.infile} {output.h5ad} {output.umap_png} {output.violin_png} {output.splitted_violin_png} "
            "{params.cell_type} {params.condition} {params.batch} {params.classes}"
        )
