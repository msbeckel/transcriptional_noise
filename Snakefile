configfile: "config.yaml"
samples = config['samples']
qc_settings = config['qc']
batches = config.get('batches', {})

rule all:
    input:
        expand(
            [
                "results/{sample}/processed.h5ad",
                "results/{sample}/processed.loom",
                "results/{sample}/qc_violin.png",
                "results/{sample}/umap.png",
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
        umap_png = "results/{sample}/umap.png",
    params:
        min_genes = qc_settings['min_genes'],
        min_cells = qc_settings['min_cells'],
        max_pct_mt = qc_settings['max_pct_mt'],
        batch = lambda wc: batches.get(wc.sample, ""),
    conda:
        "envs/scanpy.yaml"
    shell:
        (
            "python scripts/smk_preprocess_scrnaseq.py "
            "{input.infile} {output.h5ad} {output.loom} {output.qc_png} {output.umap_png} "
            "{params.min_genes} {params.min_cells} {params.max_pct_mt} {params.batch}"
        )