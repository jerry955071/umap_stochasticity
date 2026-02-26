# I/O expensive, run one sample at a time (I/O demands increase with decreasing `n`)
rule orthogonal_procrustes_per_n_embeddings:
    conda: "envs/biopy.yaml"
    threads: 40
    resources:
        mem_mb_per_cpu=4000,
        runtime=14400
    input:
        flag="outputs/UMAP/call_umap_per_sample_per_seed/{sample}/done.txt",
        param_table="outputs/UMAP/param_table.csv",
        seed_list="outputs/UMAP/seeds.txt"
    output:
        outdir=directory("outputs/AlignEmbeddings/procrustes_alignment/{sample}/per_{n}_embeddings")        
    params:
        umap_indir="outputs/UMAP/{sample}",
        n=lambda wildcards: wildcards.n
    log:
        "logs/AlignEmbeddings/procrustes_alignment/{sample}/per_{n}_embeddings.log"
    benchmark:
        "benchmarks/AlignEmbeddings/procrustes_alignment/{sample}/per_{n}_embeddings.txt"
    shell:
        """
        # run the script
        python scripts/orthogonal_procrustes_per_n_embeddings.py \
            --umap_indir {params.umap_indir} \
            --n_process {threads} \
            --n_embeddings {params.n} \
            --output_dir {output.outdir} \
            --param_table {input.param_table} \
            --seed_list {input.seed_list} \
        > {log} 2>&1
        """
