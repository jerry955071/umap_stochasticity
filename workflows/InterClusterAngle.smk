# These rules are I/O heavy, using more threads may not help
rule inter_cluster_angle_per_sample:
    container: "src/umap-learn/umap-learn_0.5.9.post2.sif"
    threads: 10
    resources:
        mem_mb_per_cpu=5000,
        runtime=240
    input:
        flag="outputs/UMAP/call_umap_per_sample_per_seed/{sample}/done.txt",
        cluster_table="outputs/Seurat/{sample}/clusters.csv",
        param_table="outputs/UMAP/param_table.csv",
        seed_list="outputs/UMAP/seeds.txt"
    output:
        directory("outputs/InterClusterAngle/{sample}")
    params:
        umap_dir="outputs/UMAP/{sample}"
    log:
        "logs/InterClusterAngle/inter_cluster_angle_per_sample/{sample}.log"
    benchmark:
        "benchmarks/InterClusterAngle/inter_cluster_angle_per_sample/{sample}.txt"
    shell:
        """
        exec python scripts/inter_cluster_angle_per_sample.py \
            --umap_indir {params.umap_dir} \
            --cluster_table {input.cluster_table} \
            --param_table {input.param_table} \
            --seed_list {input.seed_list} \
            --output_dir {output} \
            --n_process {threads} \
        > {log} 2>&1
        """

rule inter_cluster_angle_aligned_per_sample:
    container: "src/umap-learn/umap-learn_0.5.9.post2.sif"
    threads: 80
    resources:
        mem_mb_per_cpu=5000,
        runtime=480
    input:
        aligned_embeddings="outputs/AlignEmbeddings/procrustes_alignment/{sample}/per_{n}_embeddings",
        cluster_table="outputs/Seurat/{sample}/clusters.csv",
        param_table="outputs/UMAP/param_table.csv",
    output:
        directory("outputs/InterClusterAngle/{sample}/per_{n}_embeddings")
    log:
        "logs/InterClusterAngle/inter_cluster_angle_aligned_per_sample/{sample}/per_{n}_embeddings.log"
    benchmark:
        "benchmarks/InterClusterAngle/inter_cluster_angle_aligned_per_sample/{sample}/per_{n}_embeddings.txt"
    shell:
        """
        exec python scripts/inter_cluster_angle_aligned_per_sample.py \
            --n_process {threads} \
            --umap_indir {input.aligned_embeddings} \
            --cluster_table {input.cluster_table} \
            --param_table {input.param_table} \
            --output_dir {output} \
        > {log} 2>&1
        """
