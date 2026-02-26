configfile: "config/config.json"

# number of seed
n_seeds = 1000 # TODO: using 100 seeds for testing, change back to 1000 later
import random
random.seed(42)
seeds = random.sample(range(10**4, 10**5 -1), n_seeds) 

# params to test
n_neighbors=[10, 30, 50, 100, 200]
min_dist=[0.0, 0.1, 0.3, 0.5, 0.8, 0.99]
metric=["euclidean", "cosine", "correlation"]
n_epochs=[10000]
n_combinations = len(n_neighbors) * len(min_dist) * len(metric) * len(n_epochs)
rule param_table_hepochs:
    output:
        "outputs/UMAP-high-epochs/param_table.csv"
    params:
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        n_epochs=n_epochs
    run:
        import pandas as pd
        from itertools import product
        param_combinations = list(product(params.n_neighbors, params.min_dist, params.metric, params.n_epochs))
        df = pd.DataFrame(param_combinations, columns=["n_neighbors", "min_dist", "metric", "n_epochs"])
        df.index.name = "param_set"
        df.to_csv(output[0])

rule generate_seeds_hepochs:
    output:
        "outputs/UMAP-high-epochs/seeds.txt"
    params:
        n_seeds=n_seeds
    log:
        "logs/UMAP-high-epochs/generate_seeds.log"
    run:
        import random
        random.seed(42)
        seeds = random.sample(range(10**4, 10**5 -1), params.n_seeds)
        with open(output[0], "w") as f:
            for seed in seeds:
                f.write(f"{seed}\n")


def umap_runtime_hepochs(wildcards):
    multiplier = 5
    runtime_dict = {
        "ath": f"{10 * multiplier}h",
        "ptr": f"{1 * multiplier}h",
        "tma": f"{1.5 * multiplier}h",
        "lch": f"{0.7 * multiplier}h",
        "gar": f"{1 * multiplier}h",
        "osa": f"{2.5 * multiplier}h",
        "zma": f"{1.1 * multiplier}h",
        "zma-5k": f"{1.1 * multiplier}h",
        "zma-1k": f"{1.1 * multiplier}h"
    }
    return runtime_dict.get(wildcards.sample)

def umap_mem_mb_per_cpu(wildcards):
    mem_dict = {
        "ath": 5000,
        "ptr": 858.521,
        "tma": 1011.73,
        "lch": 556.733,
        "gar": 958.709,
        "osa": 5000,
        "zma": 858,
        "zma-5k": 858,
        "zma-1k": 858
    }
    return mem_dict.get(wildcards.sample)


rule run_umap_for_a_seed_hepochs:
    container: "src/umap-learn/umap-learn_0.5.9.post2.sif"
    threads: 15
    resources:
        mem_mb_per_cpu=umap_mem_mb_per_cpu,
        runtime=umap_runtime_hepochs
    input:
        rmd_output="outputs/Seurat/{sample}",
        param_table="outputs/UMAP-high-epochs/param_table.csv"
    output:
        dout=directory("outputs/UMAP-high-epochs/{sample}/seed{seed}")
    log:
        "logs/UMAP-high-epochs/seed{seed}/{sample}.log"
    benchmark:
        "benchmarks/UMAP-high-epochs/seed{seed}/{sample}.txt"
    shell:
        """
        exec python scripts/run_umap_for_a_seed.py \
            --rmd_input {input.rmd_output} \
            --param_table {input.param_table} \
            --seed {wildcards.seed} \
            --output_dir {output.dout} \
            --n_process {threads} \
        > {log} 2>&1
        """


rule call_umap_per_sample_hepochs:
    threads: 1
    resources:
        runtime=30
    input:
        param_table="outputs/UMAP-high-epochs/param_table.csv",
        seed_list="outputs/UMAP-high-epochs/seeds.txt",
        umap_dirs=lambda wildcards: [
            f"outputs/UMAP-high-epochs/{wildcards.sample}/seed{seed}" for seed in seeds
        ]
    output:
        "outputs/UMAP-high-epochs/call_umap_per_sample_per_seed/{sample}/done.txt"
    params:
        clean_dir=False
    log:
        "logs/UMAP-high-epochs/call_umap_per_sample_per_seed/{sample}.log"
    run:
        import os, sys, shutil

        # get number of param sets
        n_param_set = len([line for line in open(input.param_table)]) - 1  # minus header

        # loop through umap dirs
        ALL_SEEDS_OK = True
        with open(output[0], "w") as out_f, open(log[0], "w") as log_f:
            # 
            for umap_dir in input.umap_dirs:
                # check if number of .csv files match number of param sets
                n_files = len([f for f in os.listdir(umap_dir) if f.endswith(".csv")])
                if n_files != n_param_set:
                    print(f"Number of UMAP output files in {umap_dir} ({n_files}) does not match number of parameter sets ({n_param_set})", file=log_f)
                    ALL_SEEDS_OK = False
                    if params.clean_dir:
                        # clean up incomplete output directory
                        shutil.rmtree(umap_dir)
                        print(f"Removed incomplete output directory {umap_dir}", file=log_f)
                else:
                    print(f"UMAP output directory {umap_dir} is complete with {n_files} files.", file=log_f)
            # 
            if ALL_SEEDS_OK:
                print(f"All UMAP outputs for sample {wildcards.sample} are complete.", file=log_f)
                out_f.write("UMAP calls completed successfully.\n")
            else:
                print(f"UMAP outputs for sample {wildcards.sample} are incomplete. Please rerun the workflow.", file=log_f)
                sys.exit(1)





# rule run_umap_test:
#     container: "src/umap-learn/umap-learn_0.5.9.post2"
#     threads: 95
#     resources:
#         mem_mb=600
#     output:
#         dout=directory("outputs/UMAP/umap_test/{repeat}")
#     log:
#         "logs/UMAP/umap_test/{repeat}.log"
#     shell:
#         """
#         python scripts/run_umap_repeat.py \
#             --n_jobs {threads} \
#         2> {log}
#         1> {log}
#         """
#
# rule run_umap_for_all_seeds:
#     container: "src/umap-learn/umap-learn_0.5.9.post2"
#     threads: 5
#     resources:
#         mem_mb=600
#     input:
#         rmd_output="outputs/Seurat/{sample}",
#         param_table="outputs/UMAP/param_table.csv",
#         seed_list="outputs/UMAP/seeds.txt"
#     output:
#         dout=directory("outputs/UMAP/umap/set{param_set}/{sample}")
#     log:
#         "logs/UMAP/umap/set{param_set}/{sample}.log"
#     shell:
#         """
#         python scripts/run_umap_for_all_params.py \
#             --rmd_input {input.rmd_output} \
#             --param_table {input.param_table} \
#             --param_set {wildcards.param_set} \
#             --seed_list {input.seed_list} \
#             --output_dir {output.dout} \
#             --n_process {threads} \
#         2> {log} 
#         1> {log}
#         """
# 
# rule call_umap_per_sample_per_param_set:
#     input:
#         param_table="outputs/UMAP/param_table.csv",
#         umaps=lambda wildcards: [
#             f"outputs/UMAP/umap/set{param_set}/{wildcards.sample}" for param_set in range(n_combinations)
#         ]
#     output:
#         "outputs/UMAP/call_umap_per_sample_per_param_set/{sample}/done.txt"
#     log:
#         "logs/UMAP/call_umap_per_sample_per_param_set/{sample}.log"
#     shell:
#         """
#         echo "UMAP calls completed for sample {wildcards.sample}" > {output[0]}
#         """
