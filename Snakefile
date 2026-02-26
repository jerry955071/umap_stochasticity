include: "workflows/CellRanger.smk"
include: "workflows/Seurat.smk"
include: "workflows/UMAP.smk"
include: "workflows/UMAP-high-epochs.smk"
include: "workflows/InterClusterAngle.smk"
include: "workflows/AlignEmbeddings.smk"

configfile: "config/config.json"

wildcard_constraints:
    not_lch="|".join([s["species"] for s in config["references"] if s["species"] != "lch"])

# Custom functions used by all workflows
from typing import List
def query(d:List[dict], k:str, v:str) -> dict:
    """Return the first dictionary in a list of dictionaries where the value of key k matches v."""
    return [x for x in d if x[k] == v][0]

def query_all(d:List[dict], k:str, v:str, k_out:str) -> List[str]:
    """Return a list of values from key k_out in a list of dictionaries where the value of key k matches v."""
    return [x[k_out] for x in d if x[k] == v]

def _assembly(wildcards):
    """Get genome assembly file path by sample name"""
    species = query(config["samples"], "name", wildcards.sample)["species"]
    return query(config["references"], "species", species)["assembly"]

# samples
samples = ["ath", "gar", "lch", "osa", "ptr", "tma", "zma", "zma-5k", "zma-1k"]

rule Seurat:
    input:
        expand(
            "outputs/Seurat/{sample}",
            sample=samples
        )

rule UMAP:
    input:
        expand(
            "outputs/UMAP/call_umap_per_sample_per_seed/{sample}/done.txt",
            sample=samples
        )

rule UMAP_high_epochs:
    input:
        expand(
            "outputs/UMAP-high-epochs/call_umap_per_sample_per_seed/{sample}/done.txt",
            sample=["ptr", "gar"]
        )

rule InterClusterAngle:
    input:
        expand(
            "outputs/InterClusterAngle/{sample}",
            sample=samples
        )

rule AlignEmbeddings:
    input:
        expand(
            "outputs/InterClusterAngle/{sample}/per_{n}_embeddings",
            sample=samples,
            n=[10, 20, 50, 100]
        )


rule test_slurm:
    input:
        expand("outputs/Snakefile/test_slurm/{repeat}.txt", repeat=range(1,20))

rule test_slurm_unit:
    threads: 90
    resources:
        runtime=1
    output:
        "outputs/Snakefile/test_slurm/{repeat}.txt"
    shell:
        """
        sleep 40
        echo "Test SLURM job $(date)" > {output}
        """