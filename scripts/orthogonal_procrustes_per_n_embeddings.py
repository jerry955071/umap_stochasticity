from scipy.linalg import orthogonal_procrustes
from pathlib import Path
import concurrent.futures
from tqdm import tqdm
from typing import List
import pandas as pd
import numpy as np
import argparse, sys
import random


# -------------------------
# argument parser
# -------------------------
parser = argparse.ArgumentParser(description="Align embeddings using Procrustes analysis.")
parser.add_argument("--umap_indir", type=str, required=True)
parser.add_argument("--n_process", type=int, default=1)
parser.add_argument("--n_embeddings", type=int, required=True)
parser.add_argument("--output_dir", type=str, required=True)
parser.add_argument("--param_table", type=str, required=True)
parser.add_argument("--seed_list", type=str, required=True)

# -------------------------
# Embedding alignment
# -------------------------
class EmbeddingAlignment:
    @staticmethod
    def procrustes_alignment(Y: np.ndarray, X: np.ndarray):
        Xc = X - X.mean(axis=0)
        Yc = Y - Y.mean(axis=0)
        R, _ = orthogonal_procrustes(Yc, Xc)
        Y_aligned = Yc @ R
        rmsd = np.sqrt(((Xc - Y_aligned) ** 2).sum(axis=1).mean())
        return Y_aligned, rmsd

    @staticmethod
    def align_embeddings(embeddings: List[np.ndarray]):
        n = len(embeddings)
        mean_rmsd = np.zeros(n)
        best_aligned_embeddings = []
        best_mean_rmsd = float("inf")

        for ref_idx in range(n):
            rmsds = []
            aligned_embeddings = []
            for subject_idx in range(n):
                if ref_idx == subject_idx:
                    aligned_embeddings.append(embeddings[subject_idx])
                    continue
                aligned_embedding, rmsd = EmbeddingAlignment.procrustes_alignment(
                    embeddings[subject_idx], embeddings[ref_idx]
                )
                rmsds.append(rmsd)
                aligned_embeddings.append(aligned_embedding)

            mean_rmsd[ref_idx] = np.mean(rmsds)
            if mean_rmsd[ref_idx] < best_mean_rmsd:
                best_mean_rmsd = mean_rmsd[ref_idx]
                best_aligned_embeddings = aligned_embeddings

        ref_idx = np.argmin(mean_rmsd)
        return ref_idx, mean_rmsd, best_aligned_embeddings

# class WorkerFunctions
class WorkerFunctions:
    @staticmethod
    def call_align_embeddings(
        batch_embedding_dfs: List[pd.DataFrame],
        batch_seeds: List[int],
        batch_outdir: Path
    ):
        # load embeddings
        embeddings = [df.values for df in batch_embedding_dfs]

        # align embeddings
        ref_idx, mean_rmsd, aligned_embeddings = EmbeddingAlignment.align_embeddings(embeddings)

        # update DataFrames
        for idx, df, aligned, seed in zip(range(len(batch_embedding_dfs)), batch_embedding_dfs, aligned_embeddings, batch_seeds):
            df[["UMAP1", "UMAP2"]] = aligned
            df["seed"] = str(seed)
            df["reference"] = idx == ref_idx

        out_df = pd.concat(batch_embedding_dfs)

        # write to file
        batch_outdir.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(batch_outdir / "aligned_embeddings.csv")
        mean_rmsd_df = pd.DataFrame({"mean_rmsd": mean_rmsd, "seed": out_df["seed"].unique().tolist()})
        mean_rmsd_df.to_csv(batch_outdir / "mean_rmsd.csv", index=False)
        return

    @staticmethod
    def run_one_param(param_set, umap_indir, seeds, n_embeddings, output_dir):
        path_embeddings = [umap_indir / f"seed{seed}" / f"param_set{param_set}.csv" for seed in seeds]
        # shuffle the list to ensure random batching
        random.seed(42) 
        random.shuffle(path_embeddings)
        embedding_dfs = [pd.read_csv(p, index_col=0) for p in path_embeddings]
        batches = len(seeds) // n_embeddings
        for batch in range(batches):
            batch_embedding_dfs = embedding_dfs[batch * n_embeddings : (batch + 1) * n_embeddings]
            WorkerFunctions.call_align_embeddings(
                batch_embedding_dfs=batch_embedding_dfs,
                batch_seeds=seeds[batch * n_embeddings : (batch + 1) * n_embeddings],
                batch_outdir=output_dir / f"param_set{param_set}/batch_{batch}"
            )

# 
if __name__ == "__main__":
    args = parser.parse_args()
    umap_indir = Path(args.umap_indir)
    n_process = args.n_process
    n_embeddings = args.n_embeddings
    output_dir = Path(args.output_dir)
    param_table = pd.read_csv(args.param_table, index_col=0)
    param_sets = param_table.index.tolist()
    seeds = [i.strip() for i in open(args.seed_list, "r")]

    output_dir.mkdir(parents=True, exist_ok=True)

    # total number of batches to process (approx for writer)
    total_batches = len(param_sets) * (len(seeds) // n_embeddings)

    # run alignment in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_process) as executor:
        futures = [
            executor.submit(
                WorkerFunctions.run_one_param,
                param_set,
                umap_indir,
                seeds,
                n_embeddings,
                output_dir
            )
            for param_set in param_sets
        ]
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
            future.result()

    sys.exit(0)
