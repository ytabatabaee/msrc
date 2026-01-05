import numpy as np
import subprocess
from msrc_simulator import *

def rearrangement_distance_matrix(matrix):
    """
    matrix: L x n
    Returns:
        n x n distance matrix
    """
    L, n = matrix.shape
    D = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            diff = np.sum(matrix[:, i] != matrix[:, j])
            D[i, j] = D[j, i] = diff / L

    return D

def write_phylip_distance(D, taxa, out_file):
    n = len(taxa)
    with open(out_file, "w") as f:
        f.write(f"{n}\n")
        for i, taxon in enumerate(taxa):
            row = " ".join(f"{D[i,j]:.6f}" for j in range(n))
            f.write(f"{taxon} {row}\n")


def run_fastme(dist_file, out_tree):
    subprocess.run([
        "fastme",
        "-i", dist_file,
        "-o", out_tree,
        "--nni",  # optional refinement
        "--spr"
    ], check=True)

def infer_species_tree_fastme(matrix, taxa, prefix="msrc"):
    D = rearrangement_distance_matrix(matrix)
    dist_file = f"{prefix}.dist"
    tree_file = f"{prefix}.nwk"

    write_phylip_distance(D, taxa, dist_file)
    run_fastme(dist_file, tree_file)

    return tree_file

if __name__ == "__main__":
    # MSRC-augmented data
    matrix, taxa = augment_simphy_dataset(
        gene_tree_file="simphy_gene_trees.tre",
        k_states=5,
        rearrangement_rate=1.0
    )

    # Distance-based inference
    tree_file = infer_species_tree_fastme(matrix, taxa)
    print("FastME tree written to:", tree_file)

