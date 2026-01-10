import numpy as np
import subprocess
from msrc_simulator import *
import networkx as nx
from scipy.optimize import brentq


def msrc_expected_distance(t, rho, Ne):
    """
    Expected rearrangement mismatch under MSRC
    """
    return 0.5 * (1.0 - np.exp(-2 * rho * t) / (1 + 4 * rho * Ne))


def invert_msrc_distance(d_obs, rho, Ne, t_max=10.0):
    """
    Solve E[d(t)] = d_obs for t
    """
    if d_obs <= 0:
        return 0.0
    if d_obs >= 0.5:
        return t_max

    return brentq(
        lambda t: msrc_expected_distance(t, rho, Ne) - d_obs,
        0.0,
        t_max
    )

def robust_msrc_distance_matrix(matrix, eps=1e-8):
    L, n = matrix.shape
    D = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            d = np.sum(matrix[:, i] != matrix[:, j]) / L
            d = min(d, 0.5 - eps)
            D[i, j] = D[j, i] = -np.log(1 - 2 * d)

    return D

def msrc_corrected_distance_matrix(matrix, rho, Ne):
    L, n = matrix.shape
    D_corr = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            d_obs = np.sum(matrix[:, i] != matrix[:, j]) / L
            t_hat = invert_msrc_distance(d_obs, rho, Ne)
            D_corr[i, j] = D_corr[j, i] = t_hat

    return D_corr


def hamming_distance_matrix(matrix):
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

def infer_species_tree_fastme_msrc(
    matrix,
    taxa,
    rho,
    Ne,
    prefix="msrc_corr"
):
    D = msrc_corrected_distance_matrix(matrix, rho, Ne)
    dist_file = f"{prefix}.dist"
    tree_file = f"{prefix}.nwk"

    write_phylip_distance(D, taxa, dist_file)
    run_fastme(dist_file, tree_file)

    return tree_file


def infer_species_tree_fastme(matrix, taxa, prefix="msrc"):
    D = hamming_distance_matrix(matrix)
    dist_file = f"{prefix}.dist"
    tree_file = f"{prefix}.nwk"

    write_phylip_distance(D, taxa, dist_file)
    run_fastme(dist_file, tree_file)

    return tree_file

def dcj_distance_matrix(genomes_by_locus, taxa):
    """
    genomes_by_locus: list of dicts
        genomes_by_locus[locus][taxon] = adjacency set
    """

    n = len(taxa)
    L = len(genomes_by_locus)
    D = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            dist = 0.0
            for locus in genomes_by_locus:
                dist += dcj_distance(
                    locus[taxa[i]],
                    locus[taxa[j]]
                )
            D[i, j] = D[j, i] = dist / L

    return D


if __name__ == "__main__":
    # MSRC-augmented data
    matrix, taxa = augment_simphy_dataset(
        gene_tree_file="truegenetrees_50",
        k_states=5,
        rearrangement_rate=1.0
    )

    # Distance-based inference
    tree_file = infer_species_tree_fastme_msrc(matrix, taxa, rho=1.0, Ne=200000)
    print("FastME tree written to:", tree_file)

