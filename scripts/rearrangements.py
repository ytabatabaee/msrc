import random
import numpy as np

def identity_adjacencies(k):
    """
    Circular genome: 1-2-3-...-k-1
    """
    adj = set()
    for i in range(1, k):
        adj.add(frozenset((f"{i}h", f"{i+1}t")))
    adj.add(frozenset((f"{k}h", "1t")))
    return adj

def dcj_operation(adjacencies):
    """
    Apply a random DCJ to a genome represented by adjacencies.
    """
    adjacencies = set(adjacencies)

    a1, a2 = random.sample(list(adjacencies), 2)
    adjacencies.remove(a1)
    adjacencies.remove(a2)

    x1, x2 = tuple(a1)
    y1, y2 = tuple(a2)

    if random.random() < 0.5:
        adjacencies.add(frozenset((x1, y1)))
        adjacencies.add(frozenset((x2, y2)))
    else:
        adjacencies.add(frozenset((x1, y2)))
        adjacencies.add(frozenset((x2, y1)))

    return adjacencies

def dcj_distance(adj1, adj2):
    """
    Compute DCJ distance between two genomes given as adjacency sets.
    Assumes same gene content.
    """
    G = nx.Graph()

    for adj in adj1:
        a, b = adj
        G.add_edge(a, b, color="red")

    for adj in adj2:
        a, b = adj
        G.add_edge(a, b, color="blue")

    cycles = nx.number_connected_components(G)
    total_adjacencies = len(adj1)

    return (total_adjacencies - cycles) // 2


def adjacency_matrix_to_binary_concat(genomes_by_locus, taxa):
    all_adjs = set()
    for locus in genomes_by_locus:
        for adjset in locus.values():
            all_adjs.update(adjset)

    all_adjs = sorted(all_adjs)
    adj_index = {adj: i for i, adj in enumerate(all_adjs)}

    n = len(taxa)
    M = len(all_adjs) * len(genomes_by_locus)
    matrix = np.zeros((n, M), dtype=int)

    col = 0
    for locus in genomes_by_locus:
        for adj in all_adjs:
            for i, taxon in enumerate(taxa):
                if adj in locus[taxon]:
                    matrix[i, col] = 1
            col += 1

    return matrix
