import re
import dendropy
import numpy as np

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text
            for text in re.split(r'(\d+)', s)]

def sort_taxa_naturally(taxa):
    return sorted(taxa, key=natural_sort_key)

def read_taxa_from_gene_trees(gene_tree_file):
    """
    Reads all taxa appearing in a SimPhy-style gene tree file
    (one Newick tree per line).
    Returns a sorted list of taxon labels.
    """
    taxa = set()

    with open(gene_tree_file) as f:
        for line in f:
            newick = line.strip()
            if not newick:
                continue

            tree = dendropy.Tree.get(
                data=newick,
                schema="newick",
                preserve_underscores=True
            )

            for leaf in tree.leaf_node_iter():
                taxa.add(leaf.taxon.label)

    return sort_taxa_naturally(taxa)

def create_matrix(matrix, taxa, out_file):
    header = "\t".join(taxa)
    np.savetxt(
        out_file,
        matrix,
        fmt="%d",
        delimiter="\t",
        header=header,
        comments=""
    )

def write_phylip_concatenated(matrix, taxa, out_file):
    """
    Writes a concatenated character matrix in PHYLIP format.
    """
    L, n = matrix.shape
    with open(out_file, "w") as f:
        f.write(f"{n} {L}\n")
        for j, taxon in enumerate(taxa):
            states = "".join(str(matrix[i,j]) for i in range(L))
            f.write(f"{taxon} {states}\n")

def make_tnt_taxon_map(taxa):
    """
    taxa: list of original taxon labels (strings or ints)
    returns:
        forward: original -> TNT-safe
        reverse: TNT-safe -> original
    """
    forward = {}
    reverse = {}

    for t in taxa:
        t_str = str(t)
        safe = f"t{t_str}"
        forward[t_str] = safe
        reverse[safe] = t_str

    return forward, reverse

def write_tnt_matrix_with_map(matrix, taxa, out_file, taxon_map):
    """
        matrix: n x M binary numpy array
        taxa: list of taxon labels
        taxon_map: map of integer taxon labels to string labels
    """
    n, M = matrix.shape
    print(n, M)
    print(len(taxa))

    with open(out_file, "w") as f:
        f.write("xread\n")
        f.write("'MSRC concatenated rearrangement matrix'\n")
        f.write(f"{n} {M}\n")

        for i, taxon in enumerate(taxa):
            label = taxon_map[str(taxon)][:10]  # TNT-safe
            seq = "".join(str(int(x)) for x in matrix[i])
            f.write(f"{label} {seq}\n")

        f.write(";\n")

def relabel_newick_tree(tree_file, reverse_map, out_file):
    with open(tree_file) as f:
        newick = f.read()

    for safe, orig in reverse_map.items():
        newick = newick.replace(safe, orig)

    with open(out_file, "w") as f:
        f.write(newick)




