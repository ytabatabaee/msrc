from ete3 import Tree
import numpy as np
import random
import dendropy
import re
import argparse

class RearrangementStateSpace:
    def __init__(self, k):
        self.states = list(range(k))  # e.g. {0,1,2}
        self.k = k

    def sample_initial(self):
        return random.choice(self.states)

class RearrangementCTMC:
    def __init__(self, state_space, rate):
        self.state_space = state_space
        self.rate = rate  # λ

    def simulate_branch(self, start_state, branch_length):
        """
        Simulate CTMC along a branch.
        Returns final state.
        """
        t = 0.0
        state = start_state

        while t < branch_length:
            wait = np.random.exponential(1.0 / self.rate)
            t += wait
            if t >= branch_length:
                break
            state = self._jump(state)

        return state

    def _jump(self, current):
        choices = [s for s in self.state_space.states if s != current]
        return random.choice(choices)

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text
            for text in re.split(r'(\d+)', s)]

def sort_taxa_naturally(taxa):
    return sorted(taxa, key=natural_sort_key)

def simulate_rearrangements_on_gene_tree(
    gene_tree_newick,
    state_space,
    rate,
    root_state=None
):
    """
    Simulate rearrangement evolution on a single gene tree.
    Returns dict: taxon -> state
    """

    tree = Tree(gene_tree_newick, format=1)
    ctmc = RearrangementCTMC(state_space, rate)

    # Assign root state
    if root_state is None:
        root_state = state_space.sample_initial()

    tree.add_feature("state", root_state)

    # Traverse tree
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        parent = node.up
        node.state = ctmc.simulate_branch(
            parent.state,
            node.dist
        )

    # Collect tip states
    tip_states = {
        leaf.name: leaf.state
        for leaf in tree.iter_leaves()
    }

    return tip_states

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


def augment_simphy_dataset(
    gene_tree_file,
    k_states=3,
    rearrangement_rate=1.0,
    strict=True
):
    """
    Augments SimPhy gene trees with MSRC rearrangements.

    gene_tree_file: file with one Newick gene tree per line
    Returns:
        matrix: (L x n) numpy array
        taxa: ordered list of taxa
    """

    # Infer taxa automatically
    taxa = read_taxa_from_gene_trees(gene_tree_file)

    state_space = RearrangementStateSpace(k_states)
    dataset = []

    with open(gene_tree_file) as f:
        for locus_id, line in enumerate(f):
            newick = line.strip()
            if not newick:
                continue

            tip_states = simulate_rearrangements_on_gene_tree(
                newick,
                state_space,
                rearrangement_rate
            )

            if strict:
                missing = set(taxa) - set(tip_states)
                if missing:
                    raise ValueError(
                        f"Locus {locus_id}: missing taxa {missing}"
                    )

            locus_vector = [tip_states[t] for t in taxa]
            dataset.append(locus_vector)

    return np.array(dataset), taxa


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MSRC simulator")
    parser.add_argument("-g", "--genetrees", type=str, required=True,
                        help="MSC gene trees in newick format")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output rearrangement character matrix")
    parser.add_argument("-k", "--kstates", type=int, required=False, default=3,
                        help="Size of rearrangement state space")
    parser.add_argument("-r", "--rate", type=float, required=False, default=1.0,
                        help="Rearrangement rate")
    args = parser.parse_args()

    matrix, taxa = augment_simphy_dataset(
        gene_tree_file=args.genetrees,
        k_states=args.kstates,
        rearrangement_rate=args.rate
    )

    create_matrix(matrix, taxa, args.output)

