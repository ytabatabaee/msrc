import itertools
from collections import Counter, defaultdict

def locus_quartet_support(states, taxa, quartet):
    """
    states: array of length n (one locus)
    quartet: tuple of 4 taxon indices
    Returns:
        supported split as frozenset({frozenset(left), frozenset(right)})
        or None if uninformative
    """
    idx = list(quartet)
    vals = [states[i] for i in idx]

    counts = Counter(vals)
    if sorted(counts.values()) != [2, 2]:
        return None

    groups = defaultdict(list)
    for taxon_idx, val in zip(idx, vals):
        groups[val].append(taxon_idx)

    parts = list(groups.values())
    return frozenset([
        frozenset(parts[0]),
        frozenset(parts[1])
    ])

def extract_quartet_counts(matrix, taxa):
    """
    matrix: L x n rearrangement matrix
    taxa: list of taxa labels
    Returns:
        dict mapping quartet (4 indices) -> Counter of supported splits
    """
    n = len(taxa)
    quartet_counts = {}

    for quartet in itertools.combinations(range(n), 4):
        counts = Counter()

        for locus_states in matrix:
            split = locus_quartet_support(locus_states, taxa, quartet)
            if split is not None:
                counts[split] += 1

        quartet_counts[quartet] = counts

    return quartet_counts

def majority_quartets(quartet_counts):
    """
    Returns:
        dict: quartet -> winning split
    """
    inferred = {}

    for quartet, counts in quartet_counts.items():
        if not counts:
            continue
        inferred[quartet] = counts.most_common(1)[0][0]

    return inferred

def format_quartets(plurality, taxa):
    formatted = []
    for quartet, split in plurality.items():
        (a, b), (c, d) = [list(part) for part in split]
        formatted.append(
            f"({taxa[a]},{taxa[b]})|({taxa[c]},{taxa[d]})"
        )
    return formatted

def write_astral_quartets(
    quartet_counts,
    taxa,
    out_file,
    min_weight=1
):
    """
    quartet_counts: dict from extract_quartet_counts
    taxa: ordered list of taxon labels
    """

    with open(out_file, "w") as f:
        for quartet, counts in quartet_counts.items():
            if not counts:
                continue

            for split, weight in counts.items():
                if weight < min_weight:
                    continue

                parts = [list(part) for part in split]
                a, b = [taxa[i] for i in parts[0]]
                c, d = [taxa[i] for i in parts[1]]

                f.write(f"{a} {b} | {c} {d} {weight}\n")

