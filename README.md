# msrc

**msrc** provides simulation and analysis utilities for studying genome rearrangement characters under a multi-species rearrangement coalescent (MSRC) setting. The repository includes scripts for simulating locus-level rearrangement states from gene trees, building concatenated character matrices, and running distance- or quartet-based downstream analyses.

## Requirements

msrc is implemented in Python 3 and uses the following dependencies:

- [Python 3.x](https://www.python.org)
- [NumPy](https://numpy.org)
- [DendroPy](https://dendropy.org)
- [ETE Toolkit](http://etetoolkit.org)
- [NetworkX](https://networkx.org)
- [SciPy](https://scipy.org)

Some workflows also require:

- [FastME](https://www.atgc-montpellier.fr/fastme/) (available on `PATH`)

If you have Python and pip, install Python dependencies with:

```bash
python -m pip install numpy dendropy ete3 networkx scipy
```

## Usage Instructions
The core simulation script is:

```bash
python /home/runner/work/msrc/msrc/scripts/msrc_simulator.py \
  -g <gene-tree-file> \
  -o <output-file> \
  [-k <kstates>] \
  [-r <rate>]
```

**Arguments**
```
Required
 -g,  --genetrees    input gene trees in Newick format (one tree per line)
 -o,  --output       output rearrangement character matrix path

Optional
 -k,  --kstates      size of rearrangement state space / number of genes [default: 3]
 -r,  --rate         rearrangement rate [default: 1.0]
```

The script reads gene trees, simulates rearrangement evolution per locus, and writes a concatenated matrix (currently exported as `concat.phy` in TNT-style formatting in the default code path).

## Example

### Simulating rearrangement data from gene trees
```bash
python /home/runner/work/msrc/msrc/scripts/msrc_simulator.py \
  -g /path/to/gene_trees.tre \
  -o /path/to/output_matrix.txt \
  -k 20 \
  -r 1.0
```

### Running distance-based analysis

```bash
python /home/runner/work/msrc/msrc/scripts/distance_based.py
```

The `distance_based.py` script includes:
- raw Hamming-style distances on simulated matrices,
- an MSRC-corrected distance transform,
- FastME tree inference from distance matrices.

## Repository Contents

- `/home/runner/work/msrc/msrc/scripts/msrc_simulator.py`  
  Main simulator for rearrangement data generation from Newick gene trees.
- `/home/runner/work/msrc/msrc/scripts/rearrangements.py`  
  Rearrangement primitives (adjacency representation, DCJ operation, binary matrix conversion).
- `/home/runner/work/msrc/msrc/scripts/utils.py`  
  Utilities for parsing taxa and writing matrix/tree outputs.
- `/home/runner/work/msrc/msrc/scripts/distance_based.py`  
  Distance-based inference and MSRC-corrected branch length estimation helpers.
- `/home/runner/work/msrc/msrc/scripts/distance-based.py`  
  Earlier distance-based workflow variant.
- `/home/runner/work/msrc/msrc/scripts/quartet-based.py`  
  Quartet extraction and support-count utilities.
- `/home/runner/work/msrc/msrc/theory/`  
  Theory notes and model derivations (`msrc.pdf`, etc.).

## Notes and Caveats

- Scripts are research-oriented and currently not packaged as an installable Python module.
- Some scripts include hard-coded example filenames in their `__main__` blocks; update paths before running on new datasets.
- Input gene trees are expected in Newick format with one tree per line and taxa labels preserved.
- Model background and notation are documented in `/home/runner/work/msrc/msrc/theory/msrc.pdf`.
