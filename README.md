# msrc

Multi-species rearrangement coalescent (MSRC) model utilities for simulating genome rearrangement characters from gene trees and running basic downstream tree-inference workflows.

## Repository structure

- `/home/runner/work/msrc/msrc/scripts/msrc_simulator.py`  
  Main simulation script for generating rearrangement data across loci from Newick gene trees.
- `/home/runner/work/msrc/msrc/scripts/rearrangements.py`  
  Rearrangement primitives (adjacency representation, DCJ operation, matrix conversion).
- `/home/runner/work/msrc/msrc/scripts/utils.py`  
  I/O helpers for taxa extraction and matrix export.
- `/home/runner/work/msrc/msrc/scripts/distance_based.py`  
  Distance-based inference utilities including MSRC-corrected distances.
- `/home/runner/work/msrc/msrc/scripts/distance-based.py`  
  Older distance-based variant using Hamming/DCJ-style summaries.
- `/home/runner/work/msrc/msrc/scripts/quartet-based.py`  
  Quartet support extraction utilities.
- `/home/runner/work/msrc/msrc/theory/`  
  PDFs with theoretical notes and model details.

## Requirements

Python 3 with:

- `numpy`
- `dendropy`
- `ete3`
- `networkx`
- `scipy`

Some distance-based workflows also call the external binary:

- `fastme` (must be available on `PATH`)

Example install:

```bash
python -m pip install numpy dendropy ete3 networkx scipy
```

## Input format

Most scripts expect a gene-tree file with:

- one Newick tree per line
- taxon labels preserved (underscores supported)

## Basic usage

Run MSRC simulation on gene trees:

```bash
python /home/runner/work/msrc/msrc/scripts/msrc_simulator.py \
  -g /path/to/gene_trees.tre \
  -o /path/to/output_matrix.txt \
  -k 20 \
  -r 1.0
```

Run distance-based workflow:

```bash
python /home/runner/work/msrc/msrc/scripts/distance_based.py
```

## Notes

- Scripts are research-oriented and currently not packaged as an installable Python module.
- Some script defaults assume local filenames (e.g., example gene-tree files), so update paths/arguments for your dataset.
- Model background and notation are documented in `/home/runner/work/msrc/msrc/theory/msrc.pdf`.
