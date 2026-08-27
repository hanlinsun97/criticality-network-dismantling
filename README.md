# Criticality and universality in network dismantling

This repository contains the code for

- Lorenzo Cirigliano, Claudio Castellano, Minsuk Kim, Filippo Radicchi, Hanlin Sun, Criticality and universality in network dismantling, [Journal citation].
- [Preprint (arXiv)]

- BibTeX entry:

  ```bibtex
  % To be added
  ```

# How to Run the Code

Compile the code with the provided `Makefile`:

```bash
make
```

To generate one Erdős–Rényi network and run one percolation realization:

```text
./percolate --er <N> <avg_k> <graph_seed> <score> <a> <removal_seed> <results.csv> <sequence.txt>
```

Here, `score` is `nbc` or `2CDC`. The bias parameter `a` may be a finite number, `inf`
for maximum-score removal, or `-inf` for minimum-score removal. 

Example:

```bash
./percolate --er 1000 4 1 nbc inf 2 results.csv sequence.txt
```

To use an existing unweighted, undirected simple network:

```text
./percolate <network.txt> <score> <a> <removal_seed> <results.csv> <sequence.txt>
```

The network file must use zero-based node IDs and the format:

```text
N M
u_1 v_1
u_2 v_2
...
u_M v_M
```

`results.csv` contains the columns `removed,p,GC,2C`. Here, `removed` is the
cumulative number of edges already removed at that row, `p` is the
fraction of the edges removed, `GC` is the fraction of nodes in the giant component, and `2C` is the fraction of nodes in the 2-core.
`sequence.txt` contains the complete edge-removal order as
`step edge_id u v`, where `edge_id` identifies the edge in
the original network, and `u v` are its endpoints.

# License

This code is released under the MIT License. See [LICENSE](LICENSE).
