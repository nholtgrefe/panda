# Supplementary data and scripts
This folder contains the datasets and scripts supporting the experiments in the paper:
> **PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks**.
> *Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, Jannik Schestag*.

The results in the paper and in this repository were all obtained using version `1.0.0` of `phypanda`.

## Simulations
The phylogenetic networks for the simulation experiment in the paper were generated with the `R` script `exp1_simulate_networks.R`, using the `SciPhyNetwork` package [1]. The script exhaustively uses the function `sim.bdh.taxa.ssa` to simulate a birth-death-hybridization model until there are 100 $n$-leaf level-$l$ networks for each $n \in \{20, 50, 100, 200\}$ and $l \in \{0, \ldots, 15\}$. The speciation rate $\lambda$ was set at 1.0 and the extinction rate $\mu$ at 0.2, as in the `SciPhyNetwork` documentation. The hybridization rates $\nu$ were sampled uniformly at random from the interval $[0, 2/n]$, since this ensured the level of the networks were roughly in the right range.

The 6400 resulting simulated networks are in the file `exp1_simulated_networks.csv`. The file contains a header line and one line for each network. Each line contains the following information: id-number, number of leaves, level, hybridization rate $\nu$, extinction rate $\mu$, speciation rate $\lambda$, the `eNewick` string of the network. 

The Python script that applies PaNDA to solve MAPPD on these networks is in the file `exp1_solve_MAPPD.py`. The script solves MAPPD on each network for three different values of $k_1 := n/10$, $k_2 := n/2$ and $k_3 := n$, and tracks the computation times.

The results of the experiment are in the file `exp1_results.csv`, containing a header line, and one line per network. Each line contains the following information: id-number (corresponding to the input networks), number of leaves,	scanwidth, level, number of reticulation vertices, time to parse the eNewick string, time to initialize the scanwidth computation, time to compute the scanwidth, time to create a tree-extension, time to solve MAPPD for k1, k2 and k3, computed diversity score for k1, k2, and k3. All computations were performed on laptop with Debian 12 on a single core of an Intel® Core™ i5-1245U Processor.

## Biological data
The `eNewick` representation of the Xiphophorus network in Figure 5 of the paper, originally from [1], is in the file `xiphophorus_network.txt`. The Python script that applies PaNDA to this network is in the file `experiment2.py`.

## References

1. J. A. Justison, C. Solís-Lemus, and T. A. Heath. SiPhyNetwork: An R package for simulating phylogenetic networks. Methods in Ecology and Evolution, 14(7):1687–1698, 2023. DOI:[10.1111/2041-210X.14116]( https://doi.org/10.1111/2041-210X.14116).
2. P. Bastide, C. Solís-Lemus, R. Kriebel, K. William Sparks, and C. Ané. Phylogenetic Comparative Methods on Phylogenetic Networks with Reticulations. Systematic Biology, 67(5):800–820, 2018. DOI: [10.1093/sysbio/syy033](https://doi.org/10.1093/sysbio/syy033).
