# Supplementary data and scripts
This folder contains the biological datasets, generated networks, simulated sequences, python scripts and results of the simulation experiments of the paper:
> **PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks**.
> *Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, Jannik Schestag*.

The results in the paper and in this repository were all obtained using version `1.0.0` of `phypanda`.

## Simulations
The phylogenetic networks for the simulation experiment in the paper were generated with the `R` script `simulate_networks.R`, using the `SciPhyNetwork` package [2]. The simulated networks are in the file `simulated_networks.txt` (using the `eNewick` format). The Python script that applies PaNDA on these networks is in the file `experiment1.py`. A table containing the results of this experiment is in the file `experiment1_results.txt`.

## Xiphophorus
The `eNewick` representation of the Xiphophorus network in Figure 5 of the paper, originally from [1], is in the file `xiphophorus_network.txt`. The Python script that applies PaNDA to this network is in the file `experiment2.py`.

## References

1. P. Bastide, C. Solís-Lemus, R. Kriebel, K. William Sparks, and C. Ané. Phylogenetic Comparative Methods on Phylogenetic Networks with Reticulations. Systematic Biology, 67(5):800–820, 2018. DOI: [10.1093/sysbio/syy033](https://doi.org/10.1093/sysbio/syy033).
2. J. A. Justison, C. Solís-Lemus, and T. A. Heath. SiPhyNetwork: An R package for simulating phylogenetic networks. Methods in Ecology and Evolution, 14(7):1687–1698, 2023. DOI:[10.1111/2041-210X.14116]( https://doi.org/10.1111/2041-210X.14116).
