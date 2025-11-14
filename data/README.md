# Supplementary data and scripts
This folder contains the biological datasets, generated networks, simulated sequences, python scripts and results of the simulation experiments of the paper:
> **PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks**.
> *Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, Jannik Schestag*.

The results in the paper and in this repository were all obtained using version `1.0.0` of `phypanda`.

## Simulations
The phylogenetic networks for the simulation experiment in the PaNDA paper were generated with the `R` script `simulate_networks.R`, using the `SciPhyNetworks` package [2]. The simulated networks are in the file `simulated_networks.txt` (using the `eNewick` format).

## Xiphophorus
The Xiphophorus network in Figure 5 of the paper, originally from [1], is in the file `xiphophorus_network.txt` (in `eNewick` format). The Python script that applies PaNDA on this network is in the file `xiphophorus.py`.

## References

1. P. Bastide, C. Solís-Lemus, R. Kriebel, K. William Sparks, and C. Ané. Phylogenetic Comparative Methods on Phylogenetic Networks with Reticulations. Systematic Biology, 67(5):800–820, 2018. DOI: [10.1093/sysbio/syy033](https://doi.org/10.1093/sysbio/syy033).
