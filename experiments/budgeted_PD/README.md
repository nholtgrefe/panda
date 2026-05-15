# Supplementary data and scripts (budgeted PD paper)

This folder contains the datasets and scripts supporting the experiments in the paper:
> **Tractable Optimization of Budgeted Phylogenetic Diversity on Networks Utilizing Node-Scanwidth**.
> *Niels Holtgrefe and Jannik Schestag*.

The results in this paper (and in this folder) were all obtained using version `2.1.0` of `phypanda`. To install this exact version, run `pip install phypanda==2.1.0`.

All computations were performed on a laptop with Debian 12 on a single core of an Intel® Core™ i5-1245U Processor.

## Dataset

- **Networks** (`networks.csv`): The non-binary networks used for the experiments in the paper were obtained from the binary networks used for the first paper. Specifically, for each network, the 10% shortest edges were contracted, while still adhering to our network condition. This procedure was performed using the script file `01_generate_networks.py`. The file contains 6400 networks with columns: `id` (5-digit network identifier, corresponding to the original network id's), `ntips` (number of leaf taxa, ranging from 20 to 200), `level` (network level), and `newick` (network in eNewick format string).

- **Taxon-costs** (`costs.csv`): The taxon-costs were sampled from a log-normal distribution using the script `02_generate_costs.py`. Each row contains: `id` (network identifier) followed by a space-separated list of taxon–cost pairs in format `(taxon_id,cost)`.

## Experiment

- **(Tree)-Extensions** (`extensions_nodewidth.csv`): The optimal node-scanwidth extensions were computed in the file `03_compute_extensions.py` and saved with columns: `id` (network identifier), `nsw` (node-scanwidth value), `comp_time_xp` (computation time in seconds), and `extension` (a comma-separated ordering of nodes forming the extension). An extension is an alternative representation of a tree-extension.

- **Algorithm benchmark** (`results_benchmark.csv`): The benchmark runs presented in the paper were performed in the script `04_run_benchmark.py`. Each row records timings and solution quality for a single network across three diversity measures (all-paths, max-tree, min-tree) and three budget levels (25%, 50%, 90% of total taxon cost). Columns include: `id`, `level`, `tot_cost` (total taxon cost), `time_net` (network parsing), `time_ext` (extension parsing), and for each measure–budget combination: `time_<measure>_b<budget>` (solve time in seconds) and `val_<measure>_b<budget>` (diversity value achieved). The min-tree measure instead includes `time_mnt_full` and `val_mnt_full` for the full-taxon computation case. An initial warm-up run was performed on a single network to compile the code with numba/JIT; this run was not counted in the reported computation times.

- **Large networks** (`results_large_networks.csv`): The experiments concerning the large 1000-taxon networks were performed in the file `05_run_large_network_experiments.py` and the results (including the networks themselves) were saved in this file.
