# Graphical User-Interface (GUI) for PaNDA

The **PaNDA-GUI** is a Graphical User-Interface for exploring, visualizing and maximizing phylogenetic diversity in phylogenetic networks. To test out the GUI, we recommend using the `eNewick` string of the Xiphophorus network in the file  [`exp2_xiphophorus_network.txt`](https://github.com/nholtgrefe/panda/blob/main/data/exp2_xiphophorus_network.txt).

<img src="https://github.com/user-attachments/assets/676c980a-cc2e-4f7e-adc3-bf41b155ef8d" alt="Sample Image" width="570" >

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/36bebbc1-b5d8-4240-b8bf-31afbfab2820" />


## Installation
The GUI is available for Windows and Linux [here](https://drive.google.com/drive/folders/1QQXV0h1XUMIXVlNTy5vgyK0DtA-QtJTd?usp=drive_link). To install and run, simply download the `.exe` file for your operating system and open/execute it once downloaded. (If preffered, experienced users can also run the GUI from the [source files](https://github.com/nholtgrefe/panda/tree/main/gui/src), although we then suggest using our Python package [`phypanda`](https://github.com/nholtgrefe/panda/tree/main/phypanda) instead.)

## Features

### File Window
- `Select file`: Allows the user to select a file out of one of the two following options:
	1.  `.fasta` or `.nexus` file containing a multiple sequence alignment;
	2.  `.txt` file containing a dense set of tf-quarnets, where each line contains one tf-quarnet in the folowing format:
		- `SQ: a b c d 1.0` for a quarnet on leaves $\{a,b,c,d\}$ with a split $ab|cd$ and weight 1.0;
		- `4C: a b c d 1.0` for a quarnet on leaves $\{a,b,c,d\}$ with a four-cycle $a,b,c,d$, the leaf $a$ below the reticulation and weight 1.0.
- `File information`: Displays information on the chosen file, such as a list of taxa and the length of the sequence alignment.

<img src="https://github.com/user-attachments/assets/2ab886b6-00c6-482d-9fb4-7d5ccd02b6f1" alt="Sample Image" width="800" >

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/67788c32-64d7-4f13-8e2e-b153c7a80e45" />

### Network Window

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/6f738211-1334-4f77-a1c2-7d457c88538f" />


### Selection

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/fcab97c9-6b04-4893-8137-9583472f6bdc" />

### View

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/f039b2b0-098c-4946-ba8c-8d8f161af0cc" />


### Alg


<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/03eea5a6-ed9c-481d-a1de-958a833c84e6" />

### Algorithm Window
- `Reconstruct networks`: Runs the **Squirrel** algorithm (and possibly the **$\delta$-heuristic**) to reconstruct a phylogenetic level-1 network from the set of tf-quarnets or the sequence alignment.
- `Algorithm settings`: Allows the user to:
	- specify an outgroup (which results in rooted networks)
	- change the $\lambda$-value for the $\delta$-heuristic,
	- change the maximum number of leaves for which the Travelling Salesman Problem is solved optimally,
	- choose whether the quarnet-weights should be used.
- `Algorithm status`: Window to show the progress while the algorithm is running.

<img src="https://github.com/user-attachments/assets/f3f4480e-167c-4df4-9abe-f7147e8d76a7" alt="Sample Image" width="800" >

### Network Window
- `Reconstructed networks`: Table with the reconstructed networks and corresponding (weighted) consistency scores. The table will continuously be updated while the algorithm is running. The table also allows the user to sort the networks according to their scores.
- `eNewick save options`: Allows the user to either save the selected network or all networks into a `.txt` file using the `eNewick` format. In case no outgroup was specified, the semi-directed networks are rooted randomly to obtain an `eNewick` string.
- `Network visualizer`: Basic visualization window to depict the selected network and its corresponding `eNewick` string. The buttons below can be used for navigation, zooming and saving the network as an image. The visualization is not optimized for larger networks, in which case we recommend an external program (such as [Dendroscope 3](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/dendroscope/)) to view the networks in more detail.

<img src="https://github.com/user-attachments/assets/df1775ad-b442-448f-be77-5665ee5634ca" alt="Sample Image" width="800" >

## Citation
If you use the **PaNDA-GUI**, please cite the corresponding paper:

> **PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks**.
> *Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, Jannik Schestag.*
> bioRxiv, 2025. doi: [10.1101/2025.11.14.688467](https://www.biorxiv.org/content/10.1101/2025.11.14.688467)
