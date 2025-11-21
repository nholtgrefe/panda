# Graphical User-Interface (GUI) for PaNDA

The **PaNDA-GUI** is a Graphical User-Interface for exploring, visualizing and maximizing phylogenetic diversity in phylogenetic networks. To test out the GUI, we recommend using the `eNewick` string of the Xiphophorus network in the file  [`exp2_xiphophorus_network.txt`](https://github.com/nholtgrefe/panda/blob/main/data/exp2_xiphophorus_network.txt).

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/36bebbc1-b5d8-4240-b8bf-31afbfab2820" />

## Installation
The GUI is available for Windows and Linux [here](https://drive.google.com/drive/folders/1QQXV0h1XUMIXVlNTy5vgyK0DtA-QtJTd?usp=drive_link). To install and run, simply download the `.exe` file for your operating system and open/execute it once downloaded. (If preffered, experienced users can also run the GUI from the [source files](https://github.com/nholtgrefe/panda/tree/main/gui/src), although we then suggest using our Python package [`phypanda`](https://github.com/nholtgrefe/panda/tree/main/phypanda) instead.)

## Features

### eNewick input
Displays the `eNewick` string of the network under consideration. The four buttons from left to right:
- `Load from file`: loads a network (in `eNewick`) from a file
- `Paste eNewick`: paste or type the `eNewick` of a network in a textbox
- `Copy eNewick`: copy the displayed `eNewick` string
- `Clear All`: clears the current network

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/67788c32-64d7-4f13-8e2e-b153c7a80e45" />

### Network information
Shows the number of taxa of the network, the total all-paths phylogenetic diversity (i.e. the phylogenetic diversity of the set of all taxa), the level of the network, and the number of reticulation vertices in the network.

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/6f738211-1334-4f77-a1c2-7d457c88538f" />

### Manual taxa selection

Shows a table with the 

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/fcab97c9-6b04-4893-8137-9583472f6bdc" />

### Network viewer

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/f039b2b0-098c-4946-ba8c-8d8f161af0cc" />


### MAPPD Solver

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/03eea5a6-ed9c-481d-a1de-958a833c84e6" />

## Citation
If you use the **PaNDA-GUI**, please cite the corresponding paper:

> **PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks**.
> *Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, Jannik Schestag.*
> bioRxiv, 2025. doi: [10.1101/2025.11.14.688467](https://www.biorxiv.org/content/10.1101/2025.11.14.688467)
