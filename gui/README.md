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

Shows a table with all taxa of the network, allowing the user to manually select taxa. The user can also use the buttons `Select All` and `Deselect All` to select or deselect all taxa at once. The other two buttons allow the user to go to the previous or next selection. The total PD of the selected taxa is displayed below the table.

The second column of the table shows the marginal PD-improvement of each taxon, i.e., the increase in phylogenetic diversity when adding that taxon to the current selection (if the taxon is not selected), or the decrease in phylogenetic diversity when removing that taxon from the current selection (if the taxon is selected). The third column shows the individual PD of each taxon, i.e., phylogenetic diversity of the taxon alone. The table can also be sorted by clicking on the column headers.

The selected taxa will also be highlighted in the network viewer. Similarly, selecting taxa in the network viewer will update the selection in the table.

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/fcab97c9-6b04-4893-8137-9583472f6bdc" />

### Network viewer

Shows the phylogenetic network. The user can zoom in and out using the mouse wheel, and can pan the network by clicking and dragging. The user can also select taxa by clicking on them (this will also update the selection in the table on the left). Selected taxa are highlighted in thick red, with the paths to the root (i.e., the ones contributing to the PD-score) also highlighted in red. The total PD of the selected taxa is in the top left of the network viewer. Edge weights can be toggled on and off using the checkbox in the bottom left.

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/f039b2b0-098c-4946-ba8c-8d8f161af0cc" />


### MAPPD Solver

The middle box allows the user to run the MAPPD solver to maximize phylogenetic diversity for a chosen value of _k_ (number of taxa to select). The bottom box displays the optimal set of _k_ taxa to be selected, with the corresponding PD-score. The network viewer and manual selection table are also updated accordingly. 

The box on the top allows the user to choose how to compute the tree-extension used in the MAPPD algorithm, using one of the four algorithms from the paper 'Exact and Heuristic Computation of the Scanwidth of Directed Acyclic Graphs' by 
Niels Holtgrefe, Leo van Iersel & Mark Jones. The algorithm either uses an optimal algorithm to find a tree-extension, a 'cut-splitting' heuristic, a greedy heuristic, or simulated annealing (which allows further settings to be adjusted in the box next to it). The box also displays the scanwidth of the tree-extension used when the MAPPD solver is run. After a first run of the MAPPD solver, the program will reuses the same tree-extension unless the user toggles the corresponding checkbox.

<img width="1586" height="630" alt="Screenshot from 2025-11-21 20-28-13" src="https://github.com/user-attachments/assets/03eea5a6-ed9c-481d-a1de-958a833c84e6" />

## Citation
If you use the **PaNDA-GUI**, please cite the corresponding paper:

> **PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks**.
> *Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, Jannik Schestag.*
> bioRxiv, 2025. doi: [10.1101/2025.11.14.688467](https://www.biorxiv.org/content/10.1101/2025.11.14.688467)
