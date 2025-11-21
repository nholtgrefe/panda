
import sys, math, random
import networkx as nx
import numpy as np

import pyqtgraph as pg
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QTextEdit, QPushButton, QSpinBox, QComboBox, QListWidget, QFileDialog,
    QCheckBox, QDialog, QDialogButtonBox, QLineEdit, QGroupBox, QGridLayout, 
    QSizePolicy, QStyle, QDoubleSpinBox, QTableWidgetItem, QTableWidget,
    QAbstractItemView, QHeaderView
)
from PyQt5.QtGui import QFont, QColor
from PyQt5.QtCore import Qt, QSize, pyqtSignal

from phypanda import solve_MAPPD, DirectedNetwork, DAG


def phylo_layout(
    G,
    layer_gap=1.5,
    leaf_gap=1.0,
    trials=2000,
    seed=None,
    direction="LR",   # "TD" = top-down, "LR" = left-right
    x_scale=1.5,     
    y_scale=1.0,    
):
    """
    Layout a phylogenetic network (DAG) using a tree-backbone heuristic,
    ensuring monotone direction (root on top or left).

    Parameters
    ----------
    G : nx.DiGraph
        Phylogenetic network (must be a DAG).
    layer_gap : float
        Spacing between hierarchical layers (vertical if TD, horizontal if LR).
    leaf_gap : float
        Spacing between leaves within a layer.
    trials : int
        Number of random child orderings to try for optimization.
    seed : int or None
        Random seed for reproducibility.
    direction : {"TD", "LR"}
        Layout direction:
            "TD" = top-down (root on top, edges go downward)
            "LR" = left-right (root on left, edges go rightward)
    x_scale : float
        Scaling factor applied to the x coordinates.
    y_scale : float
        Scaling factor applied to the y coordinates.
    """
    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("Graph must be a DAG")

    rng = random.Random(seed)

    # --- Step 1: Build a tree backbone ---
    roots = [n for n in G.nodes if G.in_degree(n) == 0]
    if not roots:
        raise ValueError("No root found (no node with in-degree 0).")
    root = roots[0]

    topo_order = list(nx.topological_sort(G))
    topo_index = {n: i for i, n in enumerate(topo_order)}

    tree_edges = []
    for node in topo_order:
        preds = list(G.predecessors(node))
        if preds:
            parent = min(preds, key=lambda p: topo_index[p])
            tree_edges.append((parent, node))

    T = nx.DiGraph(tree_edges)

    # --- Step 2: Recursive layout helper (x only) ---
    def layout_tree(node, depth, x_offset):
        children = list(T.successors(node))
        if not children:
            return {node: (x_offset[0], depth)}, x_offset[0] + leaf_gap

        pos = {}
        child_xs = []
        for ch in children:
            child_pos, x_next = layout_tree(ch, depth + 1, x_offset)
            pos.update(child_pos)
            child_xs.append(pos[ch][0])
            x_offset[0] = x_next

        mean_x = np.mean(child_xs)
        pos[node] = (mean_x, depth)
        return pos, x_offset[0]

    # --- Step 3: Crossing count heuristic ---
    def count_crossings(pos, edges):
        xs = [(pos[u][0], pos[v][0]) for u, v in edges if u in pos and v in pos]
        count = 0
        for i in range(len(xs)):
            x1u, x1v = xs[i]
            for j in range(i + 1, len(xs)):
                x2u, x2v = xs[j]
                if (x1u - x2u) * (x1v - x2v) < 0:
                    count += 1
        return count

    # --- Step 4: Optimization loop (for x-ordering) ---
    best_pos, best_score = None, float("inf")
    node_children = {n: list(T.successors(n)) for n in T.nodes}

    for _ in range(trials):
        for n in node_children:
            rng.shuffle(node_children[n])

        T_tmp = nx.DiGraph()
        for u in node_children:
            for v in node_children[u]:
                T_tmp.add_edge(u, v)

        pos, _ = layout_tree(root, 0, [0.0])
        extra_edges = [e for e in G.edges if e not in tree_edges]
        score = count_crossings(pos, tree_edges + extra_edges)

        if score < best_score:
            best_pos, best_score = pos, score

    pos = best_pos

    # --- Step 5: Assign depth coordinates ---
    depths = {}
    for node in topo_order:
        preds = list(G.predecessors(node))
        if preds:
            depths[node] = max(depths[p] for p in preds) + 1
        else:
            depths[node] = 0
    max_depth = max(depths.values())

    # --- Step 6: Final coordinate mapping ---
    if direction.upper() == "TD":
        # Root at top, children below
        for n in pos:
            x, _ = pos[n]
            y = (max_depth - depths[n]) * layer_gap
            pos[n] = (x * x_scale, y * y_scale)

    elif direction.upper() == "LR":
        # Root on left, children to the right
        for n in pos:
            y, _ = pos[n]  # reuse previous 'x' as vertical coordinate
            x = depths[n] * layer_gap
            pos[n] = (x * x_scale, y * y_scale)

    else:
        raise ValueError("direction must be 'TD' or 'LR'")

    return pos



class RotatedTextItem(pg.TextItem):
    def __init__(self, text, angle=0, **kwargs):
        super().__init__(text, **kwargs)
        self.angle = angle

    def paint(self, painter, *args):
        painter.save()
        # Translate to center before rotation
        b = self.boundingRect()
        painter.translate(b.center())
        painter.rotate(self.angle)
        painter.translate(-b.center())
        super().paint(painter, *args)
        painter.restore()


class NetworkViewer(QWidget):
    leafNodeClicked = pyqtSignal(str)  # emits the node name when a leaf is clicked

    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        layout = QVBoxLayout()

        # --- Graph view ---
        self.graph_view = pg.GraphicsLayoutWidget()
        
        # set background to white
        self.graph_view.setBackground('w')

        #self.bg_image = QImage("/.../panda_face.png")
        #self.bg_brush = pg.QtGui.QBrush(self.bg_image)
        #self.graph_view.setBackground(self.bg_brush)

        self.plot = self.graph_view.addPlot()
        self.plot.setAspectLocked()
        self.plot.hideAxis('left')
        self.plot.hideAxis('bottom')
        layout.addWidget(self.graph_view)

        # --- Checkbox to toggle edge weights ---
        self.weight_checkbox = QCheckBox("Show edge weights")
        self.weight_checkbox.stateChanged.connect(self.toggle_edge_weights)
        layout.addWidget(self.weight_checkbox)

        self.setLayout(layout)

        # --- State tracking ---
        self.edge_items = {}          # (u, v) → (line, arrow)
        self.edge_weight_labels = {}  # (u, v) → TextItem
        self.node_items = {}          # node → Scatter point
        self.label_items = {}         # node → TextItem

    def draw_network(self):
        """Draws the full phylogenetic network."""
        self.plot.clear()
        self.edge_items.clear()
        self.edge_weight_labels.clear()
        self.node_items.clear()
        self.label_items.clear()

        self.graph_view.setBackground('w')

        self.overlay_label = QLabel("PD: -", self.graph_view)
        self.overlay_label.setStyleSheet("color: black; font-weight: bold;")
        self.overlay_label.move(10, 10)  # 10 px from top-left of the widget
        self.overlay_label.setFixedSize(100, 20)  # width=120px, height=40px
        self.overlay_label.show()

        G = self.main_window.network
        pos = phylo_layout(G, trials=100)
        self.node_positions = pos

        # --- Draw edges ---
        for u, v in G.edges():
            x0, y0 = pos[u]
            x1, y1 = pos[v]

            # Draw line
            pen = pg.mkPen('black', width=2)
            line = pg.PlotDataItem([x0, x1], [y0, y1], pen=pen)
            self.plot.addItem(line)

            # Draw arrow
            dx = x1 - x0
            dy = y1 - y0
            angle = 180 - math.degrees(math.atan2(dy, dx))
            arrow = pg.ArrowItem(pos=(x1, y1), angle=angle, brush='black', headLen=12, tipAngle=30)
            self.plot.addItem(arrow)
            self.edge_items[(u, v)] = (line, arrow)

            # --- Edge weight label ---
            weight = G.weight(u, v)
            mid_x, mid_y = (x0 + x1) / 2, (y0 + y1) / 2
            label = pg.TextItem(f"{weight:.2f}", color=(80, 80, 80), anchor=(0.5, 0.5))
            label.setFont(QFont("Arial", 8))
            label.setPos(mid_x, mid_y)
            self.edge_weight_labels[(u, v)] = label
            if self.weight_checkbox.isChecked():
                self.plot.addItem(label)

        # --- Prepare node visuals ---
        xs, ys, brushes, sizes, names = [], [], [], [], []
        for node, (x, y) in pos.items():
            xs.append(x)
            ys.append(y)
            names.append(node)
            if node in G.leaves:
                brushes.append(pg.mkBrush(0, 0, 0))     # black fill
                sizes.append(13)
            else:
                brushes.append(pg.mkBrush(255, 255, 255))  # white fill
                sizes.append(10)

        # --- Draw nodes ---
        self.node_plot = pg.ScatterPlotItem(xs, ys, size=sizes, brush=brushes, data=names)
        self.node_plot.sigClicked.connect(self.on_node_clicked)
        self.plot.addItem(self.node_plot)

        # Track each node’s point for later updates
        for spot in self.node_plot.points():
            node = spot.data()
            self.node_items[node] = spot

        # --- Draw leaf labels ---
        for node, (x, y) in pos.items():
            if node not in G.leaves:
                continue

            # Compute average direction for placement
            edges = list(G.in_edges(node)) + list(G.out_edges(node))
            avg_dx = avg_dy = 0.0
            for u, v in edges:
                if u in pos and v in pos:
                    x0, y0 = pos[u]
                    x1, y1 = pos[v]
                    if node == u:
                        avg_dx += (x1 - x0)
                        avg_dy += (y1 - y0)
                    else:
                        avg_dx += (x0 - x1)
                        avg_dy += (y0 - y1)
            if avg_dx == 0 and avg_dy == 0:
                dx, dy = 0, -0.3
            else:
                norm = (avg_dx**2 + avg_dy**2)**0.5
                dx, dy = -avg_dx / norm * 0.4, -avg_dy / norm * 0.4

            label = RotatedTextItem(node, anchor=(0.0, 0.4), color=(40, 40, 40), angle=0)
            label.setFont(QFont("Arial", 9))#, QFont.Bold))
            label.setPos(x + dx, y + dy)
            label.setZValue(10)
            self.plot.addItem(label)
            self.label_items[node] = label

    # ----------------------------------------------------------
    # --- Utility / UI methods
    # ----------------------------------------------------------

    def toggle_edge_weights(self, state):
        """Show or hide edge weight labels."""
        for label in self.edge_weight_labels.values():
            if state == 2:
                self.plot.addItem(label)
            else:
                self.plot.removeItem(label)

    def highlight_edges(self, highlight_edges):
        """Highlight a given list/set of edges."""
        for (u, v), (line, arrow) in self.edge_items.items():
            if (u, v) in highlight_edges:
                line.setPen(pg.mkPen('red', width=4))
                arrow.setBrush(pg.mkBrush('red'))
            else:
                line.setPen(pg.mkPen('black', width=2))
                arrow.setBrush(pg.mkBrush('black'))

    def on_node_clicked(self, plot, points):
        """Emit signal when a leaf node is clicked."""
        for p in points:
            node = p.data()
            if node in self.main_window.network.leaves:
                self.leafNodeClicked.emit(node)

    # ----------------------------------------------------------
    # --- Highlight multiple leaves
    # ----------------------------------------------------------

    def highlight_leaves(self, target_leaves: list[str]):
        """Highlight multiple leaf nodes and reset all others."""
        target_set = set(target_leaves)

        for node, spot in self.node_items.items():
            if node in self.main_window.network.leaves:
                if node in target_set:
                    spot.setBrush(pg.mkBrush(255, 0, 0))  # bright red
                    spot.setSize(16)
                    if node in self.label_items:
                        label = self.label_items[node]
                        label.setColor((200, 0, 0))
                        f = label.textItem.font()
                        f.setBold(True)
                        label.setFont(f)
                else:
                    spot.setBrush(pg.mkBrush(0, 0, 0))  # default black
                    spot.setSize(13)
                    if node in self.label_items:
                        label = self.label_items[node]
                        label.setColor((40,40,40))
                        f = label.textItem.font()
                        f.setBold(False)
                        label.setFont(f)
            else:
                # internal nodes
                spot.setBrush(pg.mkBrush(255, 255, 255))
                spot.setSize(10)

    def clear_scene(self):
        """Clear all items from the scene."""
        self.plot.clear()
        self.edge_items.clear()
        self.edge_weight_labels.clear()
        self.node_items.clear()
        self.label_items.clear()

        self.graph_view.setBackground(self.bg_brush)

        # remove overlay label if exists
        if hasattr(self, 'overlay_label'):
            self.overlay_label.deleteLater()
            del self.overlay_label





class InputDataFrame(QWidget):
    FIXED_WIDTH = 325  # same fixed width as SelectorFrame

    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.setFixedWidth(self.FIXED_WIDTH)
        self.initUI()

    def initUI(self):
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(4)

        style = self.style()  # use standard system icons

        # ===========================
        # Newick input section
        # ===========================
        newick_box = QGroupBox("eNewick Input")
        newick_box.setStyleSheet("QGroupBox { font-weight: bold; }")
        newick_layout = QVBoxLayout(newick_box)
        newick_layout.setContentsMargins(6, 6, 6, 6)
        newick_layout.setSpacing(4)

        # --- Read-only Newick text box ---
        self.newick_input = QTextEdit()
        self.newick_input.setReadOnly(True)
        self.newick_input.setMaximumHeight(60)
        self.newick_input.setPlaceholderText("Newick string will appear here...")
        newick_layout.addWidget(self.newick_input)

        # --- Buttons row (icon only) ---
        button_layout = QHBoxLayout()

        self.load_button = QPushButton()
        self.load_button.setIcon(style.standardIcon(QStyle.SP_DialogOpenButton))
        self.load_button.setToolTip("Load from File")

        self.paste_button = QPushButton()
        self.paste_button.setIcon(style.standardIcon(QStyle.SP_FileDialogDetailedView))
        self.paste_button.setToolTip("Paste eNewick")

        self.copy_button = QPushButton()
        self.copy_button.setIcon(style.standardIcon(QStyle.SP_FileDialogListView))
        self.copy_button.setToolTip("Copy Newick")

        self.clear_all_button = QPushButton()
        self.clear_all_button.setIcon(style.standardIcon(QStyle.SP_DialogCloseButton))
        self.clear_all_button.setToolTip("Clear All")


        for btn in [self.load_button, self.paste_button, self.copy_button, self.clear_all_button]:
            btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            btn.setIconSize(QSize(20, 20))
            button_layout.addWidget(btn)

        newick_layout.addLayout(button_layout)
        main_layout.addWidget(newick_box)

        # ===========================
        # Metadata section
        # ===========================
        self.meta_group = QGroupBox("Network information")
        self.meta_group.setStyleSheet("QGroupBox { font-weight: bold; }")
        self.meta_layout = QGridLayout(self.meta_group)
        self.meta_layout.setContentsMargins(6, 6, 6, 6)
        self.meta_layout.setSpacing(4)

        def create_value_label():
            label = QLabel("-")
            #label.setStyleSheet("color: gray;")
            return label

        def create_heading(text):
            label = QLabel(text)
            label.setStyleSheet("font-weight: bold;")
            return label

        self.level_label = create_value_label()
        self.reticulations_label = create_value_label()
        self.weight_label = create_value_label()
        self.nrtaxa_label = create_value_label()

        # Left column
        self.meta_layout.addWidget(create_heading("# taxa:"), 0, 0, alignment=Qt.AlignRight)
        self.meta_layout.addWidget(self.nrtaxa_label, 0, 1)
        self.meta_layout.addWidget(create_heading("total PD:"), 1, 0, alignment=Qt.AlignRight)
        self.meta_layout.addWidget(self.weight_label, 1, 1)

        # Right column
        self.meta_layout.addWidget(create_heading("level:"), 0, 2, alignment=Qt.AlignRight)
        self.meta_layout.addWidget(self.level_label, 0, 3)
        self.meta_layout.addWidget(create_heading("# rets:"), 1, 2, alignment=Qt.AlignRight)
        self.meta_layout.addWidget(self.reticulations_label, 1, 3)

        main_layout.addWidget(self.meta_group)

        # ===========================
        # Connect buttons
        # ===========================
        self.load_button.clicked.connect(self.load_newick_file)
        self.paste_button.clicked.connect(self.open_paste_dialog)
        self.copy_button.clicked.connect(self.copy_newick)
        self.clear_all_button.clicked.connect(self.clear_all)


    # ---------------------------
    #  Actions
    # ---------------------------
    def load_newick_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open Newick File", "", "Text Files (*.txt *.nwk)")
        if path:
            with open(path, 'r') as f:
                newick = f.read()
                self.newick_input.setText(newick)
                self.update_network(newick)

    def open_paste_dialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Enter eNewick String")
        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(10, 10, 10, 10)

        input_line = QLineEdit()
        input_line.setPlaceholderText("Paste or type eNewick string and press Enter")
        input_line.setMinimumWidth(400)  # make the popup wider
        layout.addWidget(input_line)

        def on_enter():
            newick = input_line.text().strip()
            if newick:
                self.newick_input.setText(newick)
                self.update_network(newick)
                dialog.accept()

        input_line.returnPressed.connect(on_enter)
        dialog.exec_()

    def copy_newick(self):
        """Copy the Newick string to clipboard."""
        newick_text = self.newick_input.toPlainText()
        if newick_text:
            clipboard = self.main_window.clipboard()
            clipboard.setText(newick_text)

    def update_network(self, newick):
        self.main_window.network = DirectedNetwork()
        self.main_window.network.load_from_enewick(newick)

        # Clear previous data
        self.main_window.panda_frame.output_species.clear()
        self.main_window.selector_frame.selected.clear()
        self.main_window.selector_frame.update_taxa()

        # Draw full network
        self.main_window.network_viewer.draw_network()
        self.update_metadata()

        maximum = len(self.main_window.network.leaves)
        self.main_window.panda_frame.k_input.setMaximum(maximum)
        self.main_window.panda_frame.k_input.setValue(5 if maximum >= 5 else maximum)

    def update_metadata(self):
        network = self.main_window.network
        self.level_label.setText(f"{network.level()}")
        self.reticulations_label.setText(f"{len(network.reticulation_nodes())}")
        self.weight_label.setText(f"{network.total_weight():.2f}")
        self.nrtaxa_label.setText(f"{len(network.leaves)}")

    def clear_all(self):
        """Completely resets the application state."""
        self.newick_input.clear()
        self.main_window.network = None
        self.main_window.panda_frame.clear()
        self.main_window.selector_frame.clear()
        self.main_window.network_viewer.clear_scene()

        for lbl in [self.level_label, self.reticulations_label, self.weight_label, self.nrtaxa_label]:
            lbl.setText("-")
        






class SelectorFrame(QWidget):
    FIXED_WIDTH = 325  # match width of InputDataFrame and other frames

    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window

        # --- selection tracking ---
        self.selected = set()
        self.history = []
        self.history_index = -1
        
        # set fixed width but resizable height
        self.setFixedWidth(self.FIXED_WIDTH)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)

        # initialize UI
        self.initUI()

        self.main_window.network_viewer.leafNodeClicked.connect(self.on_leaf_clicked)


    def initUI(self):
        """Initialize all widgets and layout."""
        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.setSpacing(0)

        # ===========================
        # Main GroupBox
        # ===========================
        group_box = QGroupBox("Manual taxa selection")
        group_box.setStyleSheet("QGroupBox { font-weight: bold; }")
        group_box.setFixedWidth(self.FIXED_WIDTH)
        group_box.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)

        group_layout = QVBoxLayout()
        group_layout.setContentsMargins(6, 6, 6, 6)
        group_layout.setSpacing(4)
        group_box.setLayout(group_layout)

        outer_layout.addWidget(group_box)

        # ===========================
        # Top row: icon buttons
        # ===========================
        button_layout = QHBoxLayout()
        style = self.style()  # use standard system icons

        self.prev_button = QPushButton()
        self.prev_button.setIcon(style.standardIcon(QStyle.SP_ArrowLeft))
        self.prev_button.setToolTip("Go to previous selection")

        self.next_button = QPushButton()
        self.next_button.setIcon(style.standardIcon(QStyle.SP_ArrowRight))
        self.next_button.setToolTip("Go to next selection")

        self.select_all_button = QPushButton()
        self.select_all_button.setIcon(style.standardIcon(QStyle.SP_DialogYesButton))
        self.select_all_button.setToolTip("Select all taxa")

        self.deselect_all_button = QPushButton()
        self.deselect_all_button.setIcon(style.standardIcon(QStyle.SP_DialogNoButton))
        self.deselect_all_button.setToolTip("Deselect all taxa")

        for btn in [self.select_all_button, self.deselect_all_button, self.prev_button, self.next_button]:
            btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            btn.setIconSize(QSize(20, 20))
            button_layout.addWidget(btn)

        group_layout.addLayout(button_layout)

        # --- connect button actions ---
        self.prev_button.clicked.connect(self.go_previous)
        self.next_button.clicked.connect(self.go_next)
        self.select_all_button.clicked.connect(self.select_all)
        self.deselect_all_button.clicked.connect(self.deselect_all)

        # ===========================
        # TABLE
        # ===========================

        self.taxa_table = QTableWidget()
        self.taxa_table.setColumnCount(3)
        self.taxa_table.setHorizontalHeaderLabels(["Taxon", "mPD", "iPD"])
        self.taxa_table.horizontalHeaderItem(1).setToolTip("Marginal Phylogenetic Diversity: increase/decrease in PD if this taxon is added/removed to/from the current selection")
        self.taxa_table.horizontalHeaderItem(2).setToolTip("Individual Phylogenetic Diversity: PD when only this taxon is saved")

        # Selection behavior
        self.taxa_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.taxa_table.setSelectionMode(QAbstractItemView.MultiSelection)
        self.taxa_table.setSortingEnabled(True)

        # Visual tweaks to make it list-like
        self.taxa_table.setShowGrid(False)
        self.taxa_table.verticalHeader().setVisible(False)
        self.taxa_table.setAlternatingRowColors(False)
        self.taxa_table.setStyleSheet("""
            QTableWidget {
                border: none;
                font-size: 9pt;
            }
            QTableWidget::item {
                padding: 1px;
            }
        """)

        # Column widths
        self.taxa_table.horizontalHeader().setStretchLastSection(False)
        self.taxa_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.taxa_table.setColumnWidth(1, 55)
        self.taxa_table.setColumnWidth(2, 45)

        # Row height
        self.taxa_table.verticalHeader().setDefaultSectionSize(22)

        self.taxa_table.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.taxa_table.setFixedWidth(self.FIXED_WIDTH - 20)

        self.taxa_table.itemSelectionChanged.connect(self.on_selection_changed_table)
        group_layout.addWidget(self.taxa_table)



        
        # ===========================
        # PD display (bottom)
        # ===========================

        # Line 1: PD and # taxa saved
        pd_row = QWidget()
        pd_layout = QHBoxLayout(pd_row)
        pd_layout.setContentsMargins(0, 0, 0, 0)

        self.pd_label = QLabel("<b>PD:</b> -")
        self.taxa_saved_label = QLabel("<b># taxa selected:</b> -")


        pd_layout.addWidget(self.taxa_saved_label)
        pd_layout.addStretch()
        pd_layout.addWidget(self.pd_label)

        group_layout.addWidget(pd_row)


        # initialize button states
        self.update_nav_buttons()

    # ---------------------------
    # Core Methods
    # ---------------------------
    def update_taxa(self):
        """Refresh taxa table from the network."""
        self.taxa_table.setRowCount(0)
        marg_pd_dict = self.main_window.network.marginal_diversities(self.selected)

        for taxon in sorted(self.main_window.network.leaves):
            row = self.taxa_table.rowCount()
            self.taxa_table.insertRow(row)

            value_mpd = marg_pd_dict[taxon]
            value_ipd = self.main_window.network.diversity({taxon})

            # Column 0: Taxon name
            taxon_item = QTableWidgetItem(taxon)
            taxon_item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)

            # Column 1: MPD value (two decimals, numeric sort)
            mpd_item = QTableWidgetItem()
            mpd_item.setData(Qt.EditRole, float(round(value_mpd, 2)))  # numeric value for sorting
            mpd_item.setData(Qt.DisplayRole, f"{value_mpd:.2f}")       # formatted display
            mpd_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            mpd_item.setForeground(QColor("gray"))
            mpd_item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)

            # Column 2: IPD value (two decimals, numeric sort)
            ipd_item = QTableWidgetItem()
            ipd_item.setData(Qt.EditRole, float(round(value_ipd, 2)))  # numeric sort
            ipd_item.setData(Qt.DisplayRole, f"{value_ipd:.2f}")       # display format
            ipd_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            ipd_item.setForeground(QColor("gray"))
            ipd_item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)

            # Insert into table
            self.taxa_table.setItem(row, 0, taxon_item)
            self.taxa_table.setItem(row, 1, mpd_item)
            self.taxa_table.setItem(row, 2, ipd_item)
        


    def on_selection_changed_table(self):
        """Triggered when user changes the selection in the table."""
        current_selection = {
            self.taxa_table.item(row, 0).text()
            for row in range(self.taxa_table.rowCount())
            if self.taxa_table.item(row, 0).isSelected()
        }
        if current_selection != self.selected and self.main_window.network is not None:
            self.selected = current_selection
            self.record_selection()
            self.update_subnetwork()
            self.fill_marginal_pd()


    def fill_marginal_pd(self):
        """Compute marginal PD for each taxon in the current selection."""
        marg_pd_dict = self.main_window.network.marginal_diversities(self.selected)

        for row in range(self.taxa_table.rowCount()):
            taxon = self.taxa_table.item(row, 0).text()
            marg_pd = marg_pd_dict[taxon]

            item = self.taxa_table.item(row, 1)
            # Update numeric value for sorting
            item.setData(Qt.EditRole, float(round(marg_pd, 2)))
            # Update display text with two decimals
            item.setData(Qt.DisplayRole, f"{marg_pd:.2f}")


    def record_selection(self):
        """Record the current selection in history, trimming future states."""
        if self.history_index >= 0 and self.history[self.history_index] == self.selected:
            return
        if len(self.selected) == 0:
            return
        self.history = self.history[:self.history_index + 1]
        self.history.append(set(self.selected))
        self.history_index += 1
        self.update_nav_buttons()

    def update_nav_buttons(self):
        """Enable/disable navigation buttons based on history state."""
        self.prev_button.setEnabled(self.history_index > 0)
        self.next_button.setEnabled(self.history_index < len(self.history) - 1)

    def go_previous(self):
        """Go back to the previous selection."""
        if self.history_index > 0:
            self.history_index -= 1
            self.apply_selection(self.history[self.history_index])
            self.update_nav_buttons()

    def go_next(self):
        """Go forward to the next selection."""
        if self.history_index < len(self.history) - 1:
            self.history_index += 1
            self.apply_selection(self.history[self.history_index])
            self.update_nav_buttons()

    def apply_selection(self, selection, history_update=False):
        """Apply a recorded selection to the list widget."""
        self.taxa_table.blockSignals(True)
        self.taxa_table.clearSelection()
        for i in range(self.taxa_table.rowCount()):
            item = self.taxa_table.item(i, 0)
            if item.text() in selection:
                self.taxa_table.selectRow(i)
        self.taxa_table.blockSignals(False)

        self.selected = set(selection)
        self.update_subnetwork()
        self.fill_marginal_pd()

        if history_update:
            self.record_selection()

    def on_leaf_clicked(self, node_name):
        new_selection = self.selected ^ {node_name}
        self.apply_selection(new_selection, history_update=True)


    # ---------------------------
    # Selection utilities
    # ---------------------------
    def select_all(self):
        self.taxa_table.selectAll()
        self.on_selection_changed_table()

    def deselect_all(self):
        self.taxa_table.clearSelection()
        self.on_selection_changed_table()

    # ---------------------------
    # Subnetwork + PD
    # ---------------------------
    def update_subnetwork(self):
        """Update network view and PD label."""
        _, highlight_edges = self.main_window.network.unreduced_subnetwork(self.selected)
        self.main_window.network_viewer.highlight_edges(highlight_edges)
        self.main_window.network_viewer.highlight_leaves(self.selected)

        pd_value = self.main_window.network.diversity(self.selected)
        self.pd_label.setText(f"<b>PD</b>: {pd_value:.2f}")
        self.taxa_saved_label.setText(f"<b># taxa selected:</b> {len(self.selected)}")

        self.main_window.network_viewer.overlay_label.setText(f"PD: {pd_value:.2f}")

    def clear(self):
        # clear selection and history
        self.selected.clear()
        self.history.clear()
        self.history_index = -1
        self.update_nav_buttons()

        # clear table
        self.taxa_table.clear()
        self.pd_label.setText("<b>PD:</b> -")
        self.taxa_saved_label.setText("<b># taxa selected:</b> -")

        # Rewrite table titles
        self.taxa_table.setHorizontalHeaderLabels(["Taxon", "mPD", "iPD"])



class SimulatedAnnealingDialog(QDialog):
    def __init__(self, parent=None, initial_settings=None):
        super().__init__(parent)
        self.setWindowTitle("Simulated Annealing Settings")
        self.setModal(True)

        # Use provided settings or defaults
        settings = initial_settings or {
            "max_iter": 100,
            "p_in": 0.9,
            "p_stop": 0.01,
            "init_ext": "cut splitting"
        }

        layout = QVBoxLayout(self)

        # Max iterations
        layout.addWidget(QLabel("Maximum iterations:"))
        self.iterations_spin = QSpinBox()
        self.iterations_spin.setRange(1, 10000)
        self.iterations_spin.setValue(settings["max_iter"])
        layout.addWidget(self.iterations_spin)

        # Initial acceptance probability
        layout.addWidget(QLabel("Initial acceptance probability (p_in):"))
        self.p_in_spin = QDoubleSpinBox()
        self.p_in_spin.setDecimals(3)
        self.p_in_spin.setRange(0.0, 1.0)
        self.p_in_spin.setSingleStep(0.01)
        self.p_in_spin.setValue(settings["p_in"])
        layout.addWidget(self.p_in_spin)

        # Final acceptance probability
        layout.addWidget(QLabel("Final acceptance probability (p_stop):"))
        self.p_stop_spin = QDoubleSpinBox()
        self.p_stop_spin.setDecimals(3)
        self.p_stop_spin.setRange(0.0, 1.0)
        self.p_stop_spin.setSingleStep(0.01)
        self.p_stop_spin.setValue(settings["p_stop"])
        layout.addWidget(self.p_stop_spin)

        # Initial extension method
        layout.addWidget(QLabel("Initial extension method:"))
        self.init_ext_selector = QComboBox()
        self.init_ext_selector.addItems(["cut splitting", "greedy"])
        self.init_ext_selector.setCurrentText(settings["init_ext"])
        layout.addWidget(self.init_ext_selector)

        # OK / Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_settings(self):
        return {
            "max_iter": self.iterations_spin.value(),
            "p_in": self.p_in_spin.value(),
            "p_stop": self.p_stop_spin.value(),
            "init_ext": self.init_ext_selector.currentText()
        }





class PandaFrame(QWidget):
    FIXED_WIDTH = 270

    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.saved_extension = None
        self.saved_scanwidth = None
        
        self.sa_settings = {
            "max_iter": 100,
            "p_in": 0.9,
            "p_stop": 0.01,
            "init_ext": "cut splitting"
        }

        self.init_GUI()

    def clear(self):
        self.saved_extension = None
        self.saved_scanwidth = None
        
        self.sa_settings = {
            "max_iter": 100,
            "p_in": 0.9,
            "p_stop": 0.01,
            "init_ext": "cut splitting"
        }

        # clear the output box
        self.output_species.clear()
        self.pd_label.setText("<b>PD:</b> -")
        self.taxa_saved_label.setText("<b># taxa selected:</b> -")

        # clear the settings box
        self.update_extension_ui()


    def clear_settings_box(self):
        self.saved_extension = None
        self.saved_scanwidth = None
        self.update_extension_ui()

    def init_GUI(self):
        self.setFixedWidth(self.FIXED_WIDTH)
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(6)

        panda_group = QGroupBox()
        panda_group.setTitle("PANDA")
        panda_group.setStyleSheet("QGroupBox { font-weight: bold; }")
        panda_layout = QVBoxLayout(panda_group)
        panda_layout.setContentsMargins(6, 6, 6, 6)
        panda_layout.setSpacing(6)

        panda_layout.addWidget(self.init_settings_box())
        panda_layout.addWidget(self.init_algorithm_box())
        panda_layout.addWidget(self.init_output_box())

        main_layout.addWidget(panda_group)

        self.run_button.clicked.connect(self.run_panda)

    def init_settings_box(self):
        settings_box = QGroupBox("Tree extension (TE)")
        self.settings_layout = QVBoxLayout(settings_box)

        # Tree Extension Method + SA Settings Button
        self.method_row = QWidget()
        method_layout = QHBoxLayout(self.method_row)
        method_layout.setContentsMargins(0, 0, 0, 0)

        self.method_label = QLabel("Method:")
        self.method_selector = QComboBox()
        self.method_selector.addItems([
            "optimal", "simulated annealing", "cut splitting", "greedy"
        ])

        self.sa_settings_button = QPushButton()
        self.sa_settings_button.setIcon(self.style().standardIcon(self.style().SP_FileDialogDetailedView))
        self.sa_settings_button.setToolTip("Simulated Annealing Settings")
        self.sa_settings_button.setFixedSize(24, 24)  # Optional: make it small and square
        self.sa_settings_button.setEnabled(False)
        self.sa_settings_button.clicked.connect(self.open_sa_dialog)


        method_layout.addWidget(self.method_label)
        method_layout.addWidget(self.method_selector)
        method_layout.addWidget(self.sa_settings_button)
        self.settings_layout.addWidget(self.method_row)

        # Extension status area
        self.extension_status_widget = QWidget()
        self.extension_status_layout = QVBoxLayout(self.extension_status_widget)
        self.extension_status_layout.setContentsMargins(0, 0, 0, 0)

        self.no_extension_label = QLabel("No TE saved.")
        self.no_extension_label.setStyleSheet("color: gray;")

        self.use_saved_checkbox = QCheckBox()
        self.use_saved_checkbox.setChecked(True)
        self.use_saved_checkbox.stateChanged.connect(self.update_extension_ui)

        self.extension_status_layout.addWidget(self.no_extension_label)
        self.extension_status_layout.addWidget(self.use_saved_checkbox)

        self.settings_layout.addWidget(self.extension_status_widget)

        self.method_selector.currentTextChanged.connect(self.on_method_changed)

        self.update_extension_ui()  # Initial state

        return settings_box

    def init_algorithm_box(self):
        algo_box = QGroupBox("MAPPD")
        algo_layout = QVBoxLayout(algo_box)

        # Number of species input on one line
        k_row = QWidget()
        k_layout = QHBoxLayout(k_row)
        k_layout.setContentsMargins(0, 0, 0, 0)

        k_label = QLabel("<b># taxa to select</b> (k):")
        self.k_input = QSpinBox()
        self.k_input.setMinimum(1)
        self.k_input.setMaximum(1000)
        self.k_input.setValue(5)
        self.k_input.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)


        k_layout.addWidget(k_label)
        k_layout.addWidget(self.k_input)
        k_layout.addStretch()

        algo_layout.addWidget(k_row)

        # Run button + spinner on same line
        run_row = QWidget()
        run_layout = QHBoxLayout(run_row)
        run_layout.setContentsMargins(0, 0, 0, 0)

        self.run_button = QPushButton("Solve MAPPD")
        self.loading_spinner = QLabel()
        self.loading_spinner.setPixmap(QApplication.style().standardIcon(QApplication.style().SP_BrowserReload).pixmap(16, 16))
        self.loading_spinner.setVisible(False)

        run_layout.addWidget(self.run_button)
        run_layout.addWidget(self.loading_spinner)
        #run_layout.addStretch()

        algo_layout.addWidget(run_row)

        return algo_box

    def init_output_box(self):
        output_box = QGroupBox("Taxa to select")
        output_layout = QVBoxLayout(output_box)

        # Selected species list (read-only)
        self.output_species = QListWidget()
        self.output_species.setSelectionMode(QListWidget.NoSelection)
        output_layout.addWidget(self.output_species)

        # Line 1: PD and # taxa saved (bold, aligned)
        pd_row = QWidget()
        pd_layout = QHBoxLayout(pd_row)
        pd_layout.setContentsMargins(0, 0, 0, 0)

        self.pd_label = QLabel("<b>PD:</b> -")
        self.taxa_saved_label = QLabel("<b># taxa selected</b>: -")

        pd_layout.addWidget(self.taxa_saved_label)
        pd_layout.addStretch()
        pd_layout.addWidget(self.pd_label)

        output_layout.addWidget(pd_row)

        return output_box

    def on_method_changed(self, method):
        if self.saved_extension and self.use_saved_checkbox.isChecked():
            self.sa_settings_button.setEnabled(False)
        else:
            self.sa_settings_button.setEnabled(method == "simulated annealing")

    
    def open_sa_dialog(self):
        dialog = SimulatedAnnealingDialog(self, initial_settings=self.sa_settings)
        if dialog.exec_() == QDialog.Accepted:
            self.sa_settings = dialog.get_settings()


    def update_extension_ui(self):
        has_extension = self.saved_extension is not None

        self.no_extension_label.setVisible(not has_extension)
        self.use_saved_checkbox.setVisible(has_extension)

        if has_extension:
            use_saved = self.use_saved_checkbox.isChecked()
            sw_text = f"Use saved TE (scanwidth = {self.saved_scanwidth})"
            self.use_saved_checkbox.setText(sw_text)

            # Method row enabled/disabled
            self.method_selector.setEnabled(not use_saved)
            self.method_label.setStyleSheet(f"color: {'gray' if use_saved else 'black'};")
            self.method_selector.setStyleSheet(f"color: {'gray' if use_saved else 'black'};")

            # SA button
            self.sa_settings_button.setEnabled(not use_saved and self.method_selector.currentText() == "simulated_annealing")
        else:
            self.method_selector.setEnabled(True)
            self.method_label.setStyleSheet("color: black;")
            self.method_selector.setStyleSheet("color: black;")
            self.sa_settings_button.setEnabled(self.method_selector.currentText() == "simulated_annealing")

    def run_panda(self):
        self.loading_spinner.setVisible(True)
        QApplication.processEvents()

        newick = self.main_window.input_frame.newick_input.toPlainText().strip()
        k = self.k_input.value()
        method = self.method_selector.currentText()

        if not newick:
            self.result_output.setText("Please enter a Newick string.")
            self.loading_spinner.setVisible(False)
            return

        self.main_window.network = DirectedNetwork()
        self.main_window.network.load_from_enewick(newick)
        network = self.main_window.network

        tree_extension = None
        if self.saved_extension and self.use_saved_checkbox.isChecked():
            tree_extension = self.saved_extension
        else:

            tree_extension, scanw = self._compute_tree_extension(method)
            self.saved_extension = tree_extension
            self.saved_scanwidth = scanw
            self.update_extension_ui()

        result = solve_MAPPD(network, k, tree_extension=tree_extension)

        saved_species = result[1]

        self.output_species.clear()
        for species in sorted(saved_species):
            self.output_species.addItem(species)

        self.main_window.input_frame.update_metadata()
        self.main_window.selector_frame.selected = set(result[1])
        #self.main_window.selector_frame.update_taxa()

        self.pd_label.setText(f"<b>PD:</b> {result[0]:.2f}")
        self.taxa_saved_label.setText(f"<b># taxa saved:</b> {len(saved_species)}")
        self.main_window.selector_frame.apply_selection(saved_species, history_update=True)

        self.loading_spinner.setVisible(False)
        self.use_saved_checkbox.setChecked(True)

    def _compute_tree_extension(self, method):
        network = self.main_window.network
        G = DAG(nx.DiGraph(network))

        if method == "optimal":
            res = G.optimal_scanwidth()
        elif method == "cut splitting":
            res = G.cut_splitting_heuristic()
        elif method == "greedy":
            res = G.greedy_heuristic()
        elif method == "simulated annealing":
            alg = "cut" if self.sa_settings["init_ext"] == "cut splitting" else "greedy"
            res = G.simulated_annealing(
                    reduced=True, verbose=False,
                    max_iter=self.sa_settings["max_iter"],
                    p_in=self.sa_settings["p_in"],
                    p_stop=self.sa_settings["p_stop"],
                    init_ext=alg
                )
        else:
            raise ValueError(f"Unknown tree extension method: {method}")
        
        scanwidth = res[0]
        extension = res[1]

        tree_extension = extension.canonical_tree_extension()

        return tree_extension, scanwidth



class PhyPandaMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks")
        self.setGeometry(100, 100, 1600, 900)
        self.network = None

        # Create frames
        self.network_viewer = NetworkViewer(self)
        self.input_frame = InputDataFrame(self)
        self.panda_frame = PandaFrame(self)
        self.selector_frame = SelectorFrame(self)

        # Column 1: Input + Selector stacked vertically
        column1 = QVBoxLayout()
        column1.setContentsMargins(0, 0, 0, 0)  # remove gaps
        column1.setSpacing(0)

        column1.addWidget(self.input_frame, 0)    # fixed height, no stretch
        column1.addWidget(self.selector_frame, 1) # fills remaining vertical space

        # Column 2: Network viewer
        column2 = QVBoxLayout()
        column2.addWidget(self.network_viewer)

        # Column 3: PANDA frame
        column3 = QVBoxLayout()
        column3.addWidget(self.panda_frame)

        # Combine columns into main horizontal layout
        main_layout = QHBoxLayout()
        main_layout.addLayout(column1, stretch=1)
        main_layout.addLayout(column2, stretch=2)
        main_layout.addLayout(column3, stretch=1)

        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PhyPandaMainWindow()

    # Set initial size: width=1200, height=800
    window.resize(1250, 600)

    window.show()

    app.setStyleSheet("""
    QWidget {
        background-color: #f5f5f5;   /* soft, modern background */
        font-family: "Segoe UI", "Arial", sans-serif;
        color: #333;
    }

    QGroupBox {
        border: 1px solid #ccc;
        border-radius: 5px;
        margin-top: 13px;
        background-color: #fafafa;
    }

    QPushButton:hover {
        background-color: #45a049;
    }

    QLineEdit, QTextEdit {
        border: 1px solid #ccc;
        border-radius: 4px;
        padding: 4px;
        background-color: white;
    }

    QHeaderView::section {
        background-color: #e0e0e0;
        padding: 4px;
        border: 1px solid #ccc;
    }
    """)
    
    sys.exit(app.exec_())