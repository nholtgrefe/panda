import physquirrel as psq
import networkx as nx
import matplotlib.pyplot as plt

from itertools import combinations
from collections import defaultdict
import math

from .scanwidth import DAG, TreeExtension

def powerset(s):
    """Yield subsets of a set s as sets, one at a time."""
    s_list = tuple(s)
    for r in range(len(s_list) + 1):
        for subset in combinations(s_list, r):
            yield set(subset)


class DirectedNetwork(psq.DirectedNetwork):
    """
    Class for directed networks, subclass of nx.DiGraph. The argument 
    'leaves' is optional input for nodes that need to be assigned as leaves.
        self.leaves: set of leaves of the network.
    """
    def __init__(self, leaves=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.weights = {(u,v): 1.0 for (u,v) in self.edges}

    def add_edge(self, u, v, weight=1.0):
        """Adds edge from u to v."""
        super().add_edge(u, v)
        self.weights[(u,v)] = weight
    
    def add_edges_from(self, edges, weights=None):
        """Adds edges."""
        if weights is None or len(weights) != len(edges):
            weights = [1.0 for edge in edges]
        for i, edge in enumerate(edges):
            u = edge[0]
            v = edge[1]
            self.add_edge(u, v, weight = weights[i])
            
    def remove_edge(self, u, v):
        """Removes edge from u to v."""
        super().remove_edge(u, v)
        del self.weights[(u,v)]
        
    def remove_edges_from(self, edges):
        """Removes edges."""
        for edge in edges:
            self.remove_edge(edge[0], edge[1])
    
    def remove_node(self, v):
        """Removes node v from the network."""
        super().remove_node(v)
        if v in self.leaves:
            self.leaves.remove(v)
        self.weights = {k: val for k, val in self.weights.items() if v not in k}
            
    def remove_nodes_from(self, nodes):
        """Removes all nodes in 'nodes' from the network."""
        for v in nodes:
            self.remove_node(v)

    def copy(self):
        """Returns a copy of the network."""
        N = DirectedNetwork()
        N.add_edges_from(self.edges)
        N.add_nodes_from(self.nodes)
        N.add_leaves_from(self.leaves)
        N.weights = self.weights
        return N
    
    def clear(self):
        """Clear the whole network."""
        super().clear()
        self.weights = dict()
    
    def weight(self, u,v):
        """Returns the weight of the edge (u,v)."""
        if (u,v) not in self.edges:
            raise ValueError
        return self.weights[(u,v)]
    
    def total_weight(self):
        """Returns the total weight of all edges in the network."""
        return sum(self.weights.values())

    def suppress_degree2_node(self, v):
        """Suppresses a degree-2 node v (taking weight into account). If this results in parallel edges, their weights are summed."""
        if self.in_degree(v) != 1 or self.out_degree(v) != 1:
            raise ValueError("Node is not degree-2.")
        
        u = list(self.predecessors(v))[0]
        w = list(self.successors(v))[0]

        new_weight = self.weight(u,v) + self.weight(v,w)
        if (u,w) in self.edges:
            # Parallel edge
            self.weights[(u,w)] += new_weight
        else:
            self.add_edge(u, w, weight=new_weight)
        
        self.remove_node(v)
    
    def degree2_nodes(self):
        """Returns a list of all degree-2 nodes."""
        return [v for v in self.nodes if self.in_degree(v) == 1 and self.out_degree(v) == 1]

    def is_binary(self):
        """Checks if the network is binary."""
        for v in self.nodes():
            if self.out_degree(v) > 2:
                return False
            if self.in_degree(v) > 2:
                return False
        return True
    
    def subnetwork(self, leaves):
        """Returns the subnetwork induced by the given leaves. 
        Suppresses resulting degree-2 nodes and parallel edges (taking weight into account)."""
        if not set(leaves).issubset(self.leaves):
            raise ValueError("Given leaves are not all in the network.")
        
        nodes_to_keep = set(leaves)
        for leaf in leaves:
            nodes_to_keep = nodes_to_keep | set(nx.ancestors(self, leaf))
        
        N = self.copy()
        N.remove_nodes_from([v for v in self.nodes if v not in nodes_to_keep])

        # Suppress degree-2 nodes
        while True:
            degree_two_nodes = N.degree2_nodes()
            if len(degree_two_nodes) == 0:
                break
            N.suppress_degree2_node(degree_two_nodes[0])

        return N
    
    def unreduced_subnetwork(self, leaves):
        """Returns the nodes and edges in the subnetwork induced by the given leaves, without suppressing degree-2 nodes."""
        if not set(leaves).issubset(self.leaves):
            raise ValueError("Given leaves are not all in the network.")
        
        nodes_to_keep = set(leaves)
        for leaf in leaves:
            nodes_to_keep = nodes_to_keep | set(nx.ancestors(self, leaf))
        
        edges_to_keep = [(u,v) for (u,v) in self.edges if u in nodes_to_keep and v in nodes_to_keep]

        return nodes_to_keep, edges_to_keep
    
    def diversity(self, leaves):
        """Returns the total weight of the subnetwork induced by the given leaves."""
        _, sub_edges = self.unreduced_subnetwork(leaves)
        return sum(self.weights[(u,v)] for (u,v) in sub_edges)
    
    def marginal_diversities(self, saved_leaves):
        """Returns a dictionary mapping each leaf to its marginal diversity contribution, given some set
        of leaves that are saved. For a leaf in saved_leaves, the marginal diversity is defined as its decrease when
        removing it from the set."""
        total_div = self.diversity(saved_leaves)
        marg_div = {}
        for leaf in self.leaves:
            if leaf in saved_leaves:
                div_minus = self.diversity(set(saved_leaves) - {leaf})
                marg_div[leaf] = div_minus - total_div
            else:
                div_plus = self.diversity(set(saved_leaves) | {leaf})
                marg_div[leaf] = div_plus - total_div
        return marg_div

    def greedy_max_pd(self, k):
        """Computes a greedy solution for MAX-PD with parameter k."""
        if k < 0 or k > len(self.leaves):
            raise ValueError("Invalid value for k.")
        
        saved_leaves = set()
        for _ in range(k):
            marg_div = self.marginal_diversities(saved_leaves)
            best_leaf = max(marg_div, key=marg_div.get)
            saved_leaves.add(best_leaf)
        
        total_div = self.diversity(saved_leaves)
        return total_div, sorted(list(saved_leaves))

    def is_ultrametric(self, tol=1e-6):
        """Checks if the network is ultrametric, i.e. if all leaves have the same distance from the root, within a tolerance tol.
        root = self.root_node()"""
        root = self.root_node()
        dists = []
        for leaf in self.leaves:
            length = nx.shortest_path_length(self, source=root, target=leaf, weight=lambda u,v,d: self.weight(u,v))
            dists.append(length)
        return max(dists) - min(dists) <= tol

    def blobs(self, include_trivial=True, include_leaves=True):
        """Returns a list of node-sets, each of which make up a blob of the network.
        include_trivial/include_leaves indicate whether the trivial/leaf blobs should
        be included."""
        if include_trivial == False: include_leaves = False
        blobs = []
        visited = set()
        
        for component in list(nx.biconnected_components(self.to_undirected())):
            if len(component) > 2:
                blobs.append(component)
                visited = visited | component
                
        for v in self.nodes:
            if v not in visited:
                if v not in self.leaves and include_trivial == True:
                    blobs.append({v})
                elif v in self.leaves and include_leaves == True:
                    blobs.append({v})

        return blobs
    
    def level(self):
        """Returns the level of the network."""
        max_lev = 0
        reticulations = self.reticulation_nodes()
        for blob in self.blobs(include_trivial=False):
            nr_blob_rets = len([v for v in blob if v in reticulations])
            max_lev = max(max_lev, nr_blob_rets)
        return max_lev


    def binary_resolution(self, weights=1.0):
        """Returns a binary resolution of the network. Assigns the given weights to the new edges."""
        high_indegree_nodes = []
        high_outdegree_nodes = []
        
        for v in self.nodes():
            if self.out_degree(v) > 2:
                high_outdegree_nodes.append(v)
            if self.in_degree(v) > 2:
                high_indegree_nodes.append(v)
        
        net = self.copy()
        
        for u in high_outdegree_nodes:
            out_neighbors = list(net.successors(u))[1:] # The first out_neighbor stays
            out_weights = {v:self.weight(u,v) for v in out_neighbors}
            net.remove_edges_from((u, v) for v in out_neighbors)
            
            # Create catterpillar
            last = u
            for i, v in enumerate(out_neighbors):
                if i < len(out_neighbors) - 1:
                    cat_node = psq.utils.id_generator()
                    net.add_node(cat_node)
                    net.add_edge(last, cat_node, weight=weights)
                    net.add_edge(cat_node, v, weight=out_weights[v])
                    last = cat_node
                else:
                    net.add_edge(last, v, weight=out_weights[v])
                    
        
        for v in high_indegree_nodes:
            in_neighbors = list(net.predecessors(v))[1:] # The first in_neighbor stays
            in_weights = {u:self.weight(u,v) for u in in_neighbors}
            net.remove_edges_from((u, v) for i in in_neighbors)
            
            # Create catterpillar
            last = v
            for i, u in enumerate(in_neighbors):
                if i < len(in_neighbors) - 1:
                    cat_node = psq.utils.id_generator()
                    net.add_node(cat_node)
                    net.add_edge(last, cat_node, weight=weights)
                    net.add_edge(u, cat_node, weight=in_weights[u])
                    last = cat_node
                else:
                    net.add_edge(u, last, weight=in_weights[u])

        return net

    def visualize(self, layout='dot', title=None, leaflabels=True, internal_labels=False, font_size=12):
        """Visualization function with several layout-options: ['dot', 'kamada', 'neato', 'twopi', 'circo'].
        If pygraphviz is not installed, use 'kamada'. Optional title can be given to the plot.
        Label printing can be turned off."""
        if layout == 'kamada':
            pos = nx.kamada_kawai_layout(self)
        else:
            pos = nx.drawing.nx_agraph.graphviz_layout(self, prog=layout)
        
        fig, ax = plt.subplots(figsize=(9, 7), dpi=200)

        
        if internal_labels == True:
            int_node_size = 300
        else:
            int_node_size = 50
            
        nx.draw_networkx_nodes(self, pos, nodelist=self.internal_nodes(), node_size=int_node_size, node_color='white',
                           edgecolors='black', alpha=1) #, label="Internal")
        nx.draw_networkx_nodes(self, pos, nodelist=self.leaves, node_size=300, node_color='white',
                           edgecolors='white', alpha=1) #, label="Leaf")
        nx.draw_networkx_edges(self, pos, edgelist=self.edges, edge_color='black', width=1,
                               node_size=200, alpha=0.85, arrows=True, arrowstyle='->', arrowsize=11)    
        if leaflabels == True:
            nx.draw_networkx_labels(self, pos, labels={leaf:leaf for leaf in self.leaves}, font_size=font_size)   
        if internal_labels == True:
            nx.draw_networkx_labels(self, pos, labels={node:node for node in self.internal_nodes()}, font_size=font_size)   


        if title == None: title = "Directed Network"
        plt.title(title)
        plt.show()
        return ax

    def load_from_enewick(self, enewick_string, weights=True):
        """Clears the network and loads the network from the given enewick_string."""
        self.clear()

        import phylox
        from phylox.newick_parser import dinetwork_to_extended_newick, extended_newick_to_dinetwork
        from phylox.constants import LENGTH_ATTR

        net = extended_newick_to_dinetwork(enewick_string, internal_labels=True)
        #mapping = {leaf: net.nodes[leaf].get("label") for leaf in net.leaves if net.nodes[leaf].get("label") is not None}
        mapping = {node: net.nodes[node].get("label") for node in net.nodes if net.nodes[node].get("label") is not None}

        for u, v, data in net.edges(data=True):
            w = data.get(LENGTH_ATTR)
            if w is None: w = 1.0
                
            if u in mapping.keys(): u = mapping[u]
            if v in mapping.keys(): v = mapping[v]
            
            self.add_edge(u,v, weight=w)

        self.add_leaves_from([v for v in self.nodes if self.out_degree(v) == 0])

############################################################

class DP_instance():
    """Class to run the scanwidth FPT-algorithm for MAX-PD. Takes as input
    a binary network, a parameter k and a tree extension.
        self.table := the dynamic programming table
        self.GW := table that stores the sets of edges in each scanwidth_bag GW_v
        """
    
    def __init__(self, network, tree_extension):
        
        self.network = network
        self.tree_extension = tree_extension
        self.minus_infinity = -2 * sum(network.weights.values()) - 1

        self.GW = self._initialize_GW()

        # Nested dictionary [vertex][l][phi]., first value = diversity score, second is pointer for backtracking
        self.table = defaultdict(lambda: defaultdict(dict))

        self.m = self._max_edge_offspring_count()

    def _initialize_GW(self):
        """Computes all GW_v scanwidth bags."""
        res = {v:None for v in self.network.nodes}
        
        for v in self.network.nodes:
            sink = nx.descendants(self.tree_extension.tree, v)
            sink.add(v)
            
            res[v] = set()
            for (u,w) in self.network.edges:
                if u not in sink and w in sink:
                    res[v].add((u,w))
                    
        return res

    def _max_edge_offspring_count(self):
        """Returns the maximum value m_x over all leaves x of the network,
        where m_x = |{e: x \in off(e)}|"""
        
        # indegree map
        indeg = dict(self.network.in_degree())

        # memoization
        edge_count = {}

        def count_edges(u):
            if u in edge_count:
                return edge_count[u]
            total = indeg[u]
            for v in self.network.predecessors(u):  # all incoming neighbors
                total += count_edges(v)
            edge_count[u] = total
            return total

        best_count = 0
        for leaf in self.network.leaves:
            c = count_edges(leaf)
            if c > best_count:
                best_count = c

        return best_count
    
    def _fill_dp_table(self, k):
        """Bottom-up dynamic programming version that fills self.table."""

        # Nested dictionary [vertex][phi][l]. First entry is the value, second is the set of species
        self.table = defaultdict(lambda: defaultdict(dict))

        for v in nx.dfs_postorder_nodes(self.tree_extension.tree):
            if v in self.network.leaves:
                self._process_leaf_node(v, k)
            else:
                self._process_non_leaf_node(v, k)
    
    def _process_leaf_node(self, v, k):
        """Processes a leaf node v and returns the (key, value) pairs to be written to the table, given the parameter k."""

        # scanwidth-bag of v
        GW_v = self.GW.get(v)

        for phi in powerset(GW_v):
            phi_frozen = frozenset(phi)
            for l in range(k + 1):
                if l == 0:
                    dp = 0 if len(phi) == 0 else self.minus_infinity
                    pointer = None
                else:
                    dp = sum(self.network.weight(u, v) for (u, v) in phi)
                    pointer = (0,)

                # Write to table
                self.table[v][l][phi_frozen] = (dp, pointer)
    
    def _process_non_leaf_node(self, v, k):
        """Processes a non-leaf node v and write values to the table, given the parameter k."""

        # scanwidth-bag of v
        GW_v = self.GW.get(v)

        # delta_in and delta_out
        delta_in_v = set(self.network.in_edges(v))
        delta_out_v = set(self.network.out_edges(v))

        # children of v in the tree-extension
        children = list(self.tree_extension.tree.successors(v))
        GW_children = [self.GW.get(child) for child in children]
        
        for phi in powerset(GW_v):
            # Precompute necessary sets in outer loop to limit nr. of set-operations #
            phi_frozen = frozenset(phi)
            phi_len = len(phi)
            
            # Used for heuristic pruning
            bound = math.ceil(phi_len / self.m)

            # val3 equals Omega(phi \cap \delta_in (v))
            val3 = sum(self.network.weight(u, v) for (u, v) in phi & delta_in_v)

            # phi_psi_subsets contains all combinations of phi \cup psi with psi \subseteq \delta_out (v)
            phi_psi_subsets = [phi | psi for psi in powerset(delta_out_v)]

            # Create list containing one list for each child x of v, containing
            # all combinations of (phi \cup psi) \cap GW_x
            phi_psi_GW_subsets = list()
            for i in range(len(children)):
                GW = GW_children[i]
                phi_psi_GW_subsets.append([])
                for phi_psi in phi_psi_subsets:
                    phi_psi_GW_subsets[i].append(frozenset(phi_psi & GW))

            for l in range(k + 1):

                ## Base case 2: l = 0 ##
                if l == 0:
                    dp = 0 if phi_len == 0 else self.minus_infinity
                    pointer = None

                # Heuristic pruning, based on trivial hitting set lower bound, i.e. if l < bound, we can never have that
                # |off(e) \cap A| > 0 for all e \in phi
                elif l  < bound:
                    dp = self.minus_infinity
                    pointer = None

                else:
                    dp = self.minus_infinity - 1
                    pointer = None

                    ## Case: v has one child ##
                    if len(children) == 1:
                        u = children[0]

                        for S in phi_psi_GW_subsets[0]:
                            # DP[u, (phi \cup psi) \cap GW_u, l]
                            val1 = self.table.get(u).get(l).get(S)[0]
                            val = val1 + val3
                            if val >= dp:
                                dp = val
                                pointer = (1, (u, l, S))

                    ## Case: v has two children ##
                    elif len(children) == 2:
                        u, w = children

                        for l_prime in range(l + 1):
                            # Iterate through all combinations of (phi \cup psi) \cap GW_x (with x \in {u, w\})
                            # Note: this loop works since both lists have the same set
                            for i in range(len(phi_psi_GW_subsets[0])):
                                # DP[u, (phi \cup psi) \cap GW_u, l']
                                S1 = phi_psi_GW_subsets[0][i]
                                val1 = self.table.get(u).get(l_prime).get(S1)[0]

                                # DP[w, (phi \cup psi) \cap GW_w, l-l']
                                S2 = phi_psi_GW_subsets[1][i]
                                val2 = self.table.get(w).get(l - l_prime).get(S2)[0]

                                val = val1 + val2 + val3
                                if val >= dp:
                                    dp = val
                                    pointer = (2, (u, l_prime, S1), (w, l-l_prime, S2))

                # Write to table
                self.table[v][l][phi_frozen] = (dp, pointer)
                

    def _backtrack_solution(self, v, l, phi_frozen):
        """Uses backtracking to get the solution set corresponding to the value for
        DP[v][l][phi_frozen]."""
        
        pointer = self.table.get(v).get(l).get(phi_frozen)[1]

        # no saving / invalid options -> return empty set
        if pointer is None:
            return set()
        
        # v is a leaf -> return {v}
        elif pointer[0] == 0:
            return {v}
        
        # v has one child -> return solutions of that child
        elif pointer[0] == 1:
            child, l_child, phi_child = pointer[1]
            return self._backtrack_solution(child, l_child, phi_child)
        
        #v has two children -> return union of solutions
        elif pointer[0] == 2:
            (u, l1, phi_u) = pointer[1]
            (w, l2, phi_w) = pointer[2]
            sol1 = self._backtrack_solution(u, l1, phi_u)
            sol2 = self._backtrack_solution(w, l2, phi_w)
            return sol1 | sol2

    def _solve_MAPPD(self, k, include_solution=True):
        """Solves MAPPD with parameter k."""

        self._fill_dp_table(k)
        
        root = self.network.root_node()
        pd = self.table[root][k][frozenset()][0]
        
        if not include_solution:
            return pd
        
        solution = self._backtrack_solution(root, k, frozenset())

        return pd, sorted(list(solution))
        


def solve_MAPPD(network, k, tree_extension="optimal_XP", include_solution=True, **kwargs):
    """Solves MAPPD on the input network (either an enewick string or network instance), with the parameter k (i.e. saving k species).
    The dynamic programming is done over a tree_extension, either given as part of the input, or
     computed by one of the following methods. (Non-binary networks are first transformed with a binary resolution, hence
    given tree_extensions are then disregarded). The solution set can either be returned or not, at the cost of some additional time
    due to backtracking.
    - "optimal_XP": optimal XP-algorithm
    - "cut_splitting": cut-splitting heuristic
    - "greedy": greedy heuristic
    - "simulated_annealing": simulated annealing heuristic
    Keywords are passed to the tree extension functions. """
    
    if type(network) == str:
        newick = network
        network = DirectedNetwork()
        network.load_from_enewick(newick)

    if k < 0 or k > len(network.leaves):
        raise ValueError("Invalid value for k.")
    
    if not network.is_binary():
        network = network.binary_resolution(weights=0.0)
        
        # If network is non-binary, a tree-extension of the original network can not be used.
        if isinstance(tree_extension, TreeExtension):
            tree_extension = "optimal_XP"

    if not isinstance(tree_extension, TreeExtension):
        G = DAG(nx.DiGraph(network))
        
        if tree_extension == "optimal_XP":
            res = G.optimal_scanwidth(**kwargs)
        elif tree_extension == "cut_splitting":
            res = G.cut_splitting_heuristic(**kwargs)
        elif tree_extension == "greedy":
            res = G.greedy_heuristic(**kwargs)
        elif tree_extension == "simulated_annealing":
            res = G.simulated_annealing(reduced=True, **kwargs)
        
        scanwidth = res[0]
        extension = res[1]
        tree_extension = extension.canonical_tree_extension()
    
    algo = DP_instance(network, tree_extension)
    res = algo._solve_MAPPD(k, include_solution)
        
    return res
