import networkx as nx
from networkx.algorithms.flow import shortest_augmenting_path

import numpy as np

import itertools, time, os, random, sys

#from plotting_module import plot_graph, plot_tree_extension

#random.seed(42)

class DAG:
    """Class for directed acyclic graphs to compute and visualize scanwidth."""
    
    def __init__(self, graph=None):
        """Takes as optional input a graph: either an nx.DiGraph, or a text file in
        edge-list format (i.e. each line contains an arc of the graph)."""
        
        if graph == None: # Empty graph
            self.graph = nx.DiGraph()
            
        elif type(graph) == str: # Initialize with graph from textfile
            self.graph = nx.DiGraph()
                     
            edge_list = []
            with open(graph, 'r') as infile:
                data = infile.readlines()
                for i in data:
                    edge = i.split()
                    edge_list.append((edge))
                    
            self.graph.add_edges_from(edge_list)
        
        else: # Initialize with given graph
            self.graph = graph
        
        self.infinity = self.graph.number_of_edges() + 1
        
        self._memory = True
        self._opt_start_time = 0
        self.func_durations = []

        
    def save_file(self, file_name):
        """Saves the graph in edge-list format at a specified location file_name."""
        
        if os.path.exists(file_name) == True:
            raise ValueError("File already exists.")
        
        f = open(file_name, "w+")
        for edge in list(self.graph.edges()):
            f.write(f"{edge[0]} {edge[1]}\n")
        f.close()
    
    
    # def plot(self, figsize=None):
    #     """Plots the graph."""
    #     plot_graph(self.graph, figsize)


    def optimal_scanwidth(self, reduced=True, method=5, memory=True):
        """Computes the scanwidth of the DAG, and also returns the corresponding
        extension. If reduced = True, we first reduce the graph.
        Method = 1, uses an exhaustive search on all possible extensions of G.
        Method = 2, uses recursive algorithm (3-partition).
        Method = 3, uses dynamic programming without component splitting.
        Method = 4, uses dynamic programming with component splitting.
        Method = 5, uses dynamic XP-program, with increasing k.
        If memory = True, we use extra space to increase speed."""
        
        if reduced == True:
            best_sw, best_ext = self._reduce("optimal", reduced=False, method=method, memory=memory)
            if best_sw is None: return None, None
            
            return best_sw, best_ext
        
        self._memory = memory
        global rpsw_table

        if method == 1:
            extension_list = [ext[::-1] for ext in nx.all_topological_sorts(self.graph)]
            SW_sigma_list = []
            
            for sigma in extension_list:
                ext = Extension(self.graph, sigma)
                SW_sigma = ext.scanwidth()
                SW_sigma_list.append(SW_sigma)

            best_sw = min(SW_sigma_list)
            best_sigma_index = SW_sigma_list.index(best_sw)
            best_ext = Extension(self.graph, extension_list[best_sigma_index])
            
            return best_sw, best_ext
    
        elif method in [2, 3, 4]:
            rpsw_table = {}
            
            if method == 4:
                sw, sigma = self._restricted_partial_scanwidth(self.graph.nodes(), cs=True)
            elif method == 2:
                sw, sigma = self._partial_scanwidth(set(), set(self.graph.nodes()), set())
            elif method == 3:
                sw, sigma = self._restricted_partial_scanwidth(self.graph.nodes(), cs=False)
                        
            ext = Extension(self.graph, sigma)
            
            return sw, ext
            
            
        elif method == 5:
            rpsw_table = {}
                        
            for i in range(1, self.infinity + 1):
                res = self.optimal_k_scanwidth(i)
                    
                if res == False:
                    # Delete infinity values from table, they have to be recalculated.
                    for key in list(rpsw_table):
                        if rpsw_table[key][0] == self.infinity:
                            del rpsw_table[key]                 
                else:                    
                    sw, ext = res
                    return sw, ext
                    
        
    def optimal_k_scanwidth(self, k, clear_table=False):
        """Solves the fixed parameter version of scanwidth with rpsw_XP algorithm. If we want to call
        this function by itself, clear_table must be True."""
        
        if clear_table == True:
            global rpsw_table
            rpsw_table = {}
            
        sw, sigma = self._restricted_partial_scanwidth(self.graph.nodes(), k)
        
        if sw == self.infinity:
            return False
        else:
            ext = Extension(self.graph, sigma)
            #assert ext.scanwidth() == sw
            
            return sw, ext
        
    def _restricted_partial_scanwidth(self, vertices, k=None, cs=True):
        """If k == None, outputs the rpsw and corresponding extension for the set vertices.
        Else, only does this if rpsw<=k, otherwise  it outputs 'infinity'.
        cs = True, uses component splitting."""
                
        # For the case that k = None, so no fixed parameter
        if k is None: k = self.infinity
        # For the fixed parameter version, we must use CS
        if k != self.infinity: assert cs == True
        
        vertex_list = list(vertices) # Takes on the role of W
        subgraph = self.graph.subgraph(vertex_list) # Takes on the role of G[W]      
        roots = [v for v in subgraph.nodes() if subgraph.in_degree(v) == 0] ##
        
        # Check table
        if self._memory == True:
            key = repr(sorted(roots))
            if key in rpsw_table.keys():
                return rpsw_table[key]

        delta_in_W = self._delta_in(vertices)
        
        # Initialize
        rpsw = self.infinity
        sigma = []
                
        # If |W| = 1 and delta_in(W) <= k
        if len(vertex_list) == 1 and delta_in_W <= k:
            rpsw = delta_in_W
            sigma = [vertex_list[0]]
        
        else:
            components = [comp for comp in nx.weakly_connected_components(subgraph)] ##
            
            # If G[W] is weakly disconnected and CS is True
            if cs == True and len(components) > 1:
                rpsw_list = []
                for U_i in components:
                    rpsw_i, sigma_i = self._restricted_partial_scanwidth(U_i, k, cs)
                    rpsw_list.append(rpsw_i)
                    sigma = sigma + sigma_i
                    
                rpsw = max(rpsw_list)
                        
            # If |W| > 1 and delta_in(W) <= k
            elif len(vertex_list) > 1 and delta_in_W <= k:
                # Minimize over the roots of the subgraph
                for rho in roots:
                    new_vertices = vertex_list.copy()
                    new_vertices.remove(rho)

                    rpsw1, sigma1 = self._restricted_partial_scanwidth(new_vertices, k, cs)
                    
                    if cs == True:
                        rpsw2 = delta_in_W
                    else:
                        component = self._find_component(components, rho)
                        rpsw2 = self._delta_in(component)
                        
                    rpsw_prime = max(rpsw1, rpsw2)
                    
                    if rpsw_prime < rpsw:
                        rpsw = rpsw_prime
                        sigma = sigma1 + [rho]

        # Save to table
        if self._memory == True: rpsw_table[key] = (rpsw, sigma)

        return rpsw, sigma
    
    def _partial_scanwidth(self, L, W, R):
        """Returns the psw and corresponding extension for the ordered 3-partition L, W, R."""

        # Initialize
        psw = self.infinity
        sigma = []

        # If |W| = 1
        if len(W) == 1:
            L_U_W = W.union(L)
            (w, ) = W # unpacks the unique element in W
            subgraph = self.graph.subgraph(L_U_W) # Takes on the role of G[L U W]
            components = [comp for comp in nx.weakly_connected_components(subgraph)]
            component = self._find_component(components, w)
                        
            psw = self._delta_in(component)
        
            sigma = [w]
        
        # If |W| > 1
        elif len(W) > 1:
            size = len(W) // 2
            for W_prime in map(set, itertools.combinations(W, size)):
                
                #self._check_timeout()
                    
                # Check if W' sqsubseteq W
                check = len([(u,v) for u in W_prime for v in W.difference(W_prime) if (u,v) in self.graph.edges()])
                if check == 0:
                    psw_1, sigma_1 = self._partial_scanwidth(L, W_prime, R.union(W.difference(W_prime)))
                    psw_2, sigma_2 = self._partial_scanwidth(L.union(W_prime), W.difference(W_prime), R)
                    psw_prime = max(psw_1, psw_2)
                    if psw_prime < psw:
                        psw = psw_prime
                        sigma = sigma_1 + sigma_2
        
        return psw, sigma
        
    
    def cut_splitting_heuristic(self, reduced=True):
        """Heuristic to compute the scanwidth, using recursive cut-splitting."""
        
        if reduced == True:
            sw, ext = self._reduce("cut", reduced=False)
            return sw, ext
        
        sigma = self._recursive_cut_splitting()
        ext = Extension(self.graph, sigma)
        sw = ext.scanwidth()
        
        return sw, ext
    
    def _recursive_cut_splitting(self):
        """Recursively splits the graph at a DAG-cut and merges the rest."""
        
        edge_cut, source_set, sink_set = self._minimum_DAG_cut(non_trivial=True)

        # If there is no non-trivial cut anymore, just pick any ordering
        if edge_cut is None:
            return list(reversed(list(nx.topological_sort(self.graph))))
        
        subgraph_top = self.graph.subgraph(source_set)
        subgraph_down = self.graph.subgraph(sink_set)
        subgraph_top = subgraph_top.copy()
        subgraph_down = subgraph_down.copy()        
        
        # Add super leaf and super root
        super_leaf = "super_leaf_" + str(random.getrandbits(64))
        super_root = "super_root_" + str(random.getrandbits(64))
        subgraph_top.add_node(super_leaf, merged=True)
        subgraph_down.add_node(super_root, merged=True)
                
        # Add edges to super nodes (or add weight)
        for (u,v) in edge_cut:
            
            if 'weight' not in self.graph[u][v].keys():
                weight = 1
            else:
                weight = self.graph[u][v]['weight']
                        
            if subgraph_top.has_edge(u, super_leaf):
                subgraph_top[u][super_leaf]['weight'] += weight
            else:
                subgraph_top.add_edge(u, super_leaf, weight=weight)
                
            if subgraph_down.has_edge(super_root, v):
                subgraph_down[super_root][v]['weight'] += weight
            else:
                subgraph_down.add_edge(super_root, v, weight=weight)            
        
        top = DAG(subgraph_top)
        down = DAG(subgraph_down)

        # Recurse
        sigma_top = top._recursive_cut_splitting()
        sigma_down = down._recursive_cut_splitting()
        sigma_top.remove(super_leaf)
        sigma_down.remove(super_root)
                    
        return sigma_down + sigma_top
        
    
    def _minimum_DAG_cut(self, non_trivial=True):
        """Finds the smallest DAG cut. If non-trivial is True, it only looks for 
        cuts s.t. both sides of the cut do not consists of a single merged node. Outputs 
        the edges in the cut, and the node sets of both sides."""
        
        aux_graph = self.graph.copy()
        
        # Give all edges weight of 1
        for (u,v) in aux_graph.edges():
            if 'weight' not in aux_graph[u][v].keys():
                aux_graph[u][v]['weight'] = 1
                
        # Add reverse edges with weight infinity
        aux_graph.add_weighted_edges_from(list((v, u, self.infinity) for (u,v) in aux_graph.edges))

        roots = [v for v in self.graph.nodes() if self.graph.in_degree(v) == 0]
        leafs = [v for v in self.graph.nodes() if self.graph.out_degree(v) == 0]
        
        cut = self.infinity
        source_set = []
        sink_set = []
        size_diff = self.infinity + 1 # We prefer balanced seperators if there is a tie
        
        # Find all DAG cuts
        if non_trivial == False:
            for root in roots:
                for leaf in leafs:
                    a, (b, c) = nx.minimum_cut(aux_graph, root, leaf, capacity='weight',flow_func=shortest_augmenting_path)
                    d = abs(len(b) - len(c))
                    if a < cut or (a == cut and d < size_diff):
                        cut, source_set, sink_set, size_diff = a, b, c, d
                
        elif non_trivial == True:
            
            for root in roots:
                # If the root was merged, it is a trivial cut
                if 'merged' in self.graph.nodes[root].keys():
                    root_childs = list(self.graph.successors(root))
                else:
                    root_childs = [root]
                

                for leaf in leafs:
                    # If the leaf was merged, it is a trivial cut
                    if 'merged' in self.graph.nodes[leaf].keys():
                        leaf_parents = list(self.graph.predecessors(leaf))
                    else:
                        leaf_parents = [leaf]
                    
        
                    for root_child in root_childs:
                        for leaf_parent in leaf_parents:
                            
                            if root_child == leaf_parent: continue # Not a possible cut
                            a, (b, c) = nx.minimum_cut(aux_graph, root_child, leaf_parent, capacity='weight', flow_func=shortest_augmenting_path)
                            d = abs(len(b) - len(c))

                            if a < cut or (a == cut and d < size_diff):
                                cut, source_set, sink_set, size_diff = a, b, c, d
        
        if cut == self.infinity:
            # No DAG cut has been found
            return None, None, None
        
        # Find corresponding edges
        edge_cut = []
        for (u,v) in self.graph.edges():
            if u in source_set and v in sink_set:
                edge_cut.append((u,v))
        
        return edge_cut, source_set, sink_set

    
    def simulated_annealing(self, max_iter=100, p_in=0.9, p_stop=0.01, init_ext='greedy', reduced=False, verbose=True, seed=42):
        """Uses the simulated annealing heuristic to compute the scanwidth. Parameters:
            max_iter: maximum number of iterations;
            p_in: the initial acceptance probability;
            p_stop: the final acceptance probability;
            init_ext: initial extension to be used. Can either be an object from the Extension class, or a string ('cut', 'greedy' or 'random')
                to indicate which method should be used to get the initial extension.
            reduced: whether the graph should first be reduced.
            """

        if reduced == True: # Split in blocks and do algorithm per block (thus |blocks|*max_iter iterations)
            best_sw, best_ext = self._reduce("annealing", max_iter=max_iter, p_in=p_in, p_stop=p_stop, init_ext=init_ext, reduced=False, verbose=verbose, seed=seed)            
            return best_sw, best_ext
        
        rng = np.random.RandomState(seed)
        
        if type(init_ext) == Extension:
            extension = init_ext
        elif init_ext == 'greedy':
            extension = self.greedy_heuristic()[1]
        elif init_ext == 'cut':
            extension = self.cut_splitting_heuristic()[1]
        elif init_ext == 'random':
            extension = self.random_extension()[1]
        
        tree_extension = extension.canonical_tree_extension()
        tree = tree_extension.tree
        
        # Initial scanwidth values
        sw_values = {}
        for v in self.graph.nodes:
            sw_values[v] = self._delta_in(nx.descendants(tree, v) | {v})
        
        # Algorithm initialization
        epoch_counter = 0    
        iteration_counter = 0
        best_tree = tree.copy()
        best_sw_values = sw_values.copy()
        best_sw = max(best_sw_values.values())
        epoch_length = len(self.graph.nodes())
        
        # Find out initial and stop temperature
        random_scanwidths = []
        for i in range(100):
            t = TreeExtension(self.graph, tree)
            s = t._random_neighbour_scanwidth(sw_values, rng)
            if s is not None: 
                random_scanwidths.append(max(s[0].values()))
                
        random_scanwidths = [elt - best_sw for elt in random_scanwidths if elt - best_sw > 0]
        if len(random_scanwidths) == 0: delta_f = 1
        else: delta_f = sum(random_scanwidths) / len(random_scanwidths)
        init_temp = - delta_f / np.log(p_in)
        stop_temp = - delta_f / np.log(p_stop)
        temp = init_temp
        vals = [max(sw_values.values())]
        
        while temp >= stop_temp:
            
            # Per epoch, we do |V| iterations of the following
            for i in range(epoch_length):
                
                # Find a random neighbour and its scanwidth
                t = TreeExtension(self.graph, tree)
                res = t._random_neighbour_scanwidth(sw_values, rng)
                if res is None:
                    print("Only one possible tree extension")
                    return max(sw_values.values()), extension, [max(sw_values.values()) for i in range(max_iter + 1)]
                
                new_sw_values, vertex, parent, connected_vertices = res
                    
                # Calculate metropolis acceptance criterion
                diff = max(new_sw_values.values()) - max(sw_values.values())
                metropolis = np.exp(-diff / temp)
                
                # If improvement or random acceptance, accept the new tree
                if diff < 0 or rng.random() < metropolis:
                    
                    # Update values
                    sw_values = new_sw_values
                    
                    # Create new tree
                    succ = list(tree.successors(vertex))
                    grandparents = list(tree.predecessors(parent))
                    
                    for child in succ:
                        if child in connected_vertices:
                            tree.remove_edge(vertex, child)
                            tree.add_edge(parent, child)
                            
                    tree.remove_edge(parent, vertex)
                    
                    if len(grandparents) != 0:
                        grandparent = list(tree.predecessors(parent))[0]
                        tree.remove_edge(grandparent, parent)
                        tree.add_edge(grandparent, vertex)
                    
                    tree.add_edge(vertex, parent)
                    
                    # Keep track of best tree so far
                    if max(sw_values.values()) <= best_sw:
                        best_tree = tree.copy()
                        best_sw_values = sw_values.copy()
                        best_sw = max(best_sw_values.values())

                #trex = TreeExtension(self.graph, tree)
                #assert max(sw_values.values()) == trex.scanwidth()
                iteration_counter = iteration_counter + 1
                
                if verbose and iteration_counter % 100 == 0:
                    sys.stdout.write("\r" + f"current scanwidth: {max(sw_values.values())};   iteration {iteration_counter};    best scanwidth: {best_sw}")
                    time.sleep(1e-40)
            
            vals.append(max(sw_values.values()))
            epoch_counter = epoch_counter + 1
            
            # Change temperature
            alpha = (stop_temp / init_temp )**(1/(max_iter -1))
            temp = alpha * temp            
        
        best_tree_extension = TreeExtension(self.graph, best_tree)

        best_extension = best_tree_extension.to_extension()
        return best_sw, best_extension, vals


    def greedy_heuristic(self, reduced=True):
        """Greedy heuristic to find the scanwidth, also returns corresponding extension."""
        
        if reduced == True:
            sw, ext = self._reduce("greedy", reduced=False)
            return sw, ext
        
        
        S = set(self.graph.nodes())
        T = set()
        sigma = []
        sw = 0
        components_T = [] # Speedup: keeps track of the components of G[T]
        
        while len(S) > 0:
            x = None
            sw_x = self.infinity
            sub = self.graph.subgraph(S)
            leafs = [v for v in S if sub.out_degree(v) == 0]
            
            connection_dict = {} # Keeps track of the components of G[T u {l}]
            for l in leafs:
                
                c = set(self.graph.successors(l))
                connected_vertices = set()
                for comp in components_T:
                    if not c.isdisjoint(comp): # If components has a child of l, it is conencted
                        connected_vertices = connected_vertices | comp
                connected_vertices.add(l)
                connection_dict[l] = connected_vertices
                
                s = self._delta_in(connected_vertices)
                
                if s < sw_x:
                    x = l
                    sw_x = s
            
            S.remove(x)
            T.add(x)
            sigma = sigma + [x]
            sw = max(sw, sw_x)
            
            # New components of G[T]
            components_T = [comp for comp in components_T if not comp.issubset(connection_dict[x])]
            components_T.append(connection_dict[x])
        
        ext = Extension(self.graph, sigma)

        return sw, ext


    def random_extension(self, seed=42):
        """Create a random extension of the DAG, using the greedy structure,
        and also retunrs the corresponding scanwidth."""       
        
        rng = np.random.RandomState(seed)
        
        S = set(self.graph.nodes())
        T = set()
        sigma = []
        sw = 0
        components_T = [] # Speedup: keeps track of the components of G[T]
        
        while len(S) > 0:
            sub = self.graph.subgraph(S)
            leafs = [v for v in S if sub.out_degree(v) == 0]
            
            l = rng.choice(leafs)
                
            c = set(self.graph.successors(l))
            connected_vertices = set()
            for comp in components_T:
                if not c.isdisjoint(comp): # If components has a child of l, it is conencted
                    connected_vertices = connected_vertices | comp
            connected_vertices.add(l)
            
            s = self._delta_in(connected_vertices)
            
            S.remove(l)
            T.add(l)
            sigma = sigma + [l]
            sw = max(sw, s)
            
            # New components of G[T]
            components_T = [comp for comp in components_T if not comp.issubset(connected_vertices)]
            components_T.append(connected_vertices)
        
        ext = Extension(self.graph, sigma)

        return sw, ext
    
    
    def sblocks(self):
        """Returns a list of nodesets, one for each s-block of the graph. The list is ordered,
        according to a reversed DFS in the sblock-cut-tree, starting from the rootblock. (Thus the rootblock is 
        at the end.) This allows splitting the scanwidth problem over these blocks and then attaching them 
        again in this ordering."""
        
        roots = {v for v in self.graph.nodes() if self.graph.in_degree(v) == 0}
        aux = self.graph.to_undirected()

        # Create auxilliary graph, note that if G is rooted, this does not do anything.
        for root1 in roots:
            for root2 in roots:
                # Make a root-clique
                if not (root1, root2) in aux.edges() and not (root2, root1) in aux.edges and root1 != root2:
                    aux.add_edge(root1, root2)
        
        sblock_sets = list(nx.biconnected_components(aux))
        dcut_vertices = list(nx.articulation_points(aux))
                        
        # Create sblock cut tree
        sblock_cut_tree = nx.Graph()
        for v in dcut_vertices: # add the directed cut-vertices
            sblock_cut_tree.add_node(v)
        
        rootblock_index = None
        for i, block in enumerate(sblock_sets): # add the s-blocks
            node_name = "block_"+str(i)
            if roots.issubset(block): # keep track of the rootblock
                rootblock_index = i
            sblock_cut_tree.add_node(node_name)
            for v in dcut_vertices:
                if v in block:
                    sblock_cut_tree.add_edge(v, node_name)
        
        # DFS in the sblock cut tree
        sblock_order = list(nx.dfs_preorder_nodes(sblock_cut_tree, source="block_"+str(rootblock_index)))
        
        # Delete cut-vertices from order and only save indices:
        sblock_order = [int(name[6:]) for name in sblock_order if str(name).startswith('block_')]
        
        return [sblock_sets[i] for i in sblock_order]



    def _reduce(self, func_type, **kwargs):
        """Function that first reduces the instance of the DAG, then it applies any
        scanwidth-algorithm or -heuristic on the reduced instance(s)."""
        
        sblock_sets = self.sblocks()
        sigma = []
        sw = 0

        # We find the scanwidth on each block seperately
        for sblock_set in sblock_sets:
            # This is a single edge, so we already know sw = 1
            if len(sblock_set) == 2: 
                u, v = sblock_set
                if (u, v) in self.graph.edges():
                    partial_sigma = [v, u]
                else:
                    partial_sigma = [u, v]
                sw = max(sw, 1)
            
            else:
                subgraph = self.graph.subgraph(sblock_set)
                roots = [v for v in subgraph.nodes() if subgraph.in_degree(v)==0 and subgraph.out_degree(v)==2]
                leafs = [v for v in subgraph.nodes() if subgraph.out_degree(v)==0 and subgraph.in_degree(v)==2]
                flow_nodes = [v for v in subgraph.nodes() if subgraph.out_degree(v)==1 and subgraph.in_degree(v)==1]

                # This means that the block is a directed cycle, so it has sw = 2
                if len(roots) == 1 and len(leafs) == 1 and len(flow_nodes) == len(subgraph.nodes()) - 2:
                    partial_sigma = list(reversed(list(nx.topological_sort(subgraph)))) # any ordering is okay
                    sw = max(sw, 2) 
                
                # This means the block is neither an edge or a cycle, so we apply the algorithm
                else: 
                    # First contract flow-edges and save the contractions
                    history = []
                    subgraph = subgraph.copy()
                    for v in flow_nodes:
                        u = list(subgraph.predecessors(v))[0]
                        w = list(subgraph.successors(v))[0]
                        if (u,w) not in subgraph.edges():
                            subgraph.remove_node(v)
                            subgraph.add_edge(u, w)
                            history.append((w, v))
                    
                    # Run algorithm to find sigma and sw                    
                    S = DAG(subgraph)

                    if func_type == 'optimal':
                        res = S.optimal_scanwidth(**kwargs)
                    elif func_type == 'cut':
                        res = S.cut_splitting_heuristic(**kwargs)
                    elif func_type == 'greedy':
                        res = S.greedy_heuristic(**kwargs)
                    elif func_type == 'annealing':
                        res = S.simulated_annealing(**kwargs)
                    
                    block_sw, partial_ext = res
                    if block_sw is None: return None, None
                    partial_sigma = partial_ext.sigma
                    
                    sw = max(sw, block_sw)
                    
                    # Backtrack the edge-contractions to get the partial extension for the whole block
                    history.reverse()
                    for (w, v) in history:
                        i = partial_sigma.index(w)
                        partial_sigma.insert(i+1, v)
            
            # Add partial extension to the whole extension
            partial_sigma = [v for v in partial_sigma if v not in sigma]
            sigma = partial_sigma + sigma
            
        return sw, Extension(self.graph, sigma)

    def nr_of_leaves(self):
        """Returns the number of leaves of the graph."""
        return len([u for u in self.graph.nodes() if self.graph.out_degree(u) == 0])

    def _delta_in(self, vertex_set, sink=True):
        """Returns the indegree of vertex_set. Setting sink to True, if we already know
        that it is a sinkset will speed up computation."""
        res = 0
        if sink == True:
            if len(vertex_set) < len(self.graph.nodes()) / 2: # Return indegree of W
                for v in vertex_set:
                    res = res + self.graph.in_degree(v) - self.graph.out_degree(v)
            else: # Return outdegree of V / W
                for v in self.graph.nodes:
                    if v not in vertex_set:
                        res = res - self.graph.in_degree(v) + self.graph.out_degree(v)
            
        if sink == False: # If no sinkset
            for (u, v) in self.graph.edges():
                if u not in vertex_set and v in vertex_set:
                    res = res + 1            
            
        return res


    @staticmethod
    def _find_component(components, v): 
        """Returns the component in the list components, that contains v."""
        for comp in components:
            if v in comp:
                return comp


class Network(DAG):
    """Subclass of the DAG-class, specifically for networks."""
    
    def reticulation_number(self):
        """Returns the reticulation number of the network."""
        rets = [v for v in self.graph.nodes() if self.graph.in_degree(v) >= 2]
        res = 0
        for ret in rets:
            res = res + self.graph.in_degree(ret) - 1
        return res
    
    def level(self):
        """Returns the level of the network."""
        block_list = self.sblocks()
        
        rets = []
        for block in block_list:
            sub = self.graph.subgraph(block)
            sub_net = Network(sub)
            ret = sub_net.reticulation_number()
            rets.append(ret)
        
        return max(rets)
    
    def is_binary(self):
        """Checks if the network is binary."""
        for v in self.graph.nodes():
            if (self.graph.in_degree(v) == 0 and self.graph.out_degree(v) == 2) or \
                (self.graph.in_degree(v) == 0 and self.graph.out_degree(v) == 1) or \
                (self.graph.in_degree(v) == 1 and self.graph.out_degree(v) == 0) or \
                (self.graph.in_degree(v) == 2 and self.graph.out_degree(v) == 1) or \
                (self.graph.in_degree(v) == 1 and self.graph.out_degree(v) == 2):
                    continue
            else:
                return False
            
        return True
    
    def nr_of_reticulations(self):
        """Returns the number of reticulations of the network."""
        return len([v for v in self.graph.nodes() if self.graph.in_degree(v) >= 2])
    


class TreeExtension:
    """Class for a tree extension of a graph."""
    
    def __init__(self, graph, tree):
        """Takes as input the corresponding graph as an nx.DiGraph, and the
        tree extension as an nx.DiGraph."""
        self.graph = graph
        self.tree = tree
    
    
    def scanwidth(self):
        """Calculates the scanwidth of the tree extension for the DAG."""
       
        GW_v_list = []
        
        for v in self.tree.nodes():
            GW_v = self.scanwidth_at_vertex(v)
            GW_v_list.append(GW_v)
        
        return max(GW_v_list)


    def scanwidth_at_vertex(self, v):
        """Calculates the scanwidth of the tree extension at vertex v."""
                    
        left = nx.descendants(self.tree, v)
        left.add(v)
        
        GW_v = 0
        for w in left:
            GW_v = GW_v + self.graph.in_degree(w) - self.graph.out_degree(w)
            
        return GW_v
    
    def scanwidth_bag(self, v):
        """Returns the set of edges in the scanwidth bag GW_v."""
        
        sink = nx.descendants(self.tree, v)
        sink.add(v)
        
        GW_v = set()
        for (u,v) in self.graph.edges:
            if u not in sink and v in sink:
                GW_v.add((u,v))
                
        return GW_v
    
    
    def to_extension(self):
        """Returns an extension of the tree. If the tree extension is canonical, this extension
        gives the same scanwidth."""
        
        sigma = list(reversed(list(nx.topological_sort(self.tree))))
        return Extension(self.graph, sigma)
    
    
    def is_canonical(self):
        """Reports whether the tree extension is canonical."""
        
        for v in self.tree.nodes():
            left = nx.descendants(self.tree, v)
            left.add(v)
            if nx.is_weakly_connected(self.graph.subgraph(left)) == False:
                return False
            
        return True
                
    
    # def plot(self, fancy=True, figsize=None):
    #     """Plots the tree extension with the underlying graph. For larger graphs,
    #     it might be necessary to adjust the figsize for a presentable plot."""
        
    #     off = (6, 3)
    #     plot_tree_extension(self.graph, self.tree, fancy, off, figsize)
        
    
    def _random_neighbour_scanwidth(self, sw_values, rng):
        """Finds the scanwidth of a random neighbour of the tree extension.
        Takes as input a dictionary of the scanwidth values of the current tree
        and a random number generator rng. Returns the a dictionairy of new
        scanwidth-values, the two vertices that are swapped and a set of vertices
        that are connected to the parent in the graph G[1..parent]."""
        
        # Choose random vertex
        possible_choices = []
        for v in self.graph.nodes:
            pred = list(self.tree.predecessors(v))
            if len(pred) == 0: # v is the root and we can not move the root up
                continue
            elif (pred[0], v) in self.graph.edges: # If an edge from pred to vertex, we can not swap them
                continue
            else:
                possible_choices.append(v)
                
        if len(possible_choices) == 0:
            return None
        
        # Vertex to be swapped above
        vertex = rng.choice(possible_choices)
        parent = list(self.tree.predecessors(vertex))[0]
        
        # Calculate the new sw_values (only parent and vertex change)
        new_sw_values = sw_values.copy()
        
        # New value for vertex
        new_sw_values[vertex] = sw_values[parent]
        
        # New value for parent
        desc = list(nx.descendants(self.tree, parent))
        desc.remove(vertex)
        desc.append(parent)
        
        connected_vertices = {parent}
        for child in self.tree.successors(parent):
            if child != vertex:
                connected_vertices  =  connected_vertices | nx.descendants(self.tree, child) | {child}
        for child in self.tree.successors(vertex):
            sinkset = nx.descendants(self.tree, child) | {child}
            connected = False
            
            for u in sinkset:
                if (parent, u) in self.graph.edges():
                    connected = True
                    break
            
            if connected == True:
                connected_vertices = connected_vertices | nx.descendants(self.tree, child) | {child}

        new_sw_values[parent] = self._delta_in(connected_vertices)
        
        return new_sw_values, vertex, parent, connected_vertices
        
    def _delta_in(self, vertex_set, sink=True):
        """Returns the indegree of vertex_set. Setting sink to True, if we already know
        that it is a sinkset will speed up computation."""
        res = 0
        if sink == True:
            if len(vertex_set) < len(self.graph.nodes()) / 2: # Return indegree of W
                for v in vertex_set:
                    res = res + self.graph.in_degree(v) - self.graph.out_degree(v)
            else: # Return outdegree of V / W
                for v in self.graph.nodes:
                    if v not in vertex_set:
                        res = res - self.graph.in_degree(v) + self.graph.out_degree(v)
            
        if sink == False: # If no sinkset
            for (u, v) in self.graph.edges():
                if u not in vertex_set and v in vertex_set:
                    res = res + 1            
            
        return res
        
class Extension:
    """Class for an extension of a graph, starting with the leafnodes."""
    
    def __init__(self, graph, sigma):
        """Takes as input the corresponding graph as an nx.DiGraph, and a list
        of nodes for the extension. Can also read sigma from a text file."""
        self.graph = graph
        
        if type(sigma) == str: # Initialize with sigma from textfile
            self.sigma = []
            with open(sigma, 'r') as infile:
                data = infile.readlines()
                for i in data:
                    v = i.split()[0]
                    self.sigma.append(v)
                            
        else: # Initialize with given sigma
            self.sigma = sigma
            
    def save_file(self, file_name):
        """Saves the extension in a file at a specified location. Each line contains one vertex of the extension."""
        
        if os.path.exists(file_name) == True:
            raise ValueError("File already exists.")
        
        f = open(file_name, "w+")
        for v in self.sigma:
            f.write(f"{v}\n")
        f.close()
        
        
    def scanwidth(self):
        """Calculates the scanwidth of the extension sigma for the DAG."""
        
        SW_i_list = []
        
        for i in range(len(self.sigma)):
            SW_i = self.scanwidth_at_vertex_i(i)
            SW_i_list.append(SW_i)
        
        return max(SW_i_list)


    def scanwidth_at_vertex_i(self, i, position=True):
        """Calculates the size of the set SW of the extension sigma at position i. If
        position is False, we find the scanwidth at the vertex i (i.e. the node-name)."""
        if position == False:
            i = self.sigma.index(i)
        
        left = self.sigma[0:i+1]
        
        sub = self.graph.subgraph(left)
        components = [comp for comp in nx.weakly_connected_components(sub)]
        for comp in components:
            if self.sigma[i] in comp:
                connected_vertices = comp
                
        SW_i = 0
        for w in connected_vertices:
            SW_i = SW_i + self.graph.in_degree(w) - self.graph.out_degree(w)
        
        return SW_i
    
    def canonical_tree_extension(self):
        """Creates the canonical tree extension with the same scanwidth as sigma."""
        
        # Initialize
        Gamma = nx.DiGraph()
        sig = self.sigma.copy()
        rho = {node: None for node in self.graph.nodes()}
        
        while len(sig) > 0:
            v = sig[0]
            sig.remove(v)
            C = list(self.graph.successors(v))
            Gamma.add_node(v)
            rho[v] = v
            
            if len(C) > 0:
                R = set([rho[c] for c in C])
                for r in R:
                    Gamma.add_edge(v, r)
                for u in Gamma.nodes():
                    if rho[u] in R:
                        rho[u] = v
        
        tree = TreeExtension(self.graph, Gamma)
        
        return tree
    
    
    def is_extension(self):
        """Checks if sigma is indeed an extension of the graph."""
        
        seen = []
        for i in range(len(self.sigma)):
            v = self.sigma[i]
            succ = list(self.graph.successors(v))
            for w in succ:
                if w  not in seen:
                    return False
            seen.append(v)
        return True
    
    def _find_component(components, v): 
        """Returns the component in the list components, that contains v."""
        for comp in components:
            if v in comp:
                return comp

