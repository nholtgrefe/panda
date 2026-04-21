"""
Scanwidth computation module for DAGs.

This module provides classes and functions for computing scanwidth of directed
acyclic graphs (DAGs), which is used in the MAPPD algorithm.

Note
----
This module is adapted from https://github.com/nholtgrefe/ComputingScanwidth. 
The scanwidth functionality will be packaged as a separate dependency in a future 
version.

Classes
-------
ScanwidthDAG
    Class for computing scanwidth of DAGs.
Extension
    Class representing an extension (ordering) of a DAG.
TreeExtension
    Class representing a tree extension of a DAG.
"""

from typing import Any

import networkx as nx

# Global table for memoization in restricted partial scanwidth
rpsw_table: dict[str, tuple[int, list[Any]]] = {}


class ScanwidthDAG:
    """
    Class for computing scanwidth of directed acyclic graphs.
    
    Scanwidth is a graph parameter that measures the complexity of a DAG
    for certain dynamic programming algorithms. This class provides methods
    to compute optimal scanwidth and corresponding extensions.
    
    Parameters
    ----------
    graph : nx.DiGraph
        The directed acyclic graph (DAG) to compute scanwidth for.
    
    Attributes
    ----------
    graph : nx.DiGraph
        The underlying DAG.
    infinity : int
        A value larger than any possible scanwidth (used as sentinel).
    _memory : bool
        Whether to use memoization for speedup.
    
    Examples
    --------
    >>> import networkx as nx
    >>> dag = nx.DiGraph([(1, 2), (1, 3), (2, 4), (3, 4)])
    >>> sw_dag = ScanwidthDAG(dag)
    >>> sw, ext = sw_dag.optimal_scanwidth()
    >>> sw >= 0
    True
    """
    
    def __init__(self, graph: nx.DiGraph) -> None:
        """
        Initialize a ScanwidthDAG.
        
        Parameters
        ----------
        graph : nx.DiGraph
            The directed acyclic graph.
        """
        self.graph = graph
        self.infinity = self.graph.number_of_edges() + 1
        self._memory = True
    
    def optimal_scanwidth(
        self,
        reduced: bool = True,
        method: int = 5,
        memory: bool = True,
    ) -> tuple[int, 'Extension'] | tuple[None, None]:
        """
        Compute the optimal scanwidth of the DAG.
        
        Parameters
        ----------
        reduced : bool, default True
            If True, first reduce the graph by splitting into s-blocks.
        method : int, default 5
            Algorithm method:
            - 5: Dynamic XP-program with increasing k (recommended)
        memory : bool, default True
            Whether to use memoization for speedup.
        
        Returns
        -------
        tuple[int, Extension] | tuple[None, None]
            The scanwidth and corresponding extension, or (None, None) if
            computation fails.
        """
        if reduced:
            best_sw, best_ext = self._reduce("optimal", reduced=False, method=method, memory=memory)
            if best_sw is None:
                return None, None
            return best_sw, best_ext
        
        self._memory = memory
        global rpsw_table
        
        if method == 5:
            rpsw_table = {}
            for i in range(1, self.infinity + 1):
                res = self.optimal_k_scanwidth(i)
                if res is False:
                    # Delete infinity values from table
                    for key in list(rpsw_table.keys()):
                        if rpsw_table[key][0] == self.infinity:
                            del rpsw_table[key]
                else:
                    sw, ext = res
                    return sw, ext
        
        return None, None
    
    def optimal_k_scanwidth(
        self,
        k: int,
        clear_table: bool = False,
    ) -> tuple[int, 'Extension'] | bool:
        """
        Solve the fixed-parameter version of scanwidth.
        
        Parameters
        ----------
        k : int
            The scanwidth parameter to check.
        clear_table : bool, default False
            If True, clear the memoization table before computation.
        
        Returns
        -------
        tuple[int, Extension] | bool
            The scanwidth and extension if rpsw <= k, False otherwise.
        """
        if clear_table:
            global rpsw_table
            rpsw_table = {}
        
        sw, sigma = self._restricted_partial_scanwidth(self.graph.nodes(), k)
        
        if sw == self.infinity:
            return False
        else:
            ext = Extension(self.graph, sigma)
            return sw, ext
    
    def _restricted_partial_scanwidth(
        self,
        vertices: set[Any],
        k: int | None = None,
        cs: bool = True,
    ) -> tuple[int, list[Any]]:
        """
        Compute restricted partial scanwidth for a set of vertices.
        
        Parameters
        ----------
        vertices : set[Any]
            Set of vertices to compute scanwidth for.
        k : int | None, optional
            If provided, only compute if rpsw <= k, otherwise return infinity.
        cs : bool, default True
            Whether to use component splitting.
        
        Returns
        -------
        tuple[int, list[Any]]
            The scanwidth and corresponding extension ordering.
        """
        # For the case that k = None, so no fixed parameter
        if k is None:
            k = self.infinity
        # For the fixed parameter version, we must use CS
        if k != self.infinity:
            assert cs is True
        
        vertex_list = list(vertices)  # Takes on the role of W
        subgraph = self.graph.subgraph(vertex_list)  # Takes on the role of G[W]
        roots = [v for v in subgraph.nodes() if subgraph.in_degree(v) == 0]
        
        # Check table
        if self._memory:
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
            components = [comp for comp in nx.weakly_connected_components(subgraph)]
            
            # If G[W] is weakly disconnected and CS is True
            if cs and len(components) > 1:
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
                    
                    rpsw1, sigma1 = self._restricted_partial_scanwidth(set(new_vertices), k, cs)
                    
                    if cs:
                        rpsw2 = delta_in_W
                    else:
                        component = self._find_component(components, rho)
                        rpsw2 = self._delta_in(component)
                    
                    rpsw_prime = max(rpsw1, rpsw2)
                    
                    if rpsw_prime < rpsw:
                        rpsw = rpsw_prime
                        sigma = sigma1 + [rho]
        
        # Save to table
        if self._memory:
            key = repr(sorted(roots))
            rpsw_table[key] = (rpsw, sigma)
        
        return rpsw, sigma
    
    def _delta_in(self, vertex_set: set[Any], sink: bool = True) -> int:
        """
        Return the in-degree of a vertex set.
        
        Parameters
        ----------
        vertex_set : set[Any]
            Set of vertices.
        sink : bool, default True
            If True, assumes vertex_set is a sink set for optimization.
        
        Returns
        -------
        int
            The in-degree of the vertex set.
        """
        res = 0
        if sink:
            if len(vertex_set) < len(self.graph.nodes()) / 2:
                # Return in-degree of W
                for v in vertex_set:
                    res = res + self.graph.in_degree(v) - self.graph.out_degree(v)
            else:
                # Return out-degree of V / W
                for v in self.graph.nodes():
                    if v not in vertex_set:
                        res = res - self.graph.in_degree(v) + self.graph.out_degree(v)
        
        if not sink:
            # If no sink set
            for (u, v) in self.graph.edges():
                if u not in vertex_set and v in vertex_set:
                    res = res + 1
        
        return res
    
    @staticmethod
    def _find_component(components: list[set[Any]], v: Any) -> set[Any]:
        """
        Find the component containing vertex v.
        
        Parameters
        ----------
        components : list[set[Any]]
            List of connected components.
        v : Any
            Vertex to find.
        
        Returns
        -------
        set[Any]
            The component containing v.
        """
        for comp in components:
            if v in comp:
                return comp
        return set()
    
    def sblocks(self) -> list[set[Any]]:
        """
        Return s-blocks of the graph.
        
        Returns a list of node sets, one for each s-block. The list is ordered
        according to a reversed DFS in the sblock-cut-tree, starting from the
        root block.
        
        Returns
        -------
        list[set[Any]]
            List of s-block node sets, ordered by DFS in sblock-cut-tree.
        """
        roots = {v for v in self.graph.nodes() if self.graph.in_degree(v) == 0}
        aux = self.graph.to_undirected()
        
        # Create auxiliary graph (root-clique for multiple roots)
        for root1 in roots:
            for root2 in roots:
                if not (root1, root2) in aux.edges() and not (root2, root1) in aux.edges() and root1 != root2:
                    aux.add_edge(root1, root2)
        
        sblock_sets = list(nx.biconnected_components(aux))
        dcut_vertices = list(nx.articulation_points(aux))
        
        # Create sblock cut tree
        sblock_cut_tree = nx.Graph()
        for v in dcut_vertices:
            sblock_cut_tree.add_node(v)
        
        rootblock_index = None
        for i, block in enumerate(sblock_sets):
            node_name = "block_" + str(i)
            if roots.issubset(block):
                rootblock_index = i
            sblock_cut_tree.add_node(node_name)
            for v in dcut_vertices:
                if v in block:
                    sblock_cut_tree.add_edge(v, node_name)
        
        # DFS in the sblock cut tree
        if rootblock_index is not None:
            sblock_order = list(nx.dfs_preorder_nodes(sblock_cut_tree, source="block_" + str(rootblock_index)))
        else:
            sblock_order = list(nx.dfs_preorder_nodes(sblock_cut_tree))
        
        # Delete cut-vertices from order and only save indices
        sblock_order = [int(name[6:]) for name in sblock_order if str(name).startswith('block_')]
        
        return [sblock_sets[i] for i in sblock_order]
    
    def _reduce(
        self,
        func_type: str,
        **kwargs: Any,
    ) -> tuple[int, 'Extension'] | tuple[None, None]:
        """
        Reduce the instance and apply scanwidth algorithm on reduced instance(s).
        
        Parameters
        ----------
        func_type : str
            Type of function to apply ('optimal').
        **kwargs : Any
            Additional keyword arguments passed to the algorithm.
        
        Returns
        -------
        tuple[int, Extension] | tuple[None, None]
            The scanwidth and extension, or (None, None) if computation fails.
        """
        sblock_sets = self.sblocks()
        sigma = []
        sw = 0
        
        # Find scanwidth on each block separately
        for sblock_set in sblock_sets:
            # This is a single edge, so we already know sw = 1
            if len(sblock_set) == 2:
                u, v = list(sblock_set)
                if (u, v) in self.graph.edges():
                    partial_sigma = [v, u]
                else:
                    partial_sigma = [u, v]
                sw = max(sw, 1)
            
            else:
                subgraph = self.graph.subgraph(sblock_set)
                roots = [v for v in subgraph.nodes() if subgraph.in_degree(v) == 0 and subgraph.out_degree(v) == 2]
                leafs = [v for v in subgraph.nodes() if subgraph.out_degree(v) == 0 and subgraph.in_degree(v) == 2]
                flow_nodes = [v for v in subgraph.nodes() if subgraph.out_degree(v) == 1 and subgraph.in_degree(v) == 1]
                
                # This means that the block is a directed cycle, so it has sw = 2
                if len(roots) == 1 and len(leafs) == 1 and len(flow_nodes) == len(subgraph.nodes()) - 2:
                    partial_sigma = list(reversed(list(nx.topological_sort(subgraph))))
                    sw = max(sw, 2)
                
                # This means the block is neither an edge or a cycle, so we apply the algorithm
                else:
                    # First contract flow-edges and save the contractions
                    history = []
                    subgraph = subgraph.copy()
                    for v in flow_nodes:
                        u = list(subgraph.predecessors(v))[0]
                        w = list(subgraph.successors(v))[0]
                        if (u, w) not in subgraph.edges():
                            subgraph.remove_node(v)
                            subgraph.add_edge(u, w)
                            history.append((w, v))
                    
                    # Run algorithm to find sigma and sw
                    S = ScanwidthDAG(subgraph)
                    
                    if func_type == 'optimal':
                        res = S.optimal_scanwidth(**kwargs)
                    else:
                        return None, None
                    
                    if res[0] is None:
                        return None, None
                    
                    block_sw, partial_ext = res
                    partial_sigma = partial_ext.sigma
                    
                    sw = max(sw, block_sw)
                    
                    # Backtrack the edge-contractions to get the partial extension for the whole block
                    history.reverse()
                    for (w, v) in history:
                        i = partial_sigma.index(w)
                        partial_sigma.insert(i + 1, v)
            
            # Add partial extension to the whole extension
            partial_sigma = [v for v in partial_sigma if v not in sigma]
            sigma = partial_sigma + sigma
        
        return sw, Extension(self.graph, sigma)


class Extension:
    """
    Class representing an extension (ordering) of a DAG.
    
    An extension is an ordering of vertices that respects the DAG structure.
    This class is used to compute scanwidth and convert to tree extensions.
    
    Parameters
    ----------
    graph : nx.DiGraph
        The underlying DAG.
    sigma : list[Any]
        The extension ordering (list of vertices).
    
    Attributes
    ----------
    graph : nx.DiGraph
        The underlying DAG.
    sigma : list[Any]
        The extension ordering.
    """
    
    def __init__(self, graph: nx.DiGraph, sigma: list[Any]) -> None:
        """
        Initialize an Extension.
        
        Parameters
        ----------
        graph : nx.DiGraph
            The underlying DAG.
        sigma : list[Any]
            The extension ordering.
        """
        self.graph = graph
        self.sigma = sigma
    
    def scanwidth(self) -> int:
        """
        Calculate the scanwidth of this extension.
        
        Returns
        -------
        int
            The maximum scanwidth over all positions in the extension.
        """
        SW_i_list = []
        
        for i in range(len(self.sigma)):
            SW_i = self.scanwidth_at_vertex_i(i)
            SW_i_list.append(SW_i)
        
        return max(SW_i_list)
    
    def scanwidth_at_vertex_i(self, i: int) -> int:
        """
        Calculate scanwidth at position i in the extension.
        
        Parameters
        ----------
        i : int
            Position in the extension.
        
        Returns
        -------
        int
            The scanwidth at position i.
        """
        left = self.sigma[0:i+1]
        
        sub = self.graph.subgraph(left)
        components = [comp for comp in nx.weakly_connected_components(sub)]
        connected_vertices = set()
        for comp in components:
            if self.sigma[i] in comp:
                connected_vertices = comp
                break
        
        SW_i = 0
        for w in connected_vertices:
            SW_i = SW_i + self.graph.in_degree(w) - self.graph.out_degree(w)
        
        return SW_i
    
    def canonical_tree_extension(self) -> 'TreeExtension':
        """
        Create the canonical tree extension with the same scanwidth.
        
        Returns
        -------
        TreeExtension
            The canonical tree extension.
        """
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
                R = set([rho[c] for c in C if rho[c] is not None])
                for r in R:
                    Gamma.add_edge(v, r)
                for u in Gamma.nodes():
                    if rho[u] in R:
                        rho[u] = v
        
        tree = TreeExtension(self.graph, Gamma)
        
        return tree


class TreeExtension:
    """
    Class representing a tree extension of a DAG.
    
    A tree extension is a tree structure that extends the DAG and is used
    for scanwidth-based dynamic programming algorithms.
    
    Parameters
    ----------
    graph : nx.DiGraph
        The underlying DAG.
    tree : nx.DiGraph
        The tree extension (must be a tree).
    
    Attributes
    ----------
    graph : nx.DiGraph
        The underlying DAG.
    tree : nx.DiGraph
        The tree extension.
    """
    
    def __init__(self, graph: nx.DiGraph, tree: nx.DiGraph) -> None:
        """
        Initialize a TreeExtension.
        
        Parameters
        ----------
        graph : nx.DiGraph
            The underlying DAG.
        tree : nx.DiGraph
            The tree extension.
        """
        self.graph = graph
        self.tree = tree
    
    def scanwidth(self) -> int:
        """
        Calculate the scanwidth of the tree extension.
        
        Returns
        -------
        int
            The maximum scanwidth over all vertices in the tree.
        """
        GW_v_list = []
        
        for v in self.tree.nodes():
            GW_v = self.scanwidth_at_vertex(v)
            GW_v_list.append(GW_v)
        
        return max(GW_v_list)
    
    def scanwidth_at_vertex(self, v: Any) -> int:
        """
        Calculate scanwidth at vertex v.
        
        Parameters
        ----------
        v : Any
            Vertex in the tree.
        
        Returns
        -------
        int
            The scanwidth at vertex v.
        """
        left = nx.descendants(self.tree, v)
        left.add(v)
        
        GW_v = 0
        for w in left:
            GW_v = GW_v + self.graph.in_degree(w) - self.graph.out_degree(w)
        
        return GW_v
    
    def scanwidth_bag(self, v: Any) -> set[tuple[Any, Any]]:
        """
        Return the set of edges in the scanwidth bag GW_v.
        
        Parameters
        ----------
        v : Any
            Vertex in the tree.
        
        Returns
        -------
        set[tuple[Any, Any]]
            Set of edges (u, w) in the scanwidth bag.
        """
        sink = nx.descendants(self.tree, v)
        sink.add(v)
        
        GW_v = set()
        for (u, w) in self.graph.edges():
            if u not in sink and w in sink:
                GW_v.add((u, w))
        
        return GW_v

