/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;

import com.google.common.collect.ImmutableSetMultimap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import net.maizegenetics.util.Graph.GraphType;

/**
 * Builder for graphs, directed and undirected
 * @author Eli Rodgers-Melnick
 */
public class GraphBuilder<T> {
    // Stores nodes
    private final HashSet<T> nodes;
    // Stores adjacency of nodes
    private final ImmutableSetMultimap.Builder<T,T> adj;
    // Stores the weights for edges in the graph
    private final HashMap<Tuple<T,T>,Double> wts;
    // The graph type
    private final GraphType type;
    /**
     * Generic constructor
     * @param type The type of graph to construct
     */
    public GraphBuilder(GraphType type) {
        this.type = type;
        this.nodes = new HashSet();
        this.adj = new ImmutableSetMultimap.Builder();
        this.wts = new HashMap();
    }
    public Graph build() {
        if(type == GraphType.UNDIRECTED) {
            return new UndirectedGraph(nodes, adj.build(), wts);
        }
        else {
            // Make predecessor multimap
            ImmutableSetMultimap successors = adj.build();
            final ImmutableSetMultimap.Builder<T,T> pred = new ImmutableSetMultimap.Builder();
            Iterator<Map.Entry<T,T>> successor_it = successors.entries().iterator();
            while (successor_it.hasNext()) {
                Map.Entry<T,T> item = successor_it.next();
                pred.put(item.getValue(), item.getKey());
            }
            // Instantiate directed graph
            return new DirectedGraph(nodes, successors, pred.build(), wts);
        }
    }
    /**
     * Adds a single node
     * @param n A node
     * @return The builder with the added node
     */
    public GraphBuilder addNode(T n) {
        this.nodes.add(n);
        return this;
    }
    /**
     * Adds an edge to the graph with weight 1.
     * @param u The first node in the edge
     * @param v The second node in the edge
     * @return The builder with the added edge
     */
    public GraphBuilder addEdge(T u, T v) {
        return addEdge(u, v, 1.);
    }
    /**
     * Adds an edge to the graph
     * @param u The first node in the edge
     * @param v The second node in the edge
     * @param wt The weight of the edge
     * @return The builder with the added edge
     */
    public GraphBuilder addEdge(T u, T v, double wt) {
        // Check if either node needs to be added to nodes
        if (!this.nodes.contains(u)) {this.nodes.add(u);}
        if (!this.nodes.contains(v)) {this.nodes.add(v);}
        // Add edge to builder
        adj.put(u, v);
        // Add weight
        wts.put(new Tuple(u,v), wt);
        // Add reverse edge if undirected
        if (type == GraphType.UNDIRECTED) {
            adj.put(v,u);
            wts.put(new Tuple(v,u), wt);
        }
        return this;
    }
}
