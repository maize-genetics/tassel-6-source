/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;

import com.google.common.collect.ImmutableSetMultimap;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;

/**
 * Class used to represent a directed graph
 * The design is inspired by the Python networkx library (https://networkx.github.io/)
 * @author Eli Rodgers-Melnick
 */
public class DirectedGraph<T> implements Graph<T> {
    // Stores nodes
    private final HashSet<T> nodes;
    // Stores adjacency of nodes
    private final ImmutableSetMultimap<T,T> adj;
    // Store predecessor adjacency for fast lookup access
    private final ImmutableSetMultimap<T,T> pred;
    // Stores the weights for edges in the graph
    private final HashMap<Tuple<T,T>,Double> wts;
    public DirectedGraph(HashSet<T> nodes, ImmutableSetMultimap<T,T> adj, ImmutableSetMultimap<T,T> pred,
            HashMap<Tuple<T,T>, Double> wts) {
        this.nodes = nodes;
        this.adj = adj;
        this.wts = wts;
        this.pred = pred;
    }
    @Override
    public Iterator<T> nodesIter() {
        return nodes.iterator();
    }

    @Override
    public Collection<T> nodes() {
        return nodes;
    }

    @Override
    public int numberOfNodes() {
        return nodes.size();
    }

    @Override
    public boolean hasNode(T n) {
        return nodes.contains(n);
    }

    @Override
    public boolean hasEdge(T u, T v) {
        return adj.containsEntry(u, v);
    }
    /**
     * Gets the successors of the node n
     * @param n A node
     * @return The collection of successors of the node (i.e. b in a->b)
     */
    @Override
    public Collection<T> neighbors(T n) {
        return adj.get(n);
    }
    /**
     * Gets the successors of the node n
     * @param n A node
     * @return The collection of successors of the node (i.e. b in a->b)
     */    
    public Collection<T> successors(T n) {
        return adj.get(n);
    }
    /**
     * Gets the predecessors of the node n
     * @param n A node
     * @return The collection of predecessors of the node (i.e. a in a->b)
     */
    public Collection<T> predecessors(T n) {
        return pred.get(n);
    }
    @Override
    public Collection<Map.Entry<T, T>> edges() {
        return adj.entries();
    }

    @Override
    public Iterator<Map.Entry<T, T>> edgesIter() {
        return adj.entries().iterator();
    }

    @Override
    public int degree(T n) {
        return (inDegree(n) + outDegree(n));
    }
    /**
     * Gets the number of nodes with an edge coming into the query node
     * @param n The query node
     * @return The number of nodes with an edge coming into the query node
     */
    public int inDegree(T n) {
        return pred.get(n).size();
    }
    /**
     * Gets the number of nodes hit by an edge coming out of the query node
     * @param n The query node
     * @return The number of nodes hit by an edge coming out of the query node
     */
    public int outDegree(T n) {
        return adj.get(n).size();
    }
    @Override
    public int size() {
        return adj.size();
    }
    @Override
    public double size(boolean weighted) {
        if (weighted) {
            double total_size = 0.;
            Iterator<Double> it = wts.values().iterator();
            while(it.hasNext()) {
                total_size += it.next();
            }
            return total_size;
        } else {
            return (double)size();
        }
    }
}    

