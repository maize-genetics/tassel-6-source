/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

/**
 * Created by Eli Rodgers-Melnick on 7/2/2014
 * 
 * This interface is implemented by graphs, both directed and undirected.
 * The design is inspired by the Python networkx library (https://networkx.github.io/)
 * @author Eli Rodgers-Melnick
 */
public interface Graph<T>{
    public static enum GraphType {UNDIRECTED, DIRECTED};
    /**
     * Gets an iterator over the nodes
     * @return An iterator over the nodes
     */
    public Iterator<T> nodesIter();
    /**
     * Gets a list of nodes in the graph
     * @return A list of nodes in the graph
     */
    public Collection<T> nodes();
    /**
     * Gets the number of nodes in the graph
     * @return The number of nodes in the graph
     */
    public int numberOfNodes();
    /**
     * Returns true if the graph contains the node n
     * 
     * @param n node
     * @return true if the graph contains the node n
     */
    public boolean hasNode(T n);
    /**
     * Returns true if an edge exists in the graph
     * @param u The first node in the edge
     * @param v The second node in the edge
     * @return 
     */
    public boolean hasEdge(T u, T v);
    /**
     * Gets a collection of the neighbors of a node
     * @param n The node for which you want the neighbors
     * @return The neighbors of the node
     */
    public Collection<T> neighbors(T n);
    /**
     * Gets a collection of the edges
     * @return A collection of the edges in the graph
     */
    public Collection<Map.Entry<T,T>> edges();
    /**
     * Iterates over the edges
     * @return An iterator over the edges
     */
    public Iterator<Map.Entry<T,T>> edgesIter();
    /**
     * Gets the degree of a node
     * @param n A node for which you want the degree
     * @return The degree of the node
     */
    public int degree(T n);
    /**
     * Gets the number of edges in the graph
     * @return The number of edges in the graph
     */
    public int size();
    /**
     * Gets the sum of edges in the graph
     * @param weighted true if weights should be used in the calculation
     * @return the sum of edges in the graph
     */
    public double size(boolean weighted);
}
