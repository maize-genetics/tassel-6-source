/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableSetMultimap;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;

/**
 * Class used to represent an undirected graph. 
 * The design is inspired by the Python networkx library (https://networkx.github.io/)
 * @author Eli Rodgers-Melnick
 */
public class UndirectedGraph<T> implements Graph<T> {
    // Stores nodes
    private final HashSet<T> nodes;
    // Stores adjacency of nodes
    private final ImmutableSetMultimap<T,T> adj;
    // Stores the weights for edges in the graph
    private final HashMap<Tuple<T,T>,Double> wts;
    protected class UndirectedIterator implements Iterator {
        // edges that have already been added
        private final HashSet<Tuple<T,T>> added;
        // The raw iterator
        private final Iterator<Map.Entry<T,T>> it;
        // The next entry to give
        private Map.Entry<T,T> on_deck;
        public UndirectedIterator(Iterator it) {
            added = new HashSet();
            this.it = it;
            if(this.it.hasNext()) {
                on_deck = this.it.next();
            } else {
                on_deck = null;
            }
        }
        @Override
        public boolean hasNext() {
            if (on_deck != null) {
                return true;
            } else {return false;}
        }

        @Override
        public Map.Entry<T, T> next() {
            if (on_deck == null) {
                throw new NoSuchElementException();
            } else { 
                // Set the on deck to be returned
                Map.Entry<T,T> return_val = on_deck;
                // Add this to the added
                added.add(new Tuple(return_val.getKey(), return_val.getValue()));
                // Find the next edge that hasn't been added
                if (it.hasNext())
                    while (it.hasNext()) {
                        on_deck = it.next();
                        if (added.contains(new Tuple(on_deck.getValue(), on_deck.getKey()))) {
                            on_deck = null;
                        } else {break;}
                    }
                else {
                    on_deck = null;
                }
                return return_val; 
            }
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Immutable"); //To change body of generated methods, choose Tools | Templates.
        }
        
    }
    protected class UndirectedSet<E> extends HashSet<E> {
        boolean instantiated = false;
        public UndirectedSet(ImmutableSet c) {
            super(c);
            instantiated = true;
        }
        @Override
        public boolean add(E e) {
            if (instantiated) {
                throw new UnsupportedOperationException("Immutable");
            } else {
                return super.add(e);
            }

        }
        @Override
        public UndirectedIterator iterator() {
            return new UndirectedIterator(super.iterator());
        }
        @Override
        public Object[] toArray() {
            Object[] arr = new Object[this.size()];
            int i = 0;
            UndirectedIterator it = this.iterator();
            while (it.hasNext()) {
                arr[i] = it.next();
                i++;
            }
            return arr;
        }
        @Override
        public int size() {
            return super.size()/2;
        }
        @Override
        public void clear() {
            throw new UnsupportedOperationException("Immutable");
        }
        @Override
        public boolean remove(Object o) {
            throw new UnsupportedOperationException("Immutable");
        }
        @Override
        public boolean removeAll(Collection<?> c) {
            throw new UnsupportedOperationException("Immutable");
        }
        @Override
        public boolean addAll(Collection<? extends E> c) {
            if (instantiated) {
                throw new UnsupportedOperationException("Immutable");
            } else {
                return super.addAll(c);
            }
        }
        @Override
        public boolean retainAll(Collection<?> c) {
            throw new UnsupportedOperationException("Immutable");
        }
    }
    /**
     * Constructs an undirected graph
     * @param nodes Set of nodes
     * @param adj Multimap of adjacency between nodes
     * @param wts Hashmap of an edge (represented by a Map Entry) and weights
     */
    public UndirectedGraph(HashSet<T> nodes, ImmutableSetMultimap<T,T> adj, HashMap<Tuple<T,T>, Double> wts) {
        this.nodes = nodes;
        this.adj = adj;
        this.wts = wts;
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

    @Override
    public Collection<T> neighbors(T n) {
        return adj.get(n);
    }

    @Override
    public Collection<Map.Entry<T, T>> edges() {
        return new UndirectedSet(adj.entries());
    }

    @Override
    public Iterator<Map.Entry<T, T>> edgesIter() {
        return new UndirectedIterator(adj.entries().iterator());
    }

    @Override
    public int degree(T n) {
        return adj.get(n).size();
    }

    @Override
    public int size() {
        return adj.size()/2;
    }
    @Override
    public double size(boolean weighted) {
        if (weighted) {
            double total_size = 0.;
            Iterator<Map.Entry<T,T>> it = edgesIter();
            while (it.hasNext()) {
                Map.Entry<T,T> entry = it.next();
                total_size += wts.get(new Tuple(entry.getKey(), entry.getValue()));
            }
            return total_size;
        } else {
            return (double) size();
        }
    }
}
