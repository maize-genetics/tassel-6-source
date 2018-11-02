/*
 *  FilterList
 * 
 *  Created on Aug 5, 2015
 */
package net.maizegenetics.dna.snp;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

/**
 *
 * @author Terry Casstevens
 */
public class FilterList implements List<Filter> {
    
    private final List<Filter> myFilters;
    
    public FilterList(List<Filter> filters) {
        myFilters = new ArrayList<>(filters);
    }
    
    public FilterList(Filter filter) {
        myFilters = new ArrayList<>();
        myFilters.add(filter);
    }
    
    @Override
    public int size() {
        return myFilters.size();
    }
    
    @Override
    public boolean isEmpty() {
        return myFilters.isEmpty();
    }
    
    @Override
    public boolean contains(Object o) {
        return myFilters.contains(o);
    }
    
    @Override
    public Iterator<Filter> iterator() {
        return myFilters.iterator();
    }
    
    @Override
    public Object[] toArray() {
        return myFilters.toArray();
    }
    
    @Override
    public <T> T[] toArray(T[] a) {
        return myFilters.toArray(a);
    }
    
    @Override
    public boolean add(Filter e) {
        throw new UnsupportedOperationException("Not supported. Immutatable.");
    }
    
    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("Not supported. Immutatable.");
    }
    
    @Override
    public boolean containsAll(Collection<?> c) {
        return myFilters.containsAll(c);
    }
    
    @Override
    public boolean addAll(Collection<? extends Filter> c) {
        throw new UnsupportedOperationException("Not supported. Immutatable.");
    }
    
    @Override
    public boolean addAll(int index, Collection<? extends Filter> c) {
        throw new UnsupportedOperationException("Not supported. Immutatable.");
    }
    
    @Override
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException("Not supported. Immutatable.");
    }
    
    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("Not supported. Immutatable.");
    }
    
    @Override
    public void clear() {
        throw new UnsupportedOperationException("Not supported. Immutatable.");
    }
    
    @Override
    public Filter get(int index) {
        return myFilters.get(index);
    }
    
    @Override
    public FilterSite set(int index, Filter element) {
        throw new UnsupportedOperationException("Not supported. Immutatable.");
    }
    
    @Override
    public void add(int index, Filter element) {
        throw new UnsupportedOperationException("Not supported. Immutatable.");
    }
    
    @Override
    public FilterSite remove(int index) {
        throw new UnsupportedOperationException("Not supported. Immutatable.");
    }
    
    @Override
    public int indexOf(Object o) {
        return myFilters.indexOf(o);
    }
    
    @Override
    public int lastIndexOf(Object o) {
        return myFilters.lastIndexOf(o);
    }
    
    @Override
    public ListIterator<Filter> listIterator() {
        return myFilters.listIterator();
    }
    
    @Override
    public ListIterator<Filter> listIterator(int index) {
        return myFilters.listIterator(index);
    }
    
    @Override
    public List<Filter> subList(int fromIndex, int toIndex) {
        return myFilters.subList(fromIndex, toIndex);
    }

}
