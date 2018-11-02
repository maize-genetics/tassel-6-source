package net.maizegenetics.taxa;

import com.google.common.collect.ImmutableMap;

import org.apache.log4j.Logger;

import java.util.*;

/**
 * In memory immutable instance of {@link TaxaList}. Basic list of taxa
 * (samples) that are used in Alignments and other purposes.
 *
 * Use {@link TaxaListBuilder} to instantiate.
 *
 * @author Ed Buckler
 *
 */
class TaxaArrayList implements TaxaList {

    private static final Logger myLogger = Logger.getLogger(TaxaArrayList.class);
    private final List<Taxon> myTaxaList;
    private final int myNumTaxa;
    private final ImmutableMap<String, Integer> myNameToIndex;

    TaxaArrayList(TaxaListBuilder builder) {
        List<Taxon> srcList = builder.getImmutableList();
        myTaxaList = new ArrayList<Taxon>(srcList.size());
        myNumTaxa = srcList.size();
        int index = 0;
        ImmutableMap.Builder<String, Integer> nToIBuilder=new ImmutableMap.Builder<>();
        for (Taxon Taxon : srcList) {
            myTaxaList.add(Taxon);
            nToIBuilder.put(Taxon.getName(), index);
            index++;
        }
        myNameToIndex=nToIBuilder.build();
    }

    @Override
    public int numberOfTaxa() {
        return myNumTaxa;
    }

    @Override
    public String taxaName(int index) {
        return myTaxaList.get(index).getName();
    }

    @Override
    public int size() {
        return myNumTaxa;
    }

    @Override
    public int indexOf(String name) {
        Integer index=myNameToIndex.get(name);
        if(index==null) return -1;
        return index;
    }

    @Override
    public int indexOf(Taxon taxon) {
        return indexOf(taxon.getName());
    }

    @Override
    public boolean isEmpty() {
        return myTaxaList.isEmpty();
    }

    @Override
    public boolean contains(Object o) {
        return myTaxaList.contains(o);
    }

    @Override
    public Iterator<Taxon> iterator() {
        return myTaxaList.iterator();
    }

    @Override
    public Object[] toArray() {
        return myTaxaList.toArray();
    }

    @Override
    public <T> T[] toArray(T[] a) {
        return myTaxaList.toArray(a);
    }

    @Override
    public boolean add(Taxon Taxon) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        return myTaxaList.containsAll(c);
    }

    @Override
    public boolean addAll(Collection<? extends Taxon> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean addAll(int index, Collection<? extends Taxon> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public Taxon get(int index) {
        return myTaxaList.get(index);
    }

    @Override
    public Taxon set(int index, Taxon element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public void add(int index, Taxon element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public Taxon remove(int index) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public int indexOf(Object o) {
        return indexOf((Taxon) o);
    }

    @Override
    public int lastIndexOf(Object o) {
        return myTaxaList.lastIndexOf(o);
    }

    @Override
    public ListIterator<Taxon> listIterator() {
        return listIterator(0);
    }

    @Override
    public ListIterator<Taxon> listIterator(final int index) {
        return new ListIterator<Taxon>() {
            private final ListIterator<Taxon> i = myTaxaList.listIterator(index);

            public boolean hasNext() {
                return i.hasNext();
            }

            public Taxon next() {
                return i.next();
            }

            public boolean hasPrevious() {
                return i.hasPrevious();
            }

            public Taxon previous() {
                return i.previous();
            }

            public int nextIndex() {
                return i.nextIndex();
            }

            public int previousIndex() {
                return i.previousIndex();
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }

            public void set(Taxon e) {
                throw new UnsupportedOperationException();
            }

            public void add(Taxon e) {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Override
    public List<Taxon> subList(int fromIndex, int toIndex) {
        return Collections.unmodifiableList(myTaxaList.subList(fromIndex, toIndex));
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof TaxaList)) return false;
        TaxaList otherTaxa = (TaxaList) obj;
        if (numberOfTaxa() != otherTaxa.numberOfTaxa()) return false;
        Iterator<Taxon> myIter = myTaxaList.iterator();
        Iterator<Taxon> otherIter = otherTaxa.iterator();
        while (myIter.hasNext()) {
            if ( !(myIter.next().equals(otherIter.next())) ) {
                return false;
            }
        }
        
        return true;
    }
}
