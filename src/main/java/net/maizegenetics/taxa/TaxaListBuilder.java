package net.maizegenetics.taxa;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;

import java.util.*;

/**
 * A builder for creating immutable {@link TaxaList} instances.
 * Example:<pre>   {@code
 *   TaxaListBuilder tlb=new TaxaListBuilder();
 *   for (int i = 0; i < 10; i++) {
 *      Taxon at= new Taxon.Builder("Z"+i+":Line:mays:Zea")
 *           .inbreedF(0.99f)
 *           .parents("B73","B97")
 *           .pedigree("(B73xB97)S6I1")
 *           .build();
 *       tlb.add(at);
 *       }
 *   TaxaList tl=tlb.build();}</pre>
 * <p></p>
 *
 * @author Ed Buckler
 */
public class TaxaListBuilder {
    //TODO need to move union and intersection utils to the builder
    private List<Taxon> myTaxaList;
    private final HashMap<Taxon, Integer> tempLookup;

    public TaxaListBuilder() {
        myTaxaList = new ArrayList<>();
        tempLookup = new HashMap<>();
    }

    public static TaxaList getInstance(int numTaxa) {
        TaxaListBuilder builder = new TaxaListBuilder();
        for (int i = 0; i < numTaxa; i++) {
            builder.add(new Taxon("Taxa_" + i));
        }
        return builder.build();
    }

    public TaxaListBuilder add(Taxon taxon) {
        if (tempLookup.containsKey(taxon)) {
            throw new IllegalStateException("Taxon [" + taxon.getName() + "] already exists in the list.  Duplicated taxa not allowed.");
        }
        myTaxaList.add(taxon);
        tempLookup.put(taxon, myTaxaList.size() - 1);
        return this;
    }

    public TaxaListBuilder add(String taxon) {
        Taxon obj = new Taxon(taxon);
        return add(obj);
    }

    public boolean contains(Taxon taxon) {
        return tempLookup.containsKey(taxon);
    }

    public TaxaListBuilder addOrMerge(Taxon taxon) {
        if (tempLookup.containsKey(taxon)) {
            int indexOfOrig = tempLookup.get(taxon);
            Taxon orgTaxon = myTaxaList.get(indexOfOrig);
            Taxon.Builder tb = new Taxon.Builder(orgTaxon);
            for (Map.Entry<String, String> entry : taxon.getAnnotation().getAllAnnotationEntries()) {
                if (!orgTaxon.getAnnotation().isAnnotatedWithValue(entry.getKey(), entry.getValue()))
                    tb.addAnno(entry.getKey(), entry.getValue());
            }
            taxon = tb.build();
            myTaxaList.set(indexOfOrig, taxon);
            tempLookup.put(taxon, indexOfOrig);
        } else {
            myTaxaList.add(taxon);
            tempLookup.put(taxon, myTaxaList.size() - 1);
        }
        return this;
    }

    public TaxaListBuilder addAll(Collection<Taxon> taxa) {
        for (Taxon taxon : taxa) {
            add(taxon);
        }
        return this;
    }

    public TaxaListBuilder addAll(GenotypeTable a) {
        for (Taxon taxon : a.taxa()) {
            add(taxon);
        }
        return this;
    }

    public TaxaListBuilder addAll(String[] taxa) {
        for (int i = 0, n = taxa.length; i < n; i++) {
            add(taxa[i]);
        }
        return this;
    }

    public TaxaListBuilder addAll(List<String> taxa) {
        for (int i = 0, n = taxa.size(); i < n; i++) {
            add(new Taxon.Builder(taxa.get(i)).build());
        }
        return this;
    }

    public TaxaListBuilder addAll(Taxon[] taxa) {
        for (int i = 0, n = taxa.length; i < n; i++) {
            add(new Taxon.Builder(taxa[i]).build());
        }
        return this;
    }

    public TaxaListBuilder addAll(TaxaListBuilder builder) {
        for (Taxon current : builder.myTaxaList) {
            add(current);
        }
        return this;
    }

    public TaxaList build() {
        return new TaxaArrayList(this);
    }

    //Default package private method to hand the list to the instance
    List<Taxon> getImmutableList() {
        return Collections.unmodifiableList(myTaxaList);
    }

    public TaxaListBuilder sortTaxaAlphabetically(GenotypeCallTableBuilder genotypes) {
        int numTaxa = myTaxaList.size();
        if (numTaxa != genotypes.getTaxaCount()) {
            throw new IllegalArgumentException("TaxaListBuilder: sortTaxaAlphabetically: taxa list size: " + numTaxa + " doesn't match genotypes num taxa: " + genotypes.getTaxaCount());
        }
        genotypes.reorderTaxa(sortAlphabetically());
        return this;
    }

    public TaxaListBuilder sortTaxaAlphabetically() {
        sortAlphabetically();
        return this;
    }

    public int[] sortAlphabetically() {

        int numTaxa = myTaxaList.size();

        final int indicesOfSortByTaxa[] = new int[numTaxa];
        for (int i = 0; i < indicesOfSortByTaxa.length; i++) {
            indicesOfSortByTaxa[i] = i;
        }

        Swapper swapTaxa = new Swapper() {
            @Override
            public void swap(int a, int b) {
                int temp = indicesOfSortByTaxa[a];
                indicesOfSortByTaxa[a] = indicesOfSortByTaxa[b];
                indicesOfSortByTaxa[b] = temp;
            }
        };

        IntComparator compTaxa = new IntComparator() {
            @Override
            public int compare(int a, int b) {
                return myTaxaList.get(indicesOfSortByTaxa[a]).compareTo(myTaxaList.get(indicesOfSortByTaxa[b]));
            }
        };

        GenericSorting.quickSort(0, indicesOfSortByTaxa.length, compTaxa, swapTaxa);

        List<Taxon> temp = new ArrayList<>(numTaxa);
        for (int t = 0; t < numTaxa; t++) {
            temp.add(myTaxaList.get(indicesOfSortByTaxa[t]));
        }

        myTaxaList = temp;

        return indicesOfSortByTaxa;

    }

    public int numberOfTaxa() {
        return myTaxaList.size();
    }
}
