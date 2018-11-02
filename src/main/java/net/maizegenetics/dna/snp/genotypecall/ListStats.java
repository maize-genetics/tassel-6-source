/*
 *  ListStats
 * 
 *  Created on Dec 27, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

import java.lang.ref.WeakReference;
import java.util.AbstractList;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Terry Casstevens
 */
public abstract class ListStats extends AbstractList<Stats> {

    private final static Map<GenotypeCallTable, WeakReference<ListStats>> TAXA_INSTANCES = new HashMap<>();
    private final static Map<GenotypeCallTable, WeakReference<ListStats>> SITE_INSTANCES = new HashMap<>();

    protected final GenotypeCallTable myGenotype;
    private final int myNumIndices;

    ListStats(GenotypeCallTable genotype, int numIndices) {
        myGenotype = genotype;
        myNumIndices = numIndices;
    }

    public static ListStats getTaxaInstance(GenotypeCallTable genotype) {

        WeakReference<ListStats> temp = TAXA_INSTANCES.get(genotype);

        ListStats result = null;
        if (temp != null) {
            result = temp.get();
        }

        if (result == null) {
            if (genotype instanceof FilterGenotypeCallTable) {
                FilterGenotypeCallTable filter = (FilterGenotypeCallTable) genotype;
                WeakReference<ListStats> temp1 = TAXA_INSTANCES.get(filter.myBaseGenotype);
                ListStats base = null;
                if (temp1 != null) {
                    base = temp1.get();
                }
                if (base == null || filter.myTranslate.hasSiteTranslations()) {
                    result = new ListStatsTaxa(genotype);
                } else {
                    result = new ListStatsFilterTaxa(filter, base);
                }
            } else {
                result = new ListStatsTaxa(genotype);
            }
            TAXA_INSTANCES.put(genotype, new WeakReference<>(result));
        }

        return result;

    }

    public static ListStats getSiteInstance(GenotypeCallTable genotype) {

        WeakReference<ListStats> temp = SITE_INSTANCES.get(genotype);

        ListStats result = null;
        if (temp != null) {
            result = temp.get();
        }

        if (result == null) {
            if (genotype instanceof FilterGenotypeCallTable) {
                FilterGenotypeCallTable filter = (FilterGenotypeCallTable) genotype;
                WeakReference<ListStats> temp1 = SITE_INSTANCES.get(filter.myBaseGenotype);
                ListStats base = null;
                if (temp1 != null) {
                    base = temp1.get();
                }
                if (base == null || filter.myTranslate.hasTaxaTranslations()) {
                    result = new ListStatsSite(genotype);
                } else {
                    result = new ListStatsFilterSite(filter, base);
                }
            } else {
                result = new ListStatsSite(genotype);
            }
            SITE_INSTANCES.put(genotype, new WeakReference<>(result));
        }

        return result;

    }

    @Override
    public int size() {
        return myNumIndices;
    }

    public byte majorAllele(int index) {
        return get(index).majorAllele();
    }

    public double majorAlleleFrequency(int index) {
        return get(index).majorAlleleFrequency();
    }

    public byte minorAllele(int index) {
        return get(index).minorAllele();
    }

    public double minorAlleleFrequency(int index) {
        return get(index).minorAlleleFrequency();
    }

}
