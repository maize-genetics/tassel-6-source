/*
 *  ListStatsFilterTaxa
 * 
 *  Created on Dec 27, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.Translate;

/**
 *
 * @author Terry Casstevens
 */
public class ListStatsFilterTaxa extends ListStats {

    private final ListStats myBase;
    private final Translate myTranslate;
    private final Stats[] myCache;

    ListStatsFilterTaxa(FilterGenotypeCallTable genotype, ListStats base) {
        super(genotype, genotype.numberOfTaxa());
        myBase = base;
        myTranslate = genotype.myTranslate;
        myCache = new Stats[genotype.numberOfTaxa()];
    }

    @Override
    public Stats get(int index) {
        if (myCache[index] == null) {
            myCache[index] = Stats.getInstance(myBase.get(myTranslate.taxon(index)), index);
        }
        return myCache[index];
    }

}
