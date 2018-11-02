/*
 *  ListStatsFilterSite
 * 
 *  Created on Dec 27, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.Translate;

/**
 *
 * @author Terry Casstevens
 */
public class ListStatsFilterSite extends ListStats {

    private final ListStats myBase;
    private final Translate myTranslate;
    private final Stats[] myCache;

    ListStatsFilterSite(FilterGenotypeCallTable genotype, ListStats base) {
        super(genotype, genotype.numberOfSites());
        myBase = base;
        myTranslate = genotype.myTranslate;
        myCache = new Stats[genotype.numberOfSites()];
    }

    @Override
    public Stats get(int index) {
        if (myCache[index] == null) {
            myCache[index] = Stats.getInstance(myBase.get(myTranslate.site(index)), index);
        }
        return myCache[index];
    }

}
