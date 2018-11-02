/*
 *  ListStatsSite
 * 
 *  Created on Dec 27, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

/**
 *
 * @author Terry Casstevens
 */
public class ListStatsSite extends ListStats {

    private final Stats[] myCache;

    ListStatsSite(GenotypeCallTable genotype) {
        super(genotype, genotype.numberOfSites());
        myCache = new Stats[genotype.numberOfSites()];
    }

    @Override
    public Stats get(int index) {
        if (myCache[index] == null) {
            myCache[index] = myGenotype.siteStats(index);
        }
        return myCache[index];
    }

}
