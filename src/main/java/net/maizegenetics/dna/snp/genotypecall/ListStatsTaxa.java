/*
 *  ListStatsTaxa
 * 
 *  Created on Dec 27, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

/**
 *
 * @author Terry Casstevens
 */
public class ListStatsTaxa extends ListStats {

    private final Stats[] myCache;

    ListStatsTaxa(GenotypeCallTable genotype) {
        super(genotype, genotype.numberOfTaxa());
        myCache = new Stats[genotype.numberOfTaxa()];
    }

    @Override
    public Stats get(int index) {
        if (myCache[index] == null) {
            myCache[index] = myGenotype.taxonStats(index);
        }
        return myCache[index];
    }

}
