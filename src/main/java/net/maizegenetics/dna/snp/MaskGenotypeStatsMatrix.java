/*
 *  MaskGenotypeStatsMatrix
 * 
 *  Created on Jan 6, 2017
 */
package net.maizegenetics.dna.snp;

import java.util.function.BiPredicate;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.ListStats;
import net.maizegenetics.dna.snp.genotypecall.Stats;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskGenotypeStatsMatrix extends AbstractMaskMatrix {

    private final GenotypeCallTable myGenotype;
    private final ListStats myStats;
    private ListStats myTaxaStats = null;
    private final BiPredicate<Byte, Stats> myPredicate;

    MaskGenotypeStatsMatrix(GenotypeCallTable genotype, BiPredicate<Byte, Stats> predicate) {
        super(genotype.numberOfTaxa(), genotype.numberOfSites());
        myGenotype = genotype;
        myPredicate = predicate;
        myStats = ListStats.getSiteInstance(genotype);
    }

    @Override
    protected BitSet siteMask(int site) {
        BitSet result = new OpenBitSet(myNumTaxa);
        byte[] temp = myGenotype.genotypeForAllTaxa(site);
        for (int t = 0; t < myNumTaxa; t++) {
            if (myPredicate.test(temp[t], myStats.get(site))) {
                result.fastSet(t);
            }
        }
        return result;
    }

    @Override
    protected BitSet taxonMask(int taxon) {
        if (myTaxaStats == null) {
            myTaxaStats = ListStats.getTaxaInstance(myGenotype);
        }
        BitSet result = new OpenBitSet(myNumSites);
        byte[] temp = myGenotype.genotypeForAllSites(taxon);
        for (int t = 0; t < myNumTaxa; t++) {
            if (myPredicate.test(temp[t], myTaxaStats.get(taxon))) {
                result.fastSet(t);
            }
        }
        return result;
    }

    @Override
    protected boolean isMasked(int taxon, int site) {
        return myPredicate.test(myGenotype.genotype(taxon, site), myStats.get(site));
    }

}
