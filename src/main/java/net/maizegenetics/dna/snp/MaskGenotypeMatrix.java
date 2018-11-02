/*
 *  MaskGenotypeMatrix
 * 
 *  Created on Jan 8, 2017
 */
package net.maizegenetics.dna.snp;

import java.util.function.Predicate;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskGenotypeMatrix extends AbstractMaskMatrix {

    private final GenotypeCallTable myGenotype;
    private final Predicate<Byte> myPredicate;

    MaskGenotypeMatrix(GenotypeCallTable genotype, Predicate<Byte> predicate) {
        super(genotype.numberOfTaxa(), genotype.numberOfSites());
        myGenotype = genotype;
        myPredicate = predicate;
    }

    @Override
    protected BitSet siteMask(int site) {
        BitSet result = new OpenBitSet(myNumTaxa);
        byte[] temp = myGenotype.genotypeForAllTaxa(site);
        for (int t = 0; t < myNumTaxa; t++) {
            if (myPredicate.test(temp[t])) {
                result.fastSet(t);
            }
        }
        return result;
    }

    @Override
    protected BitSet taxonMask(int taxon) {
        BitSet result = new OpenBitSet(myNumSites);
        byte[] temp = myGenotype.genotypeForAllSites(taxon);
        for (int t = 0; t < myNumTaxa; t++) {
            if (myPredicate.test(temp[t])) {
                result.fastSet(t);
            }
        }
        return result;
    }

    @Override
    protected boolean isMasked(int taxon, int site) {
        return myPredicate.test(myGenotype.genotype(taxon, site));
    }

}
