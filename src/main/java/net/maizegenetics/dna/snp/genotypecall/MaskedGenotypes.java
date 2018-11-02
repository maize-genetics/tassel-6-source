/*
 *  MaskedGenotypes
 * 
 *  Created on May 8, 2015
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskedGenotypes {

    final private int myNumTaxa;
    final private int myNumSites;
    final private BitSet[] myMasks;

    public MaskedGenotypes(int numTaxa, int numSites) {
        myNumTaxa = numTaxa;
        myNumSites = numSites;
        myMasks = new BitSet[myNumTaxa];
        for (int t = 0; t < myNumTaxa; t++) {
            myMasks[t] = new OpenBitSet(myNumSites);
        }
    }

    public void set(int taxon, int site) {
        myMasks[taxon].fastSet(site);
    }

    public boolean get(int taxon, int site) {
        return myMasks[taxon].fastGet(site);
    }

}
