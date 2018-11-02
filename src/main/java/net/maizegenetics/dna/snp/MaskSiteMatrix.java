/*
 *  MaskSiteMatrix
 * 
 *  Created on May 6, 2016
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.UnmodifiableBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskSiteMatrix implements MaskMatrix {

    private final BitSet[] myBitSets;
    private final int myNumTaxa;
    private final int myNumSites;

    MaskSiteMatrix(BitSet[] bitSets, int numTaxa, int numSites) {
        myBitSets = bitSets;
        myNumTaxa = numTaxa;
        myNumSites = numSites;
    }

    @Override
    public boolean get(int taxon, int site) {
        return myBitSets[site].fastGet(taxon);
    }

    /**
     * Returns false if specified taxon is definitely not masked. Otherwise
     * returns true. A true returned doesn't necessarily mean this taxon is
     * masked. This can be used to optimize performance of masked genotype table
     * components. If false returned, checking for masking can be skipped.
     *
     * @param taxon taxon
     *
     * @return false if taxon definitely not masked, true otherwise
     */
    @Override
    public boolean isTaxonMaskedHint(int taxon) {
        for (int s = 0; s < myNumSites; s++) {
            if (myBitSets[s].fastGet(taxon)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Mask for specified taxon
     *
     * @param taxon taxon
     *
     * @return mask
     */
    @Override
    public BitSet maskForTaxon(int taxon) {
        BitSet result = new OpenBitSet(myNumSites);
        for (int s = 0; s < myNumSites; s++) {
            if (myBitSets[s].fastGet(taxon)) {
                result.fastSet(s);
            }
        }
        return result;
    }

    /**
     * Returns false if specified site is definitely not masked. Otherwise
     * returns true. A true returned doesn't necessarily mean this site is
     * masked. This can be used to optimize performance of masked genotype table
     * components. If false returned, checking for masking can be skipped.
     *
     * @param site site
     *
     * @return false if site definitely not masked, true otherwise
     */
    @Override
    public boolean isSiteMaskedHint(int site) {
        return myBitSets[site].cardinality() != 0;
    }

    /**
     * Mask for specified site
     *
     * @param site site
     *
     * @return mask
     */
    @Override
    public BitSet maskForSite(int site) {
        return UnmodifiableBitSet.getInstance(myBitSets[site]);
    }

    @Override
    public int numTaxa() {
        return myNumTaxa;
    }

    @Override
    public int numSites() {
        return myNumSites;
    }

    @Override
    public boolean isSiteOptimized() {
        return true;
    }

}
