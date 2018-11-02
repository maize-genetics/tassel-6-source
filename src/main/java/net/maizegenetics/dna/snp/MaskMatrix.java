/*
 *  MaskMatrix
 * 
 *  Created on May 6, 2016
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.util.BitSet;

/**
 *
 * @author Terry Casstevens
 */
public interface MaskMatrix {

    public boolean get(int taxon, int site);

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
    public boolean isTaxonMaskedHint(int taxon);

    /**
     * Mask for specified taxon
     *
     * @param taxon taxon
     *
     * @return mask
     */
    public BitSet maskForTaxon(int taxon);

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
    public boolean isSiteMaskedHint(int site);

    /**
     * Mask for specified site
     *
     * @param site site
     *
     * @return mask
     */
    public BitSet maskForSite(int site);

    public int numTaxa();

    public int numSites();

    public boolean isSiteOptimized();

}
