/*
 * BitStorage
 */
package net.maizegenetics.dna.snp.bit;

import net.maizegenetics.util.BitSet;

/**
 * Interface provides genotypes in a binary fashion to be used in rapid
 * computation. See the package descriptions
 * {@link net.maizegenetics.dna.snp.bit} for more information on how bits
 * encoding is used throughout TASSEL.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public interface BitStorage {

    /**
     * Returns sequence of true/false values indicating whether taxon at each
     * site matches a specific allele.
     *
     * @param taxon taxon
     *
     * @return sequence of true/false values.
     */
    public BitSet allelePresenceForAllSites(int taxon);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon matches a specific allele.
     *
     * @param site site
     *
     * @return sequence of true/false values.
     */
    public BitSet allelePresenceForAllTaxa(int site);

    /**
     * Returns sequence of true/false values indicating whether taxon at sites
     * (in given blocks, 64 sites per block including start block but excluding
     * end block) matches a specific allele.
     *
     * @param taxon taxon
     * @param startBlock starting block
     * @param endBlock end block
     *
     * @return sequence of true/false values.
     */
    public long[] allelePresenceForSitesBlock(int taxon, int startBlock, int endBlock);

    /**
     * Returns sequence of true/false values indicating whether taxon at each
     * site for given parent matches a specific allele (based on frequency).
     * Allele number of value 0 would be the major allele. Allele number of
     * value 1 would be the minor allele. Allele number of value 2 would be the
     * third most frequent allele value and so on.
     *
     * @param taxon taxon
     * @param firstParent true for first parent (false for second parent)
     *
     * @return sequence of true/false values.
     */
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon for given parent matches a specific allele (based on frequency).
     * Allele number of value 0 would be the major allele. Allele number of
     * value 1 would be the minor allele. Allele number of value 2 would be the
     * third most frequent allele value and so on.
     *
     * @param site site
     * @param firstParent true for first parent (false for second parent)
     *
     * @return sequence of true/false values.
     */
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent);

    /**
     * Returns sequence of true/false values indicating whether taxon at sites
     * (in given blocks, 64 sites per block including start block but excluding
     * end block) for given parent matches a specific allele (based on
     * frequency). Allele number of value 0 would be the major allele. Allele
     * number of value 1 would be the minor allele. Allele number of value 2
     * would be the third most frequent allele value and so on.
     *
     * @param taxon taxon
     * @param firstParent true for first parent (false for second parent)
     * @param startBlock starting block
     * @param endBlock end block
     *
     * @return sequence of true/false values.
     */
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, int startBlock, int endBlock);

}
