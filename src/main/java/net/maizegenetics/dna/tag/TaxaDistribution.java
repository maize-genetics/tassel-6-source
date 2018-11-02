package net.maizegenetics.dna.tag;

import com.google.common.collect.Multiset;

/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.  Instead this the general implementations use primitive arrays.
 *
 * @author Ed Buckler
 */
public interface TaxaDistribution {

    /**
     * Add the taxa to do the distribution, with an additional depth of 1
     * @param taxaNum
     * @return
     */
    TaxaDistribution increment(int taxaNum);

    /**
     * Distribution of taxa depths across all taxa.  Taxa with zero depth are included in the array
     */
    int[] depths();

    /**
     * Two arrays with the list of taxa with the tag (i.e. depth>0) int[0], and the depth of the taxa in the second array int[1]
     */
    int[][] taxaWithDepths();

    /**
     * Custom run length encoding compression that also use Snappy
     */
    byte[] encodeTaxaDepth();

    /**
     * Multiset version of the taxa depth.  This is a convenient data structure, but it is slow to create compared
     * to depths or taxaWithDepths.  Do not use this if performance is key.
     * @return
     */
    Multiset<Integer> taxaDepthMap();

    /**
     * Total depth across all taxa of the tag
     */
    int totalDepth();

    /**
     * Number of taxa with at least one read
     */
    int numberOfTaxaWithTag();


    /**
     * Number of taxa that depth is being recorded for.  Maximum taxa index = (maxTaxa-1)
     */
    int maxTaxa();

    /**
     * Estimated memory footprint of this taxa distribution.  Can be used to estimate the size of the Map containing
     * these distributions.
     */
    int memorySize();
}
