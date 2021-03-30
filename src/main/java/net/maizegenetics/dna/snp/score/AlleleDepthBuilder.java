/*
 *  AlleleDepthBuilder
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.MaskMatrix;
import net.maizegenetics.dna.snp.Translate;
import net.maizegenetics.dna.snp.TranslateBuilder;
import net.maizegenetics.dna.snp.byte2d.Byte2D;
import net.maizegenetics.dna.snp.byte2d.Byte2DBuilder;
import net.maizegenetics.taxa.TaxaList;

/**
 * @author Terry Casstevens
 */
public class AlleleDepthBuilder {

    private Byte2DBuilder[] myBuilders;
    private final int myNumSites;

    private AlleleDepthBuilder(int numTaxa, int numSites, TaxaList taxaList) {
        myBuilders = new Byte2DBuilder[AlleleDepth.NUM_ALLELE_DEPTH_TYPES];
        for (int i = 0; i < AlleleDepth.NUM_ALLELE_DEPTH_TYPES; i++) {
            myBuilders[i] = Byte2DBuilder.getInstance(numTaxa, numSites, AlleleDepth.ALLELE_DEPTH_TYPES[i], taxaList);
        }
        myNumSites = numSites;
    }

    /**
     * This returns an AlleleDepthBuilder where depths are stored in memory.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites
     * @param taxaList taxa list
     *
     * @return AlleleDepthBuilder
     */
    public static AlleleDepthBuilder getInstance(int numTaxa, int numSites, TaxaList taxaList) {
        return new AlleleDepthBuilder(numTaxa, numSites, taxaList);
    }

    /**
     * This creates a filtered AlleleDepth.
     *
     * @param base original AlleleDepth
     * @param translate translate
     *
     * @return filtered AlleleDepth
     */
    public static AlleleDepth getFilteredInstance(AlleleDepth base, Translate translate) {
        if (base instanceof FilterAlleleDepth) {
            FilterAlleleDepth filter = (FilterAlleleDepth) base;
            Translate merged = TranslateBuilder.getInstance(filter.myTranslate, translate);
            return new FilterAlleleDepth(filter.myBase, merged);
        }
        return new FilterAlleleDepth(base, translate);
    }

    public static AlleleDepth getMaskInstance(AlleleDepth base, MaskMatrix mask) {
        return new MaskAlleleDepth(base, mask);
    }

    public AlleleDepthBuilder addTaxon(int taxon, int[] values, SiteScore.SITE_SCORE_TYPE type) {
        if (myNumSites != values.length) {
            throw new IllegalArgumentException("AlleleDepthBuilder: addTaxon: number of values: " + values.length + " doesn't equal number of sites: " + myNumSites);
        }
        byte[] result = AlleleDepthUtil.depthIntToByte(values);
        myBuilders[type.getIndex()].addTaxon(taxon, result);
        return this;
    }

    public AlleleDepthBuilder addTaxon(int taxon, byte[][] values) {
        if (AlleleDepth.NUM_ALLELE_DEPTH_TYPES != values.length) {
            throw new IllegalArgumentException("AlleleDepthBuilder: addTaxon: number of alleles: " + values.length + " doesn't equals: " + AlleleDepth.NUM_ALLELE_DEPTH_TYPES);
        }
        if (myNumSites != values[0].length) {
            throw new IllegalArgumentException("AlleleDepthBuilder: addTaxon: number of values: " + values[0].length + " doesn't equal number of sites: " + myNumSites);
        }
        for (int i = 0; i < AlleleDepth.NUM_ALLELE_DEPTH_TYPES; i++) {
            myBuilders[i].addTaxon(taxon, values[i]);
        }
        return this;
    }

    /**
     * Set depth for range of sites and alleles for a taxon simultaneously.
     * First dimension of depths is number of alleles (6 for Nucleotide) and
     * second dimension is sites.
     *
     * @param taxon Index of taxon
     * @param siteOffset site offset
     * @param depths array[allele][site] of all values
     *
     * @return builder
     */
    public AlleleDepthBuilder setDepthRangeForTaxon(int taxon, int siteOffset, byte[][] depths) {

        int numAlleles = depths.length;
        if (numAlleles != AlleleDepth.NUM_ALLELE_DEPTH_TYPES) {
            throw new IllegalArgumentException("AlleleDepthBuilder: setDepthRangeForTaxon: value number of alleles: " + numAlleles + " should be: " + AlleleDepth.NUM_ALLELE_DEPTH_TYPES);
        }
        for (int a = 0; a < AlleleDepth.NUM_ALLELE_DEPTH_TYPES; a++) {
            myBuilders[a].setDepthRangeForTaxon(taxon, siteOffset, depths[a]);
        }

        return this;

    }

    public void reorderPositions(int[] newIndices) {
        for (int i = 0; i < AlleleDepth.NUM_ALLELE_DEPTH_TYPES; i++) {
            myBuilders[i].reorderPositions(newIndices);
        }
    }

    public AlleleDepth build() {
        Byte2D[] input = new Byte2D[AlleleDepth.NUM_ALLELE_DEPTH_TYPES];
        for (int a = 0; a < AlleleDepth.NUM_ALLELE_DEPTH_TYPES; a++) {
            input[a] = myBuilders[a].build();
        }
        myBuilders = null;
        return new AlleleDepth(input);
    }

}
