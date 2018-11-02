/*
 *  ReferenceProbabilityBuilder
 */
package net.maizegenetics.dna.snp.score;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.maizegenetics.dna.snp.MaskMatrix;

import net.maizegenetics.dna.snp.Translate;
import net.maizegenetics.dna.snp.TranslateBuilder;
import net.maizegenetics.dna.snp.byte2d.Byte2D;
import net.maizegenetics.dna.snp.byte2d.Byte2DBuilder;
import net.maizegenetics.taxa.TaxaList;

/**
 *
 * @author Terry Casstevens
 */
public class ReferenceProbabilityBuilder {

    private Byte2DBuilder myBuilder;
    private final int myNumSites;

    private ReferenceProbabilityBuilder(int numTaxa, int numSites, TaxaList taxaList) {
        myBuilder = Byte2DBuilder.getInstance(numTaxa, numSites, SiteScore.SITE_SCORE_TYPE.ReferenceProbablity, taxaList);
        myNumSites = numSites;
    }

    private ReferenceProbabilityBuilder(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        myBuilder = Byte2DBuilder.getInstance(writer, numSites, SiteScore.SITE_SCORE_TYPE.ReferenceProbablity, taxaList);
        myNumSites = numSites;
    }

    public static ReferenceProbabilityBuilder getInstance(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        return new ReferenceProbabilityBuilder(writer, numTaxa, numSites, taxaList);
    }

    public static ReferenceProbabilityBuilder getInstance(int numTaxa, int numSites, TaxaList taxaList) {
        return new ReferenceProbabilityBuilder(numTaxa, numSites, taxaList);
    }

    public static ReferenceProbability getFilteredInstance(ReferenceProbability base, Translate translate) {
        if (base instanceof FilterReferenceProbability) {
            FilterReferenceProbability filter = (FilterReferenceProbability) base;
            Translate merged = TranslateBuilder.getInstance(filter.myTranslate, translate);
            return new FilterReferenceProbability(filter.myBase, merged);
        }
        return new FilterReferenceProbability(base, translate);
    }

    public static ReferenceProbability getMaskInstance(ReferenceProbability base, MaskMatrix mask) {
        return new MaskReferenceProbability(base, mask);
    }

    public static ReferenceProbability getInstance(IHDF5Reader reader) {
        return new ReferenceProbability(Byte2DBuilder.getInstance(reader, SiteScore.SITE_SCORE_TYPE.ReferenceProbablity));
    }

    public ReferenceProbabilityBuilder addTaxon(int taxon, float[] values) {
        if (myNumSites != values.length) {
            throw new IllegalArgumentException("ImputeProbabilityBuilder: addTaxon: number of values: " + values.length + " doesn't equal number of sites: " + myNumSites);
        }
        try {
            myBuilder.addTaxon(taxon, SiteScoreUtil.floatToBytePercentage(values));
        } catch (Exception e) {
            throw new IllegalArgumentException("ImputeProbabilityBuilder: addTaxon: taxon number: " + taxon + ". " + e.getMessage());
        }
        return this;
    }

    public ReferenceProbability build() {
        Byte2D input = myBuilder.build();
        myBuilder = null;
        return new ReferenceProbability(input);
    }

}
