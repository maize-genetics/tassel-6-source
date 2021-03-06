/*
 *  AlleleProbabilityBuilder
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.MaskMatrix;
import net.maizegenetics.dna.snp.Translate;
import net.maizegenetics.dna.snp.TranslateBuilder;
import net.maizegenetics.dna.snp.byte2d.Byte2D;
import net.maizegenetics.dna.snp.byte2d.Byte2DBuilder;
import net.maizegenetics.taxa.TaxaList;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * @author Terry Casstevens
 */
public class AlleleProbabilityBuilder {

    private final Map<SiteScore.SITE_SCORE_TYPE, Byte2DBuilder> myBuilders = new LinkedHashMap<>();
    private final int myNumSites;

    private AlleleProbabilityBuilder(int numTaxa, int numSites, TaxaList taxaList) {
        for (int i = 0; i < AlleleProbability.ALLELE_PROBABILITY_TYPES.length; i++) {
            myBuilders.put(AlleleProbability.ALLELE_PROBABILITY_TYPES[i], Byte2DBuilder.getInstance(numTaxa, numSites, AlleleProbability.ALLELE_PROBABILITY_TYPES[i], taxaList));
        }
        myNumSites = numSites;
    }

    public static AlleleProbabilityBuilder getAlleleProbabilityInstance(int numTaxa, int numSites, TaxaList taxaList) {
        return new AlleleProbabilityBuilder(numTaxa, numSites, taxaList);
    }

    public static AlleleProbability getFilteredInstance(AlleleProbability base, Translate translate) {
        if (base instanceof FilterAlleleProbability) {
            FilterAlleleProbability filter = (FilterAlleleProbability) base;
            Translate merged = TranslateBuilder.getInstance(filter.myTranslate, translate);
            return new FilterAlleleProbability(filter.myBase, merged);
        }
        return new FilterAlleleProbability(base, translate);
    }

    public static AlleleProbability getMaskInstance(AlleleProbability base, MaskMatrix mask) {
        return new MaskAlleleProbability(base, mask);
    }

    public AlleleProbabilityBuilder addTaxon(int taxon, byte[] values, SiteScore.SITE_SCORE_TYPE type) {
        myBuilders.get(type).addTaxon(taxon, values);
        return this;
    }

    public AlleleProbabilityBuilder addTaxon(int taxon, float[] values, SiteScore.SITE_SCORE_TYPE type) {
        if (myNumSites != values.length) {
            throw new IllegalArgumentException("AlleleProbabilityBuilder: addTaxon: number of values: " + values.length + " doesn't equal number of sites: " + myNumSites);
        }
        byte[] result = SiteScoreUtil.floatToBytePercentage(values);
        myBuilders.get(type).addTaxon(taxon, result);
        return this;
    }

    public AlleleProbability build() {
        Byte2D[] input = new Byte2D[myBuilders.size()];
        int count = 0;
        for (Byte2DBuilder builder : myBuilders.values()) {
            input[count++] = builder.build();
        }
        myBuilders.clear();
        return new AlleleProbability(input);
    }

}
