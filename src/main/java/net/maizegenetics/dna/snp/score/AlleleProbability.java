/*
 *  AlleleProbability
 */
package net.maizegenetics.dna.snp.score;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class AlleleProbability implements SiteScore {

    public static final SiteScore.SITE_SCORE_TYPE[] ALLELE_PROBABILITY_TYPES = new SiteScore.SITE_SCORE_TYPE[]{
        SiteScore.SITE_SCORE_TYPE.ProbA, SiteScore.SITE_SCORE_TYPE.ProbC,
        SiteScore.SITE_SCORE_TYPE.ProbG, SiteScore.SITE_SCORE_TYPE.ProbT,
        SiteScore.SITE_SCORE_TYPE.ProbGap, SiteScore.SITE_SCORE_TYPE.ProbInsertion};

    private final Map<SITE_SCORE_TYPE, Byte2D> myValues;
    private final int myNumTaxa;
    private final int myNumSites;

    AlleleProbability(Byte2D[] values) {
        if (values.length == 0) {
            throw new IllegalArgumentException("AlleleProbability: init: no values provided.");
        }
        myValues = new HashMap<>();
        myNumTaxa = values[0].numTaxa();
        myNumSites = values[0].numSites();
        for (int i = 0; i < values.length; i++) {
            if ((myNumTaxa != values[i].numTaxa()) || (myNumSites != values[i].numSites())) {
                throw new IllegalArgumentException("AlleleProbability: init: number of taxa or sites don't match for all values.");
            }
            myValues.put(values[i].siteScoreType(), values[i]);
        }
    }

    public AlleleProbability(int numTaxa, int numSites) {
        myNumTaxa = numTaxa;
        myNumSites = numSites;
        myValues = null;
    }

    public float value(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        return SiteScoreUtil.byteToFloatPercentage(myValues.get(scoreType).valueForAllele(taxon, site));
    }

    Collection<Byte2D> byteStorage() {
        return myValues.values();
    }

    @Override
    public Set<SITE_SCORE_TYPE> siteScoreTypes() {
        return new HashSet<>(Arrays.asList(ALLELE_PROBABILITY_TYPES));
    }

    @Override
    public int numTaxa() {
        return myNumTaxa;
    }

    @Override
    public int numSites() {
        return myNumSites;
    }
}
