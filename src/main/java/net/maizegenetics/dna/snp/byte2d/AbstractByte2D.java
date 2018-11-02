/*
 *  AbstractByte2D
 */
package net.maizegenetics.dna.snp.byte2d;

import net.maizegenetics.dna.snp.score.SiteScore;

/**
 *
 * @author Terry Casstevens
 */
public abstract class AbstractByte2D implements Byte2D {

    private final SiteScore.SITE_SCORE_TYPE myScoreType;
    private final int myNumTaxa;
    private final int myNumSites;

    public AbstractByte2D(SiteScore.SITE_SCORE_TYPE scoreType, int numTaxa, int numSites) {
        myScoreType = scoreType;
        myNumTaxa = numTaxa;
        myNumSites = numSites;
    }

    @Override
    public byte[] valuesForAllSites(int taxon) {
        byte[] result = new byte[myNumSites];
        for (int site = 0; site < myNumSites; site++) {
            result[site] = valueForAllele(taxon, site);
        }
        return result;
    }

    @Override
    public byte[] valuesForAllTaxa(int site) {
        byte[] result = new byte[myNumTaxa];
        for (int taxon = 0; taxon < myNumTaxa; taxon++) {
            result[site] = valueForAllele(taxon, site);
        }
        return result;
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
    public SiteScore.SITE_SCORE_TYPE siteScoreType() {
        return myScoreType;
    }

}
