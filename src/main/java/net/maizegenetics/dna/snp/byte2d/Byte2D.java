/*
 *  Byte2D
 */
package net.maizegenetics.dna.snp.byte2d;

import net.maizegenetics.dna.snp.score.SiteScore;

/**
 * @author Terry Casstevens
 */
public interface Byte2D {

    public byte valueForAllele(int taxon, int site);

    public byte[] valuesForAllSites(int taxon);
    
    public byte[] valuesForAllTaxa(int site);

    public int numTaxa();

    public int numSites();

    public SiteScore.SITE_SCORE_TYPE siteScoreType();

}
