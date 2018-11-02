/*
 *  SiteScore
 */
package net.maizegenetics.dna.snp.score;

import java.util.Set;

/**
 *
 * @author Terry Casstevens
 */
public interface SiteScore {

    public static enum SITE_SCORE_TYPE {

        None(0), QualityScore(0), ReferenceProbablity(0), Dosage(0),
        DepthA(0), DepthC(1), DepthG(2), DepthT(3), DepthGap(4), DepthInsertion(5),
        ProbA(0), ProbC(1), ProbG(2), ProbT(3), ProbGap(4), ProbInsertion(5);

        private final int myIndex;

        SITE_SCORE_TYPE(int index) {
            myIndex = index;
        }

        public int getIndex() {
            return myIndex;
        }

    };

    /**
     * Return the site scores types.
     *
     * @return site score types.
     */
    public Set<SITE_SCORE_TYPE> siteScoreTypes();

    public int numTaxa();

    public int numSites();

}
