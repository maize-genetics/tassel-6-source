/*
 *  MaskAlleleProbability
 * 
 *  Created on May 7, 2016
 */
package net.maizegenetics.dna.snp.score;

import java.util.Collection;
import net.maizegenetics.dna.snp.MaskMatrix;
import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class MaskAlleleProbability extends AlleleProbability {

    private final AlleleProbability myProbability;
    private final MaskMatrix myMask;

    public MaskAlleleProbability(AlleleProbability probability, MaskMatrix mask) {
        super(probability.numTaxa(), probability.numSites());
        myProbability = probability;
        myMask = mask;
    }

    @Override
    public float value(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        if (myMask.get(taxon, site)) {
            return 0.0f;
        } else {
            return myProbability.value(taxon, site, scoreType);
        }
    }

    @Override
    Collection<Byte2D> byteStorage() {
        return null;
    }

    public AlleleProbability base() {
        return myProbability;
    }

    public MaskMatrix mask() {
        return myMask;
    }

}
