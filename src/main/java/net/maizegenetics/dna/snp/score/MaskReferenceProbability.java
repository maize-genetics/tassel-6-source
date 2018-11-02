/*
 *  MaskReferenceProbability
 * 
 *  Created on Dec 16, 2016
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.MaskMatrix;

/**
 *
 * @author Terry Casstevens
 */
public class MaskReferenceProbability extends ReferenceProbability {

    private final ReferenceProbability myBase;
    private final MaskMatrix myMask;

    MaskReferenceProbability(ReferenceProbability referenceProbability, MaskMatrix mask) {
        super(mask.numTaxa(), mask.numSites());
        myBase = referenceProbability;
        myMask = mask;
    }

    @Override
    public float value(int taxon, int site) {
        if (myMask.get(taxon, site)) {
            return 0.0f;
        }
        return myBase.value(taxon, site);
    }

}
