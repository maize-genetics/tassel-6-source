/*
 *  MaskAlleleDepth
 * 
 *  Created on May 6, 2016
 */
package net.maizegenetics.dna.snp.score;

import java.util.Collection;
import net.maizegenetics.dna.snp.MaskMatrix;
import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class MaskAlleleDepth extends AlleleDepth {

    private final AlleleDepth myDepth;
    private final MaskMatrix myMask;

    public MaskAlleleDepth(AlleleDepth depth, MaskMatrix mask) {
        super(depth.numTaxa(), depth.numSites());
        if (depth.numTaxa() != mask.numTaxa() || depth.numSites() != mask.numSites()) {
            throw new IllegalArgumentException("MaskAlleleDepth: init: depth and mask dimensions don't match");
        }
        myDepth = depth;
        myMask = mask;
    }

    @Override
    public int value(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        if (myMask.get(taxon, site)) {
            return 0;
        } else {
            return myDepth.value(taxon, site, scoreType);
        }
    }

    @Override
    public byte valueByte(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        if (myMask.get(taxon, site)) {
            return 0;
        } else {
            return myDepth.valueByte(taxon, site, scoreType);
        }
    }

    @Override
    Collection<Byte2D> byteStorage() {
        return null;
    }

    public AlleleDepth base() {
        return myDepth;
    }

    public MaskMatrix mask() {
        return myMask;
    }

}
