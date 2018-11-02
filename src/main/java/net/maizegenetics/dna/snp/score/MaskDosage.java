/*
 *  MaskDosage
 * 
 *  Created on Dec 16, 2016
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.MaskMatrix;

/**
 *
 * @author Terry Casstevens
 */
public class MaskDosage extends Dosage {

    private final Dosage myBase;
    private final MaskMatrix myMask;

    MaskDosage(Dosage dosage, MaskMatrix mask) {
        super(mask.numTaxa(), mask.numSites());
        myBase = dosage;
        myMask = mask;
    }

    @Override
    public byte value(int taxon, int site) {
        if (myMask.get(taxon, site)) {
            return 0;
        }
        return myBase.value(taxon, site);
    }

}
