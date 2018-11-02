/*
 *  FilterReferenceProbability
 * 
 *  Created on Dec 16, 2016
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.Translate;

/**
 *
 * @author Terry Casstevens
 */
public class FilterReferenceProbability extends ReferenceProbability {

    final ReferenceProbability myBase;
    final Translate myTranslate;

    FilterReferenceProbability(ReferenceProbability referenceProbability, Translate translate) {
        super(translate.numTaxa(), translate.numSites());
        myBase = referenceProbability;
        myTranslate = translate;
    }

    @Override
    public float value(int taxon, int site) {
        long taxonSite = myTranslate.taxonSite(taxon, site);
        if (taxonSite == -1) {
            return 0.0f;
        }
        return myBase.value((int) (taxonSite >>> 32), (int) (taxonSite & 0xFFFFFFFF));
    }

}
