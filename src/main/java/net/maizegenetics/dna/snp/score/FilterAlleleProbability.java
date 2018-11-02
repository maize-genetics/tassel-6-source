/*
 *  FilterAlleleProbability
 * 
 *  Created on Dec 16, 2016
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.Translate;

/**
 *
 * @author Terry Casstevens
 */
public class FilterAlleleProbability extends AlleleProbability {

    final AlleleProbability myBase;
    final Translate myTranslate;

    FilterAlleleProbability(AlleleProbability alleleProbability, Translate translate) {
        super(translate.numTaxa(), translate.numSites());
        if (alleleProbability instanceof FilterAlleleProbability) {
            throw new IllegalArgumentException();
        }
        myBase = alleleProbability;
        myTranslate = translate;
    }

    @Override
    public float value(int taxon, int site, SiteScore.SITE_SCORE_TYPE scoreType) {
        long taxonSite = myTranslate.taxonSite(taxon, site);
        if (taxonSite == -1) {
            return 0.0f;
        }
        return myBase.value((int) (taxonSite >>> 32), (int) (taxonSite & 0xFFFFFFFF), scoreType);
    }

}
