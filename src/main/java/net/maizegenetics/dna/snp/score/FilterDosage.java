/*
 *  FilterDosage
 * 
 *  Created on Dec 16, 2016
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.Translate;

/**
 *
 * @author Terry Casstevens
 */
public class FilterDosage extends Dosage {

    final Dosage myBase;
    final Translate myTranslate;

    FilterDosage(Dosage dosage, Translate translate) {
        super(translate.numTaxa(), translate.numSites());
        myBase = dosage;
        myTranslate = translate;
    }

    @Override
    public byte value(int taxon, int site) {
        long taxonSite = myTranslate.taxonSite(taxon, site);
        if (taxonSite == -1) {
            return 0;
        }
        return myBase.value((int) (taxonSite >>> 32), (int) (taxonSite & 0xFFFFFFFF));
    }

}
