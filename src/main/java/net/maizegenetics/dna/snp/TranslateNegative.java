/*
 *  TranslateNegative
 * 
 *  Created on Dec 10, 2016
 */
package net.maizegenetics.dna.snp;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateNegative extends Translate {

    TranslateNegative(TranslateIndex translateTaxa, TranslateIndex translateSite) {
        super(translateTaxa, translateSite);
    }

    @Override
    public long taxonSite(int taxon, int site) {
        int newTaxon = taxon(taxon);
        if (newTaxon == -1) {
            return -1;
        }
        int newSite = site(site);
        if (newSite == -1) {
            return -1;
        }
        return (long) newTaxon << 32 | newSite;
    }

}
