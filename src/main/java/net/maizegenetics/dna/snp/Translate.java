/*
 *  Translate
 * 
 *  Created on Dec 10, 2016
 */
package net.maizegenetics.dna.snp;

/**
 * This translates filtered taxa and sites. Also indicates masking, but this
 * implementation has no mask.
 *
 * @author Terry Casstevens
 */
public class Translate {

    private final TranslateIndex myTranslateTaxa;
    private final TranslateIndex myTranslateSite;

    Translate(TranslateIndex translateTaxa, TranslateIndex translateSite) {
        if (translateTaxa == null) {
            throw new IllegalArgumentException("Translate: init: translate taxa can't be null");
        }
        if (translateSite == null) {
            throw new IllegalArgumentException("Translate: init: translate site can't be null");
        }
        myTranslateTaxa = translateTaxa;
        myTranslateSite = translateSite;
    }

    public int taxon(int taxon) {
        return myTranslateTaxa.translate(taxon);
    }

    public int site(int site) {
        return myTranslateSite.translate(site);
    }

    public long taxonSite(int taxon, int site) {
        return (long) taxon(taxon) << 32 | site(site);
    }

    public boolean hasSiteTranslations() {
        return myTranslateSite.hasTranslations();
    }

    public boolean hasTaxaTranslations() {
        return myTranslateTaxa.hasTranslations();
    }

    public int[] siteTranslations() {
        return myTranslateSite.getTranslations();
    }

    public int[] taxaTranslations() {
        return myTranslateTaxa.getTranslations();
    }

    TranslateIndex translateTaxa() {
        return myTranslateTaxa;
    }

    TranslateIndex translateSite() {
        return myTranslateSite;
    }

    public int numTaxa() {
        return myTranslateTaxa.numIndices();
    }

    public int numSites() {
        return myTranslateSite.numIndices();
    }

}
