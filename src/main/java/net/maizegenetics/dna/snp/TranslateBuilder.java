/*
 *  TranslateBuilder
 * 
 *  Created on Dec 10, 2016
 */
package net.maizegenetics.dna.snp;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateBuilder {

    private TranslateBuilder() {
    }

    /**
     *
     * @param translateTaxa translate taxa
     * @param translateSite translate site
     *
     * @return Translate
     */
    public static Translate getInstance(TranslateIndex translateTaxa, TranslateIndex translateSite) {

        if (translateTaxa == null) {
            throw new IllegalArgumentException("TranslateBuilder: filter: must specify translateTaxa");
        }

        if (translateSite == null) {
            throw new IllegalArgumentException("TranslateBuilder: filter: must specify translateSite");
        }

        if (hasNegativeIndices(translateTaxa) || hasNegativeIndices(translateSite)) {
            return new TranslateNegative(translateTaxa, translateSite);
        } else {
            return new Translate(translateTaxa, translateSite);
        }

    }

    public static Translate getInstance(Translate base, Translate translate) {
        TranslateIndex translateTaxa = TranslateIndexBuilder.merge(base.translateTaxa(), translate.translateTaxa());
        TranslateIndex translateSite = TranslateIndexBuilder.merge(base.translateSite(), translate.translateSite());
        if (hasNegativeIndices(translateTaxa) || hasNegativeIndices(translateSite)) {
            return new TranslateNegative(translateTaxa, translateSite);
        } else {
            return new Translate(translateTaxa, translateSite);
        }
    }

    private static boolean hasNegativeIndices(TranslateIndex translateIndex) {

        if (!(translateIndex instanceof TranslateIndexRedirectUnordered)) {
            return false;
        }

        for (int current : translateIndex.getTranslations()) {
            if (current == -1) {
                return true;
            }
        }

        return false;

    }

}
