/*
 *  TranslateIndexRedirectUnordered
 * 
 *  Created on Dec 14, 2016
 */
package net.maizegenetics.dna.snp;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateIndexRedirectUnordered extends TranslateIndexRedirect {

    /**
     * Constructor
     *
     * @param indexRedirect redirected index indices to base index. These are
     * unordered.
     */
    TranslateIndexRedirectUnordered(int[] indexRedirect) {
        super(indexRedirect);
    }

    /**
     * Translates base index to this index. Searches all indices since not
     * ordered.
     *
     * @param index index
     * @return translated index
     */
    @Override
    public int reverseTranslateIndex(int index) {
        for (int i = 0; i < numIndices(); i++) {
            if (myIndexRedirect[i] == index) {
                return i;
            }
        }
        return -1;
    }

}
