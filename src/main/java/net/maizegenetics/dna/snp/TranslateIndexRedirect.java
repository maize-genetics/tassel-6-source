/*
 *  TranslateIndexRedirect
 * 
 *  Created on Dec 14, 2016
 */
package net.maizegenetics.dna.snp;

import java.util.Arrays;

/**
 * Translation redirects index to corresponding index in base genotype table.
 *
 * @author Terry Casstevens
 */
public class TranslateIndexRedirect extends TranslateIndex {

    protected final int[] myIndexRedirect;

    /**
     * Constructor
     *
     * @param indexRedirect redirect index to base index. Should be ordered.
     */
    TranslateIndexRedirect(int[] indexRedirect) {
        super(indexRedirect.length, true);
        myIndexRedirect = indexRedirect;
    }

    /**
     * Translates index to base index.
     *
     * @param index index
     * @return translated base index
     */
    @Override
    public int translate(int index) {
        return myIndexRedirect[index];
    }

    /**
     * Translates base index to this index. Uses binary search algorithm since
     * indices are ordered.
     *
     * @param index index
     * @return translated index
     */
    @Override
    public int reverseTranslateIndex(int index) {
        return Arrays.binarySearch(myIndexRedirect, index);
    }

}
