/*
 *  TranslateIndex
 * 
 *  Created on Dec 14, 2016
 */
package net.maizegenetics.dna.snp;

/**
 * No translation to index.
 *
 * @author Terry Casstevens
 */
public class TranslateIndex {

    private final int myNumIndices;
    private final boolean myHasTranslations;

    /**
     * Constructor
     *
     * @param numIndices number of sites
     */
    TranslateIndex(int numIndices, boolean hasTranslations) {
        myNumIndices = numIndices;
        myHasTranslations = hasTranslations;
    }

    /**
     * Translates index to base index. This class has no translation.
     *
     * @param index index
     * @return translated base index
     */
    public int translate(int index) {
        return index;
    }

    /**
     * Translates base index to this index. This class has no translation.
     *
     * @param index index
     * @return translated index
     */
    public int reverseTranslateIndex(int index) {
        return index;
    }

    /**
     * Number of indices represented by this translation. Number of base indices
     * will be the same or larger.
     *
     * @return number of indices
     */
    public int numIndices() {
        return myNumIndices;
    }

    /**
     * Indicates whether this index translation has translations. If every index
     * is translated to itself, there are no translations. If this represents a
     * subset of the total indices, there are translations even if all indices
     * in the subset translate to itself.
     *
     * @return whether this has translations
     */
    public boolean hasTranslations() {
        return myHasTranslations;
    }

    /**
     * Returns all the translated indices.
     *
     * @return all translated indices
     */
    public int[] getTranslations() {
        int[] result = new int[numIndices()];
        for (int i = 0; i < numIndices(); i++) {
            result[i] = translate(i);
        }
        return result;
    }

}
