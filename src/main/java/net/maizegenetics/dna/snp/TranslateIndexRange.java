/*
 *  TranslateIndexRange
 * 
 *  Created on Dec 14, 2016
 */
package net.maizegenetics.dna.snp;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateIndexRange extends TranslateIndex {

    private final int myRangeStart;
    private final int myRangeEnd;

    /**
     * Index Range Translation
     *
     * @param start start index (inclusive)
     * @param end end index (inclusive)
     */
    TranslateIndexRange(int start, int end) {
        super(end - start + 1, true);
        myRangeStart = start;
        myRangeEnd = end;
    }

    @Override
    public int translate(int index) {
        return index + myRangeStart;
    }

    @Override
    public int reverseTranslateIndex(int index) {
        return index - myRangeStart;
    }

}
