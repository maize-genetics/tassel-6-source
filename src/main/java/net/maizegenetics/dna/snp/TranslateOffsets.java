package net.maizegenetics.dna.snp;

/**
 * @author Terry Casstevens
 *         Created March 30, 2017
 */
public class TranslateOffsets {

    private final int[] myOffsets;
    private final int myNumIndices;

    /**
     * @param offsets offsets should start with zero and be the beginning index of each offset.
     *                The last offset is the total number of indices.
     */
    public TranslateOffsets(int[] offsets) {
        if (offsets[0] != 0) {
            throw new IllegalArgumentException("TranslateOffsets: init: first offset must be zero.");
        }
        for (int i = 1; i < offsets.length; i++) {
            if (offsets[i] <= offsets[i - 1]) {
                throw new IllegalArgumentException("TranslateOffsets: init: each offset must be larger than the previous.");
            }
        }
        myOffsets = offsets;
        myNumIndices = offsets[offsets.length - 1];
    }

    public long translate(int index) {
        for (int i = 1; i < myOffsets.length; i++) {
            if (index < myOffsets[i]) {
                return ((long) (i - 1) << 32) | (long) (index - myOffsets[i - 1]);
            }
        }
        throw new IndexOutOfBoundsException("TranslateOffsets: translate: " + index);
    }

}
