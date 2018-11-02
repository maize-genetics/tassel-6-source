package net.maizegenetics.dna.snp;

/**
 * @author Terry Casstevens
 *         Created March 30, 2017
 */
public class TranslateIndexOffsets extends TranslateIndex {

    private final int[] myOffsets;

    /**
     * Constructor.
     *
     * @param offsets offsets should start with zero and be the beginning index of each offset.
     *                The last offset is the total number of indices.
     */
    public TranslateIndexOffsets(int[] offsets) {
        super(offsets[offsets.length - 1], true);
        int last = 0;
        for (int current : offsets) {
            if (current <= last) {
                throw new IllegalArgumentException("TranslateIndexOffsets: init: offsets should be accending and not zero.");
            }
            last = current;
        }
        myOffsets = offsets;
    }

    /**
     * Translates index to base index.
     *
     * @param index index
     * @return translated base index
     */
    @Override
    public int translate(int index) {
        for (int i = 1; i < myOffsets.length; i++) {
            if (index < myOffsets[i]) {
                return i - 1;
            }
        }
        throw new IndexOutOfBoundsException("TranslateIndexOffsets: translate: " + index);
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
        return myOffsets[index];
    }

}
