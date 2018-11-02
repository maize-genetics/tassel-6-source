/*
 *  LineIndex
 * 
 *  Created on Aug 29, 2015
 */

package net.maizegenetics.dna.snp.io;

/**
 *
 * @author Terry Casstevens
 */
public class LineIndex {

    public static final int NUM_LINES_PER_INTERVAL = 10;

    private final int myMagic;
    private final char myCommentChar;
    private final int myNumHeaderLinesToSkip;
    private final int myNumLinesPerInterval;
    private final long[] myVirtualFileOffsets;

    public LineIndex(int magic, char commentChar, int numHeaderLinesToSkip, int numLinesPerInterval, long[] virtualFileOffsets) {
        myMagic = magic;
        myCommentChar = commentChar;
        myNumHeaderLinesToSkip = numHeaderLinesToSkip;
        myNumLinesPerInterval = numLinesPerInterval;
        myVirtualFileOffsets = virtualFileOffsets;
    }

    public int magicNumber() {
        return myMagic;
    }

    public char commentChar() {
        return myCommentChar;
    }

    public int numHeaderLinesToSkip() {
        return myNumHeaderLinesToSkip;
    }

    public int numLinesPerInterval() {
        return myNumLinesPerInterval;
    }

    public int numVirtualOffsets() {
        return myVirtualFileOffsets.length;
    }

    public long virtualOffset(int index) {
        return myVirtualFileOffsets[index];
    }

}

