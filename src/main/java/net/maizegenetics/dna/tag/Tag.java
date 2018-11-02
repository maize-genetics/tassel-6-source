package net.maizegenetics.dna.tag;

import net.maizegenetics.dna.BaseEncoder;

import java.io.Serializable;

/**
 * Interface for tags (sequences without N) that captures these bit encoded sequence and there length.
 * 
 * @author Ed Buckler
 */
public interface Tag {

    /*
    Sequence of the tag (A,C,G,T) - no ambiguity possible
     */
    String sequence();

    /*
    Name of the tag (or sequence).  Generally used for contig sequences, but left blank if GBS tag
     */
    String name();

    long[] seq2Bit();

    byte[] seq2BitAsBytes();

    short seqLength();

    boolean isReference();

    default String toCSVString() {
        return sequence() + "," + seqLength();
    }

    default String toReverseComplement() {
        return BaseEncoder.getReverseComplement(sequence());
    }

}
