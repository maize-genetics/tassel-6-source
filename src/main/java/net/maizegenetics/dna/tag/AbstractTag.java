package net.maizegenetics.dna.tag;

import com.google.common.collect.ComparisonChain;
import net.maizegenetics.dna.BaseEncoder;

import java.nio.ByteBuffer;
import java.util.Arrays;

/**
 * Created by edbuckler on 7/26/14.
 */
public abstract class AbstractTag implements Tag, Comparable<Tag> {
    private final short length;
    private final boolean reference;
    private final String name;

    protected AbstractTag(short length, boolean reference, String name) {
        this.length=length;
        this.reference=reference;
        this.name=name;
    }

    @Override
    public short seqLength() {
        return length;
    }

    @Override
    public boolean isReference() {
        return reference;
    }

    @Override
    public String name() {
        if(name==null) return "";
        return name;
    }

    @Override
    public String sequence() {
        return getSequenceFromLong(seq2Bit(),seqLength());
    }

    @Override
    public byte[] seq2BitAsBytes() {
        long[] seqInLong=seq2Bit();
        ByteBuffer b= ByteBuffer.allocate(seqInLong.length*8);
        for (long l : seqInLong) {
            b.putLong(l);
        }
        return b.array();
    }

    @Override
    public int compareTo(Tag o) {
        long[] t=this.seq2Bit();
        long[] to=o.seq2Bit();
        for (int i = 0; i < t.length; i++) {
            int c=Long.compare(t[0],to[1]);
            if(c!=0) return c;
        }
        return Short.compare(this.seqLength(),o.seqLength());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Tag that = (Tag) o;

        if (seqLength() != that.seqLength()) return false;
        if (!Arrays.equals(this.seq2Bit(), that.seq2Bit())) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(seq2Bit());
        result = 31 * result + (int) seqLength();
        return result;
    }

    @Override
    public String toString() {
        return "Tag{" +
                "seq=" + sequence() +
                ", length=" + seqLength() +
                ", Ref=" + isReference() +
                "}";
    }



    /**
     * Return a string representation of an array of 2-bit encoded longs.
     * @param val array of 2-bit encoded sequences
     * @return DNA sequence as a string
     */
    protected static String getSequenceFromLong(long[] val, int length) {
        StringBuilder seq = new StringBuilder();
        for (long v : val) {
            //System.out.println(BaseEncoder.getSequenceFromLong(v));
            seq.append(BaseEncoder.getSequenceFromLong(v));
        }
        return seq.toString().substring(0,length);
    }

    /**
     * @param seq A String containing a DNA sequence.
     * @return result A array of Long containing the binary representation of the sequence.
     * null if sequence length has any non-DNA characters
     */
    protected static long[] getLongArrayFromSeq(String seq) {
        final int chunkSize=32;
        int longsNeeded = (seq.length() + chunkSize - 1) / chunkSize;  //this is doing the integer version of Math.ceil
        long[] result = new long[longsNeeded];
        for (int i = 0; i < result.length; i++) {
            result[i] = BaseEncoder.getLongFromSeq(seq.substring(i * chunkSize, Math.min((i + 1) * chunkSize,seq.length())));
            if(result[i]==-1) return null;
        }
        return result;
    }
}
