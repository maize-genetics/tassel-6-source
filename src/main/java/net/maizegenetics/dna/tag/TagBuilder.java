package net.maizegenetics.dna.tag;

import net.maizegenetics.dna.BaseEncoder;

import java.nio.ByteBuffer;

/**
 * Builder for tags that optimizes the memory footprint.
 * Tags are encoded in long with 2 bits per bp.  In GBS, we generally use only record the 64 or 96bp.
 *
 * Custom classes are generated with 32, 64 or 96bp tags (1, 2 or 3 longs) otherwise an array is used.
 *
 * @author Ed Buckler
 */
public class TagBuilder {
    private final long[] seq2Bit;
    private final short length;
    private boolean isReference=false;
    private String name=null;


    private TagBuilder(long[] seq2Bit, short length) {
        this.seq2Bit=seq2Bit;
        this.length=length;
    }

    public TagBuilder reference() {
        isReference=true;
        return this;
    }

    public TagBuilder name(String name) {
        this.name=name;
        return this;
    }

    public Tag build() {
        switch (seq2Bit.length) {
            case 0: return null;
            case 1: return new Tag1Long(seq2Bit,length, isReference, name);
            case 2: return new Tag2Long(seq2Bit,length, isReference, name);
            case 3: return new Tag3Long(seq2Bit,length, isReference, name);
            default: return new TagVarLong(seq2Bit,length, isReference, name);
        }
    }

    public static TagBuilder instance(long[] seq2Bit, short length) {
        return new TagBuilder(seq2Bit,length);
    }

    public static TagBuilder instance(byte[] seq2BitInBytes, short length) {
        int seqBitLength=seq2BitInBytes.length/8;
        long[] seq2Bit=new long[seqBitLength];
        ByteBuffer bb=ByteBuffer.wrap(seq2BitInBytes);
        for (int i = 0; i < seq2Bit.length; i++) {
            seq2Bit[i]=bb.getLong();
        }
        return new TagBuilder(seq2Bit,length);
    }

    public static TagBuilder instance(String sequence) {
        long[] seq2Bit = AbstractTag.getLongArrayFromSeq(sequence);
        if (seq2Bit == null) { 
        	seq2Bit = new long[0];
        }
        return new TagBuilder(seq2Bit,(short)sequence.length());
    }

    public static TagBuilder reverseComplement(Tag tag) {
        String revSequence = BaseEncoder.getReverseComplement(tag.sequence());
        if (revSequence == null) return null;
        return instance(revSequence);
    }
}

class Tag1Long extends AbstractTag {
    //memory 8 + 8 + 1 =17 bytes
    //An array would add 12 to it
    private final long val0;

    Tag1Long(long[] val, short length, boolean reference, String name) {
        super(length,reference, name);
        val0=val[0];
    }

    @Override
    public long[] seq2Bit() {
        return new long[]{val0};
    }
}

class Tag2Long extends AbstractTag {
    //memory 8 + 16 + 1 =25 bytes
    //An array would add 12 to it
    private final long val0, val1;

    Tag2Long(long[] val, short length, boolean reference, String name) {
        super(length, reference, name);
        val0 = val[0];
        val1 = val[1];
    }

    @Override
    public long[] seq2Bit() {
        return new long[]{val0,val1};
    }
}

class Tag3Long extends AbstractTag {
    //memory 8 + 24 + 1 = 33 bytes
    //An array would add 12 to it
    private long val0, val1, val2;

    Tag3Long(long[] val, short length, boolean reference, String name) {
        super(length,reference, name);
        val0=val[0];
        val1=val[1];
        val2=val[2];
    }

    @Override
    public long[] seq2Bit() {
        return new long[]{val0,val1,val2};
    }

}

class TagVarLong extends AbstractTag {
    //memory 8 + 12 + 8*LongLen + 2 = XX bytes
    private long[] val;

    TagVarLong(long[] val, short length, boolean reference, String name) {
        super(length,reference, name);
        this.val=val;
    }

    @Override
    public long[] seq2Bit() {
        return val;
    }

}
