/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import net.maizegenetics.dna.tag.Tag;

/**
 * This class is used by RepGenAlignerPlugin to store alignment
 * info to db table tagAlignments.  It is also used when pulling
 * alignments from the DB.  The tag2chrom and tag2pos fields are
 * used to determine if the tag2 of this alignment class is
 * a reference tag.  If tag2chrom is null and tag2pos = -1, the
 * tag alignment info is a non-ref tag.  If these fields are
 * populated with good values, the tag2 alignment info is for
 * a reference tag.
 * 
 * The "alignmentPos" field indicates the position within tag2
 * where the tag1 alignment starts.  This position is adjusted
 * for any clipping of tag1 that occurred during alignment.
 * 
 * @author lcj34
 *
 */
public class AlignmentInfo implements Comparable<AlignmentInfo>{
    private final Tag tag2;
    private  final String tag2chrom; // this field is null is tag2 is NOT a reference tag
    private  final int tag2pos;
    private  final int alignmentPos;
    private final int ref_strand; // forward/plus=1, reverse/minus=0, unknown = -1,(mostly) as per Position interface
    private final String ref_genome; // reference genome from which a ref tag was derived
    private  final int myScore;

    public AlignmentInfo(Tag tag2, String chromosome, int position, int alignmentpos, int ref_strand, String ref_genome, int score) {
        this.tag2 = tag2;
        this.tag2chrom = chromosome;
        this.tag2pos = position;
        this.alignmentPos = alignmentpos; // the positions in tag2 where alignment of tag1 starts.
        this.ref_strand = ref_strand;
        this.ref_genome = ref_genome;
        this.myScore = score;
    }

    public Tag tag2() {
        return tag2;
    }
    public String tag2chrom() {
        return tag2chrom;
    }

    public int tag2pos() {
        return tag2pos;
    }
    
    public int alignmentPos() {
        return alignmentPos;
    }
    public int ref_strand() {
        return ref_strand;
    }
    public String ref_genome() {
        return ref_genome;
    }
    public  int score() {
        return myScore;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("Alignment:");
        sb.append("\tTag2:").append(tag2.sequence());
        sb.append("\tChr:").append(tag2chrom);
        sb.append("\tPos:").append(tag2pos);
        sb.append("\tAlignmentPos:").append(alignmentPos);
        sb.append("\tScore:").append(myScore);
        sb.append("\n");
        return sb.toString();
    }
    @Override
    public int compareTo(AlignmentInfo other) {
        // TODO Auto-generated method stub
        return this.myScore > other.score() ? 1: this.myScore < other.score() ? -1 : 0;
    }
}
