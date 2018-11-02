package net.maizegenetics.dna.map;

import com.google.common.collect.ComparisonChain;

/**
 * Defines the haplotypes and positions of a chromosomal segment
 *
 * @author Ed Buckler
 */
public class DonorHaplotypes implements Comparable<DonorHaplotypes>{
    private final Chromosome chromosome;
    private final int startPosition;
    private final int endPosition;
    private final int parent1index;
    private final int parent2index;

    public DonorHaplotypes(Chromosome chromosome, int startPosition, int endPosition, int parent1index, int parent2index) {
        this.chromosome=chromosome;
        this.startPosition=startPosition;
        this.endPosition=endPosition;
        this.parent1index=parent1index;
        this.parent2index=parent2index;
    }

    public static DonorHaplotypes getMergedInstance(DonorHaplotypes dh1, DonorHaplotypes dh2) {
        if(!dh1.getChromosome().equals(dh2.getChromosome())) return null;
        if(dh1.getParent1index()!=dh2.getParent1index()) return null;
        if(dh1.getParent2index()!=dh2.getParent2index()) return null;
        int startPosition=Math.min(dh1.getStartPosition(), dh2.getStartPosition());
        int endPosition=Math.max(dh1.getEndPosition(), dh2.getEndPosition());
        return new DonorHaplotypes(dh1.getChromosome(),startPosition,endPosition,dh1.getParent1index(),dh1.getParent2index());
    }

    public Chromosome getChromosome() {
        return chromosome;
    }

    public int getStartPosition() {
        return startPosition;
    }

    public int getEndPosition() {
        return endPosition;
    }

    public int getParent1index() {
        return parent1index;
    }

    public int getParent2index() {
        return parent2index;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) {return true;}
        if (!(obj instanceof DonorHaplotypes)) {return false;}
        DonorHaplotypes o=(DonorHaplotypes)obj;
        if(compareTo(o)!=0) return false;
        return true;
    }

    @Override
    public int compareTo(DonorHaplotypes o) {
        int result= ComparisonChain.start()
                .compare(chromosome, o.getChromosome())
                .compare(startPosition, o.getStartPosition())
                .compare(endPosition, o.getEndPosition())
                .compare(parent1index, o.getParent1index())
                .compare(parent2index, o.getParent2index())
                .result();
        return result;
    }

    @Override
    public String toString() {
        return "DonorHaplotypes{"+
                "chromosome="+chromosome+
                ", startPosition="+startPosition+
                ", endPosition="+endPosition+
                ", parent1index="+parent1index+
                ", parent2index="+parent2index+
                '}';
    }
}
