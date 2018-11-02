package net.maizegenetics.dna.map;

import net.maizegenetics.util.BitSet;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by zrm22 on 3/27/17.
 *
 * Interface is used to store GenomeSequences defined by a GATK generated GVCF file
 */
public interface GVCFGenomeSequence extends GenomeSequence {
    public Map<Chromosome, byte[]> getChrPosMap();
    public PositionList getGVCFPositions();
    public HashMap<Chromosome,ArrayList<ArrayList<Integer>>> getConsecutiveRegions();
    public BitSet getMaskBitSet();
    public void setMaskBitSet(BitSet newMaskBitSet);

    public BitSet getFilterBitSet();
    public void setFilterBitSet(BitSet newFilterBitSet);

    public void flipMaskBit(int index);
    public void flipFilterBit(int index);
    public void writeFASTA(String fileName);
    public HashMap<String, String> chromosomeSequenceAndStats(Chromosome chrom, int startSite, int lastSite);

    public HashMap<String,Integer> getPreviousRegionStats();
    public void resetCounters();
}
