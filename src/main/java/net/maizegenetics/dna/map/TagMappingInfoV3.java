
package net.maizegenetics.dna.map;

/**
 * Container class for information on mapping position.  

 * <p>
 * The mappingSource can be an identifier. Byte.MIN_VALUE means the tag doesn't exist (missing) in the mapping result. Other values means it exist in mapping result. The result could be either mapped or unmapped.
 * 
 * Note: BWA doesn't provide mapping score, so the rank is just based on the order provided for multiple hits
 * Note: Currently, perfectMatch of PE is unknown (Byte.MIN_VALUE);
 * 
 * @author edbuckler Fei Lu
 */
public class TagMappingInfoV3 {
    /** Chromosome as an integer, unknown = Integer.MIN_VALUE */
    public int chromosome=Integer.MIN_VALUE;  // 4 bytes
    /** Strand relative to reference genome. 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE */
    public byte strand=Byte.MIN_VALUE;   // 1 byte
    /**Chromosomal position of the barcoded end of the tag, unknown = Integer.MIN_VALUE */
    public int startPosition=Integer.MIN_VALUE;  //   // 4 bytes
    /**Chromosomal position of the common adapter end of the tag (smaller than startPosition if tag matches minus strand), inclusive, unknown = Integer.MIN_VALUE */
    public int endPosition=Integer.MIN_VALUE;  //   // 4 bytes
    /**Number of diverging bp (edit distance) from reference, unknown = Byte.MIN_VALUE*/
    public byte divergence=Byte.MIN_VALUE;
    /**If the alignment is a perfect match to the reference, y = 1, n = 0, unknown = Byte.MIN_VALUE*/
    public byte perfectMatch = Byte.MIN_VALUE;
    /**Code of mappingSource.0: Bowtie2; 1: BWA; 2: BLAST; 3: BWAMEM; 4: PE one end; 5: PE the other end; 6: Genetic Mapping; missing = Byte.MIN_VALUE*/
    public byte mappingSource = Byte.MIN_VALUE;
    /**The rank of this mapping based on the scores from one aligner, starting from 0. If there are two 0, the rank of the third one is 1*/
    public byte mappingRank = Byte.MIN_VALUE;
    /**The mapping score of this mapping, unknown = Byte.MIN_VALUE*/
    public short mappingScore = Short.MIN_VALUE;
    /**Double cross-over probability Round(Log2(P)), unknown = Byte.MIN_VALUE */
    public byte dcoP=Byte.MIN_VALUE;
    /**Genetic mapping probability Round(Log2(P)), unknown= Byte.MIN_VALUE */
    public byte mapP=Byte.MIN_VALUE;  //Round(Log2(P)), unknown Byte.MIN_VALUE
    
    /**
     * MappingSource definition
     */
    public static enum Aligner {
        Bowtie2((byte)0, "Bowtie2"),
        BWA((byte)1, "BWA"),
        Blast((byte)2, "Blast"),
        BWAMEM((byte)3,"BWAMEM"),
        PEEnd1((byte)4, "PEEnd1"),
        PEEnd2((byte)5, "PEEnd2");
        private byte value;
        private String name;
        Aligner (byte mappingSource, String alignerName) {
            value = mappingSource;
            name = alignerName;
        }
        public byte getValue () {
            return value;
        }
        public String getName () {
            return name;
        }
        public static byte getValueFromName (String name) {
            if (name.equalsIgnoreCase(Bowtie2.name)) return Bowtie2.getValue();
            if (name.equalsIgnoreCase(BWA.name)) return BWA.getValue();
            if (name.equalsIgnoreCase(Blast.name)) return Blast.getValue();
            if (name.equalsIgnoreCase(BWAMEM.name)) return BWAMEM.getValue();
            if (name.equalsIgnoreCase(PEEnd1.name)) return PEEnd1.getValue();
            if (name.equalsIgnoreCase(PEEnd2.name)) return PEEnd2.getValue();
            return Byte.MIN_VALUE;
        }
        public static Aligner getAlignerFromName (String name) {
            if (name.equalsIgnoreCase(Bowtie2.name)) return Bowtie2;
            if (name.equalsIgnoreCase(BWA.name)) return BWA;
            if (name.equalsIgnoreCase(Blast.name)) return Blast;
            if (name.equalsIgnoreCase(BWAMEM.name)) return BWAMEM;
            if (name.equalsIgnoreCase(PEEnd1.name)) return PEEnd1;
            if (name.equalsIgnoreCase(PEEnd2.name)) return PEEnd2;
            return null;
        }
        public static Aligner getAlignerFromValue (byte value) {
            if (value == Bowtie2.getValue()) return Bowtie2;
            if (value == BWA.getValue()) return BWA;
            if (value == Blast.getValue()) return Blast;
            if (value == BWAMEM.getValue()) return BWAMEM;
            if (value == PEEnd1.getValue()) return PEEnd1;
            if (value == PEEnd2.getValue()) return PEEnd2;
            return null;
        }
    }

    public TagMappingInfoV3() {   
    }
    
    public TagMappingInfoV3(int chromosome, byte strand, int startPosition, int endPosition, byte divergence, byte perfectMatch, byte mappingSource, short mappingScore) {
        this.chromosome=chromosome;
        this.strand=strand;
        this.startPosition=startPosition;
        this.endPosition=endPosition;
        this.divergence=divergence;
        this.perfectMatch = perfectMatch;
        this.mappingSource = mappingSource;
        this.mappingScore = mappingScore;
    } 
    
    public TagMappingInfoV3(int chromosome, byte strand, int startPosition, 
                int endPosition, byte divergence, byte perfectMatch, byte mappingSource, short mappingScore, byte[] variantPosOff, byte[] variantDef,
                byte dcoP, byte mapP) {
        this(chromosome, strand, startPosition, endPosition, divergence, perfectMatch, mappingSource, mappingScore);
        this.dcoP=dcoP;
        this.mapP=mapP;
    }
    
    public void setMappingRank (byte mappingRank) {
        this.mappingRank = mappingRank;
    }
    
    public void setMappingSource (byte mappingSource) {
        this.mappingSource = mappingSource;
    }
}
