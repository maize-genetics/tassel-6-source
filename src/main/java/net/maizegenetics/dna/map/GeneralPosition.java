/*
 * GeneralPosition
 */
package net.maizegenetics.dna.map;

import com.google.common.base.Joiner;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;

import java.util.Arrays;

/**
 * Provide information on a site and its annotations. This includes information
 * on position, MAF, coverage. This class is immutable.
 * <p>
 * </p>
 * The annotations are all set using the builder.
 * <p>
 * </p>
 * This class has been optimized for memory size compared to the other version.
 * Every annotated position takes about 84 bytes, even with a names and declared
 * variants. Positions with extended annotation take roughly 117 bytes, and as
 * long as the annotation are repeated extensively the take an additional 8
 * bytes each. So a site with 100 annotation would take <1000 bytes.
 *
 * @author Ed Buckler
 */
public final class GeneralPosition implements Position {

    /**
     * Locus of the site (required)
     */
    private final Chromosome myChromosome;
    /**
     * Physical position of the site (required)
     */
    private final int myPosition;
    /**
     * Insertion-position for any insertions relative to reference.
     */
    private final short myInsertionPosition;
    /**
     * Strand of the site (unknown = Byte.MIN_VALUE)
     */
    private final byte myStrand;
    /**
     * Is type Nucleotide or Text
     */
    private final boolean isNucleotide;
    /**
     * Whether the variant define the nature of the indel
     */
    private final boolean isIndel;
    private final float myMAF;
    private final float mySiteCoverage;
    private final long myAlleleValue;  //todo consider TAS-376
    /**
     * Name of the site (default = SLocus_Position)
     */
    private final byte[] mySNPIDAsBytes;
    /**
     * Define the nature of the polymorphism {"ACTAT","-"} or {"A","C","G"} or
     * {"100","103","106"}
     */
    private final GeneralAnnotation myVariantsAndAnno;
    private final int hashCode;

    /**
     * A builder for creating immutable CoreAnnotatedPosition instances.
     * AnnotatedPositions are built off a base of a CorePosition, so build it
     * first.
     * <p>
     * Example:
     * <pre>   {@code
     * Position cp= new CorePosition.Builder(new Chromosome("1"),1232).build();
     * CoreAnnotatedPosition ap= new CoreAnnotatedPosition.Builder(cp)
     *    .maf(0.05f)
     *    .ancAllele(NucleotideAlignmentConstants.C_ALLELE)
     *    .build();}</pre>
     * <p>
     * This would create nucleotide position on chromosome 1 at position 1232.
     * The MAF is 0.05 and the ancestral allele is C.
     */
    public static class Builder {

        // Required parameters
        private Chromosome myChromosome;
        private int myPosition;
        // Optional parameters - initialized to default values
        private short myInsertionPosition = 0;
        private byte myStrand = Position.STRAND_PLUS;
        private String mySNPID = null;
        private boolean isNucleotide = true;
        private boolean isIndel = false;

        //in an allele annotation objects
        private float myMAF = Float.NaN;
        private float mySiteCoverage = Float.NaN;
        private byte[] myAlleles = new byte[WHICH_ALLELE.COUNT];
        private long myAllelesAsLong;
        private final GeneralAnnotationStorage.Builder myAnnotations;

        /**
         * Constructor requires a Position before annotation of the position
         *
         * @see Position#of(String, int)
         * @see Position#builder(String, int)
         */
        public Builder(Chromosome chr, int position) {
            this(chr, position, GeneralAnnotationStorage.getBuilder());
        }

        public Builder(Chromosome chr, int position, GeneralAnnotationStorage.Builder annotationBuilder) {
            myChromosome = Chromosome.getCanonicalChromosome(chr);
            myPosition = position;
            Arrays.fill(myAlleles, GenotypeTable.UNKNOWN_ALLELE);
            myAnnotations = annotationBuilder;
        }

        /**
         * Constructor from an existing position
         */
        public Builder(Position aCorePosition) {
            this(aCorePosition.getChromosome(), aCorePosition.getPosition());
            myInsertionPosition = aCorePosition.getInsertionPosition();
            myStrand = aCorePosition.getStrand();
            mySNPID = aCorePosition.getSNPID();
            isNucleotide = aCorePosition.isNucleotide();
            isIndel = aCorePosition.isIndel();
            myMAF = aCorePosition.getGlobalMAF();
            mySiteCoverage = aCorePosition.getGlobalSiteCoverage();
            for (WHICH_ALLELE alleleType : WHICH_ALLELE.values()) {
                myAlleles[alleleType.index()] = aCorePosition.getAllele(alleleType);
            }
            myAnnotations.addAnnotations(aCorePosition.getAnnotation());

        }

        /**
         * Set Chromosome
         */
        public Builder chromosome(Chromosome val) {
            myChromosome = val;
            return this;
        }

        /**
         * Set Position in Chromosome
         */
        public Builder position(int val) {
            myPosition = val;
            return this;
        }

        /**
         * Set Insertion-position
         */
        public Builder insertionPosition(short val) {
            myInsertionPosition = val;
            return this;
        }

        /**
         * Set strand (default=1)
         */
        public Builder strand(byte val) {
            myStrand = val;
            return this;
        }

        public Builder strand(String val) {
            myStrand = Position.getStrand(val);
            return this;
        }

        /**
         * Set SNP name (default="S"+Chromosome+"_"+position)
         */
        public Builder snpName(String val) {
            mySNPID = val;
            return this;
        }

        /**
         * Set whether position is nucleotide (default=true)
         */
        public Builder nucleotide(boolean val) {
            isNucleotide = val;
            return this;
        }

        /**
         * Set whether position is indel (default=false)
         */
        public Builder indel(boolean val) {
            isIndel = val;
            return this;
        }

        /**
         * Set text definition of variants (default=null)
         */
        public Builder knownVariants(String[] val) {
            StringBuilder sb = new StringBuilder();
            for (String s : val) {
                sb.append(s).append("/");
            }
            sb.setLength(sb.length() - 1); //remove terminal /
            myAnnotations.addAnnotation("VARIANT", sb.toString());
            return this;
        }

        /**
         * Set text definition of variants (default=null)
         */
        public Builder knownVariants(String val) {
            myAnnotations.addAnnotation("VARIANT", val);
            return this;
        }

        /**
         * Set Minor Allele Frequency annotation (default=Float.NaN)
         */
        public Builder maf(float val) {
            myMAF = val;
            return this;
        }

        /**
         * Set site coverage annotation (default=Float.NaN)
         */
        public Builder siteCoverage(float val) {
            mySiteCoverage = val;
            return this;
        }

        /**
         * Set allele annotation by Allele type
         * (default=Alignment.UNKNOWN_ALLELE)
         */
        public Builder allele(WHICH_ALLELE aT, byte val) {
            myAlleles[aT.index()] = val;
            return this;
        }

        /**
         * Add non-standard annotation, boolean and String are converted to String with toString, Number are stored as
         * Number
         */
        public Builder addAnno(String key, Object value) {
            if (value instanceof Number) {
                myAnnotations.addAnnotation(key, (Number) value);
            } else {
                myAnnotations.addAnnotation(key, value.toString());
            }
            return this;
        }

        /**
         * Add non-standard annotation with key-value separated by '='
         */
        public Builder addAnno(String keyValue) {
            String[] sub = keyValue.split("=");
            if (sub.length == 1) {
                return addAnno(keyValue, "TRUE");
            }
            return addAnno(sub[0], sub[1]);
        }

        public GeneralPosition build() {
            for (int i = myAlleles.length - 1; i >= 0; i--) {
                myAllelesAsLong = (myAllelesAsLong << 4) | myAlleles[i];
            }
            if (mySNPID != null) {
                String defaultS = (new StringBuilder("S").append(myChromosome.getName()).append("_").append(myPosition)).toString();
                if (defaultS.equals(mySNPID)) {
                    mySNPID = null;
                }
            }
            return new GeneralPosition(this);
        }
    }

    private GeneralPosition(Builder builder) {
        myChromosome = builder.myChromosome;
        myPosition = builder.myPosition;
        myInsertionPosition = builder.myInsertionPosition;
        myStrand = builder.myStrand;
        if (builder.mySNPID == null) {
            mySNPIDAsBytes = null;
        } else {
            mySNPIDAsBytes = builder.mySNPID.getBytes();
        }
        isNucleotide = builder.isNucleotide;
        isIndel = builder.isIndel;

        myVariantsAndAnno = builder.myAnnotations.build();

        hashCode = calcHashCode();

        myMAF = builder.myMAF;
        mySiteCoverage = builder.mySiteCoverage;
        myAlleleValue = builder.myAllelesAsLong;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) {
            return true;
        }
        if (!(obj instanceof Position)) {
            return false;
        }
        Position o = (Position) obj;
        int result = compareTo(o);
        return result == 0;
    }

    @Override
    public int compareTo(Position o) {
        int result = myChromosome.compareTo(o.getChromosome());
        if (result != 0) {
            return result;
        }
        result = Integer.compare(myPosition, o.getPosition());
        if (result != 0) {
            return result;
        }
        result = Short.compare(myInsertionPosition, o.getInsertionPosition());
        if (result != 0) {
            return result;
        }
        result = Byte.compare(myStrand, o.getStrand());
        if (result != 0) {
            return result;
        }
        if (getActualSNPID() == null && o.getActualSNPID() == null) {
            return 0;
        }
        return getSNPID().compareTo(o.getSNPID());
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("Position");
        sb.append("\tChr:").append(getChromosome().getName());
        sb.append("\tPos:").append(getPosition());
        sb.append("\tInsertionPos:").append(getInsertionPosition());
        sb.append("\tName:").append(getSNPID());
        if (myVariantsAndAnno != null) {
            String[] variants = myVariantsAndAnno.getTextAnnotation("VARIANT");
            if ((variants != null) && (variants.length != 0)) {
                sb.append("\tVariants:").append(Joiner.on(",").join(variants));
            }
        }
        sb.append("\tMAF:").append(getGlobalMAF());
        sb.append("\tRef:").append(NucleotideAlignmentConstants.getHaplotypeNucleotide(getAllele(WHICH_ALLELE.Reference)));
        return sb.toString();
    }

    private int calcHashCode() {
        int hash = 7;
        hash = 37 * hash + this.myChromosome.hashCode();
        hash = 37 * hash + this.myPosition;
        hash = 37 * hash + this.myInsertionPosition;
        hash = 37 * hash + this.myStrand;
        if (mySNPIDAsBytes != null) {
            hash = 37 * hash + Arrays.hashCode(this.mySNPIDAsBytes);
        }
        return hash;
    }

    @Override
    public float getGlobalMAF() {
        return myMAF;
    }

    @Override
    public float getGlobalSiteCoverage() {
        return mySiteCoverage;
    }

    @Override
    public byte getAllele(WHICH_ALLELE alleleType) {
        return (byte) ((myAlleleValue >> (alleleType.index() * 4)) & 0xF);
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    @Override
    public Chromosome getChromosome() {
        return myChromosome;
    }

    @Override
    public int getPosition() {
        return myPosition;
    }

    @Override
    public short getInsertionPosition() {
        return myInsertionPosition;
    }

    @Override
    public byte getStrand() {
        return myStrand;
    }

    @Override
    public String getStrandStr() {
        return Position.getStrand(myStrand);
    }

    @Override
    public String getSNPID() {
        if (mySNPIDAsBytes == null) {
            return (new StringBuilder("S").append(getChromosome().getName()).append("_").append(myPosition)).toString();
        } else {
            return new String(mySNPIDAsBytes);
        }
    }

    @Override
    public String getActualSNPID() {
        if (mySNPIDAsBytes == null) {
            return null;
        } else {
            return new String(mySNPIDAsBytes);
        }
    }

    @Override
    public boolean isNucleotide() {
        return isNucleotide;
    }

    @Override
    public boolean isIndel() {
        return isIndel;
    }

    @Override
    public String[] getKnownVariants() {
        if (myVariantsAndAnno == null) {
            return new String[0];
        }
        if (myVariantsAndAnno.getTextAnnotation("VARIANT").length == 0) {
            return new String[0];
        }
        return myVariantsAndAnno.getTextAnnotation("VARIANT")[0].replace("[", "").replace("]", "").split("/");
    }

    @Override
    public GeneralAnnotation getAnnotation() {
        return myVariantsAndAnno;
    }
}
