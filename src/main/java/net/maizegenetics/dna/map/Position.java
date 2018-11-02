package net.maizegenetics.dna.map;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.util.GeneralAnnotation;

/**
 * Defines a genomic positions and its known variants. Includes attributes of
 * chromosome, position, sub-position, strand, name (or SNP ID), whether this
 * position is a nucleotide, or includes an indel.
 *
 * @author Ed Buckler
 */
public interface Position extends Comparable<Position> {

    public static final byte STRAND_PLUS = (byte) 1;
    public static final byte STRAND_MINUS = (byte) 0;
    public static final byte STRAND_UNKNOWN = Byte.MIN_VALUE;

    public static final String STRAND_PLUS_STR = "+";
    public static final String STRAND_MINUS_STR = "-";
    public static final String STRAND_UNKNOWN_STR = "N";

    /**
     * Create Position given chromosome number and position.
     *
     * @param chromosomeNumber chromosome number
     * @param positionWithinChromosome physical position
     *
     * @return Position
     */
    public static Position of(int chromosomeNumber, int positionWithinChromosome) {
        return new GeneralPosition.Builder(Chromosome.instance(chromosomeNumber), positionWithinChromosome).build();
    }

    /**
     * Create Position given chromosome name and position.
     *
     * @param chromosomeName chromosome name
     * @param positionWithinChromosome physical position
     *
     * @return Position
     */
    public static Position of(String chromosomeName, int positionWithinChromosome) {
        return new GeneralPosition.Builder(Chromosome.instance(chromosomeName), positionWithinChromosome).build();
    }

    /**
     * Create Position given chromosome and position.
     *
     * @param chromosome chromosome
     * @param positionWithinChromosome physical position
     *
     * @return Position
     */
    public static Position of(Chromosome chromosome, int positionWithinChromosome) {
        return new GeneralPosition.Builder(chromosome, positionWithinChromosome).build();
    }

    /**
     * Create Position Builder given chromosome number and position.
     *
     * @param chromosomeNumber chromosome number
     * @param positionWithinChromosome physical position
     *
     * @return Position Builder
     */
    public static GeneralPosition.Builder builder(int chromosomeNumber, int positionWithinChromosome) {
        return new GeneralPosition.Builder(Chromosome.instance(chromosomeNumber), positionWithinChromosome);
    }

    /**
     * Create Position Builder given chromosome name and position.
     *
     * @param chromosomeName chromosome name
     * @param positionWithinChromosome physical position
     *
     * @return Position Builder
     */
    public static GeneralPosition.Builder builder(String chromosomeName, int positionWithinChromosome) {
        return new GeneralPosition.Builder(Chromosome.instance(chromosomeName), positionWithinChromosome);
    }

    public static String getStrand(byte value) {
        switch (value) {
            case STRAND_PLUS:
                return STRAND_PLUS_STR;
            case STRAND_MINUS:
                return STRAND_MINUS_STR;
            case STRAND_UNKNOWN:
                return STRAND_UNKNOWN_STR;
            default:
                throw new IllegalStateException("Position: getStrand: unknown strand value: " + value);
        }
    }

    public static byte getStrand(String value) {
        switch (value) {
            case STRAND_PLUS_STR:
                return STRAND_PLUS;
            case STRAND_MINUS_STR:
                return STRAND_MINUS;
            case STRAND_UNKNOWN_STR:
                return STRAND_UNKNOWN;
            default:
                throw new IllegalStateException("Position: getStrand: unknown strand value: " + value);
        }
    }

    /**
     * Return the locus (generally a chromosome) of a site
     */
    public Chromosome getChromosome();

    /**
     * Return the physical position of a site
     */
    public int getPosition();

    /**
     * Return the insertion-position of a site. This will be zero for the main
     * physical position and sequentially numbered (1, 2, 3,...) for
     * any (if any) sites with the same {@link #getPosition()}.
     * This will mostly be used for insertions relative to the reference.
     */
    public short getInsertionPosition();

    /**
     * Return the strand for a site definition
     */
    public byte getStrand();

    public String getStrandStr();

    /**
     * Return the ID (name) for a site
     */
    public String getSNPID();

    /**
     * Returns SNP ID only if assigned. getSNPID() returns a default if not
     * assigned.
     *
     * @return SNP ID
     */
    public String getActualSNPID();

    /**
     * Whether the position is a nucleotide position or another marker type
     * (SSR, AFLP, RAPD, CNV, which are recoded with text states)
     */
    public boolean isNucleotide();

    /**
     * Whether the position includes indels, which would be defined in the
     * variants
     */
    public boolean isIndel();

    /**
     * Returns the nature of the polymorphism {"ACTAT","-"} or {"A","C","G"} or
     * {"100","103","106"}
     */
    public String[] getKnownVariants();

    /**
     * Return the minor allele frequency in a global scope
     */
    public float getGlobalMAF();

    /**
     * Returns the proportion of genotypes scored at a given site
     */
    public float getGlobalSiteCoverage();

    /**
     * Return the allele specified by alleleType, if unknown Alignment.Unknown
     * is return
     */
    public byte getAllele(WHICH_ALLELE alleleType);

    public GeneralAnnotation getAnnotation();

}
