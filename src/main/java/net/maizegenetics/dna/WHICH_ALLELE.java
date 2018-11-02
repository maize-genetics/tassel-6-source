package net.maizegenetics.dna;

/**
 * This defines the possible alleles.
 */
public enum WHICH_ALLELE {

    /**
     * Major Allele - Most frequent allele.
     */
    Major(0),
    /**
     * Minor Allele - Second most frequent allele.
     */
    Minor(1),
    /**
     * Global Major Allele
     */
    GlobalMajor(2),
    /**
     * Global Minor Allele
     */
    GlobalMinor(3),
    /**
     * Reference Allele
     */
    Reference(4),
    /**
     * Alternate to Reference Allele
     */
    Alternate(5),
    /**
     * High Coverage Allele
     */
    HighCoverage(6),
    /**
     * Low Coverage Allele
     */
    LowCoverage(7),
    /**
     * Remaining Minor Alleles
     */
    Minor2(8), Minor3(9), Minor4(10), Minor5(11),
    /**
     * Ancestral Allele
     */
    Ancestral(12),
    /**
     * Diploid Unknown
     */
    Unknown(13)
    ;

    private final int myIndex;
    /**
     * Count of the number of allele types
     */
    public final static int COUNT = WHICH_ALLELE.values().length;

    WHICH_ALLELE(int index) {
        myIndex = index;
    }

    /**
     * Sequential index that can be use for primitive arrays
     */
    public int index() {
        return myIndex;
    }

    private static WHICH_ALLELE[] FREQ_ALLELES = new WHICH_ALLELE[]{Major, Minor, Minor2, Minor3, Minor4, Minor5};

    public static WHICH_ALLELE[] frequencyAlleles() {
        return FREQ_ALLELES;
    }
}
