/*
 * GenotypeTableUtils
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.factor.FactorTableKt;
import org.apache.log4j.Logger;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_GENOTYPE;

/**
 * Utility methods for comparing, sorting, and counting genotypes.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class GenotypeTableUtils {

    private static final Logger myLogger = Logger.getLogger(GenotypeTableUtils.class);

    private static final byte HIGHMASK = (byte) 0x0F;

    private GenotypeTableUtils() {
        // utility class
    }

    /**
     * Returns whether diploid allele values are heterozygous. First 4 bits in
     * byte is one allele value. Second 4 bits is other allele value.
     *
     * @param diploidAllele alleles
     *
     * @return true if allele values different; false if values the same.
     */
    public static boolean isHeterozygous(byte diploidAllele) {
        return ((diploidAllele >>> 4) & 0xf) != (diploidAllele & 0xf);
    }

    /**
     * Returns whether diploid allele values are homozygous. Unknown values
     * return false.
     *
     * @param diploidAllele
     *
     * @return true if allele values are the same; false if unknown or not equal
     */
    public static boolean isHomozygous(byte diploidAllele) {
        if (diploidAllele == GenotypeTable.UNKNOWN_GENOTYPE) {
            return false;
        }
        return ((diploidAllele >>> 4) & 0xf) == (diploidAllele & 0xf);
    }

    /**
     * Combines two allele values into one diploid value. Assumed phased.
     *
     * @param a allele 1
     * @param b allele 2
     *
     * @return diploid value
     */
    public static byte getDiploidValuePhased(byte a, byte b) {
        return (byte) ((a << 4) | (HIGHMASK & b));
    }

    /**
     * Combines two allele values into one diploid value. Assumed phased.
     *
     * @param a allele 1
     * @param b allele 2
     *
     * @return diploid value
     */
    public static byte getDiploidValue(byte a, byte b) {
        return getDiploidValuePhased(a, b);
    }

    /**
     * Combines two allele values into one diploid value. In alphabetical order
     *
     * @param a allele 1
     * @param b allele 2
     *
     * @return diploid value sorted by order A < C < G < T
     */
    public static byte getUnphasedDiploidValue(byte a, byte b) {
        a = (byte) (HIGHMASK & a);
        b = (byte) (HIGHMASK & b);
        if (a < b) {
            return (byte) ((a << 4) | b);
        }
        return (byte) ((b << 4) | a);
    }

    /**
     * Combines two genotype values into one diploid value. Returns unknown if
     * either parent is heterozygous or unknown, or alleles are swapped.
     *
     * @param g1 genotype 1
     * @param g2 genotype 2
     *
     * @return diploid value
     */
    public static byte getUnphasedDiploidValueNoHets(byte g1, byte g2) {
        if ((g2 == g1) && (!isHeterozygous(g1))) {
            return g1;
        }
        if (g1 == UNKNOWN_GENOTYPE) {
            return UNKNOWN_GENOTYPE;
        }
        if (g2 == UNKNOWN_GENOTYPE) {
            return UNKNOWN_GENOTYPE;
        }
        if (isHeterozygous(g1)) {
            return UNKNOWN_GENOTYPE;
        }
        if (isHeterozygous(g2)) {
            return UNKNOWN_GENOTYPE;
        }

        return getUnphasedDiploidValue(g1, g2);
    }

    /**
     * Separates diploid allele value into it's two values.
     *
     * @param genotype diploid value
     *
     * @return separated allele values
     */
    public static byte[] getDiploidValues(byte genotype) {
        byte[] result = new byte[2];
        result[0] = (byte) ((genotype >>> 4) & 0xf);
        result[1] = (byte) (genotype & 0xf);
        if (result[0] == 0xF) result[0] = FactorTableKt.UNKNOWN_ALLELE;
        if (result[1] == 0xF) result[1] = FactorTableKt.UNKNOWN_ALLELE;
        return result;
    }

}
