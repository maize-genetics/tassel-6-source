/*
 * NucleotideAlignmentConstants
 */
package net.maizegenetics.dna.snp;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

/**
 *
 * @author terry
 */
public final class NucleotideAlignmentConstants {

    // Byte Values for Nucleotide Alleles
    public static final byte A_ALLELE = (byte) 0x0;
    public static final byte C_ALLELE = (byte) 0x1;
    public static final byte G_ALLELE = (byte) 0x2;
    public static final byte T_ALLELE = (byte) 0x3;
    public static final byte INSERT_ALLELE = (byte) 0x4;
    public static final byte GAP_ALLELE = (byte) 0x5;
    public static final byte UNDEFINED_ALLELE = (byte) 0x6;

    // Diploid Byte Values for Nucleotide Alleles
    public static final byte A_HOMOZYGOUS = (byte) 0x00;
    public static final byte C_HOMOZYGOUS = (byte) 0x11;
    public static final byte G_HOMOZYGOUS = (byte) 0x22;
    public static final byte T_HOMOZYGOUS = (byte) 0x33;
    public static final byte INSERT_HOMOZYGOUS = (byte) 0x44;
    public static final byte GAP_HOMOZYGOUS = (byte) 0x55;
    public static final byte UNDEFINED_HOMOZYGOUS = (byte) 0x66;

    // String Values for Nucleotide Alleles
    public static final String INSERT_ALLELE_STR = "+";
    public static final String GAP_ALLELE_STR = "-";
    public static final char GAP_ALLELE_CHAR = '-';
    public static final String UNDEFINED_ALLELE_STR = "X";

    public static final String[][] NUCLEOTIDE_ALLELES = new String[][]{{"A", "C", "G", "T", "+", "-",
        UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR,
        UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, GenotypeTable.RARE_ALLELE_STR, GenotypeTable.UNKNOWN_ALLELE_STR}};
    /**
     * Number of nucleotide states excluding rare and unknown.
     */
    public static final int NUMBER_NUCLEOTIDE_ALLELES = 6;
    private static final Map<String, Byte> NUCLEOTIDE_DIPLOID_HASH = new HashMap<>();

    static {
        NUCLEOTIDE_DIPLOID_HASH.put("AA", (byte) 0x00);
        NUCLEOTIDE_DIPLOID_HASH.put("AC", (byte) 0x01);
        NUCLEOTIDE_DIPLOID_HASH.put("AG", (byte) 0x02);
        NUCLEOTIDE_DIPLOID_HASH.put("AT", (byte) 0x03);
        NUCLEOTIDE_DIPLOID_HASH.put("A+", (byte) 0x04);
        NUCLEOTIDE_DIPLOID_HASH.put("A-", (byte) 0x05);
        NUCLEOTIDE_DIPLOID_HASH.put("AN", (byte) 0x0F);
        NUCLEOTIDE_DIPLOID_HASH.put("AX", (byte) 0x0F);
        NUCLEOTIDE_DIPLOID_HASH.put("AZ", (byte) 0x0E);

        NUCLEOTIDE_DIPLOID_HASH.put("CA", (byte) 0x10);
        NUCLEOTIDE_DIPLOID_HASH.put("CC", (byte) 0x11);
        NUCLEOTIDE_DIPLOID_HASH.put("CG", (byte) 0x12);
        NUCLEOTIDE_DIPLOID_HASH.put("CT", (byte) 0x13);
        NUCLEOTIDE_DIPLOID_HASH.put("C+", (byte) 0x14);
        NUCLEOTIDE_DIPLOID_HASH.put("C-", (byte) 0x15);
        NUCLEOTIDE_DIPLOID_HASH.put("CN", (byte) 0x1F);
        NUCLEOTIDE_DIPLOID_HASH.put("CX", (byte) 0x1F);
        NUCLEOTIDE_DIPLOID_HASH.put("CZ", (byte) 0x1E);

        NUCLEOTIDE_DIPLOID_HASH.put("GA", (byte) 0x20);
        NUCLEOTIDE_DIPLOID_HASH.put("GC", (byte) 0x21);
        NUCLEOTIDE_DIPLOID_HASH.put("GG", (byte) 0x22);
        NUCLEOTIDE_DIPLOID_HASH.put("GT", (byte) 0x23);
        NUCLEOTIDE_DIPLOID_HASH.put("G+", (byte) 0x24);
        NUCLEOTIDE_DIPLOID_HASH.put("G-", (byte) 0x25);
        NUCLEOTIDE_DIPLOID_HASH.put("GN", (byte) 0x2F);
        NUCLEOTIDE_DIPLOID_HASH.put("GX", (byte) 0x2F);
        NUCLEOTIDE_DIPLOID_HASH.put("GZ", (byte) 0x2E);

        NUCLEOTIDE_DIPLOID_HASH.put("TA", (byte) 0x30);
        NUCLEOTIDE_DIPLOID_HASH.put("TC", (byte) 0x31);
        NUCLEOTIDE_DIPLOID_HASH.put("TG", (byte) 0x32);
        NUCLEOTIDE_DIPLOID_HASH.put("TT", (byte) 0x33);
        NUCLEOTIDE_DIPLOID_HASH.put("T+", (byte) 0x34);
        NUCLEOTIDE_DIPLOID_HASH.put("T-", (byte) 0x35);
        NUCLEOTIDE_DIPLOID_HASH.put("TN", (byte) 0x3F);
        NUCLEOTIDE_DIPLOID_HASH.put("TX", (byte) 0x3F);
        NUCLEOTIDE_DIPLOID_HASH.put("TZ", (byte) 0x3E);

        NUCLEOTIDE_DIPLOID_HASH.put("+A", (byte) 0x40);
        NUCLEOTIDE_DIPLOID_HASH.put("+C", (byte) 0x41);
        NUCLEOTIDE_DIPLOID_HASH.put("+G", (byte) 0x42);
        NUCLEOTIDE_DIPLOID_HASH.put("+T", (byte) 0x43);
        NUCLEOTIDE_DIPLOID_HASH.put("++", (byte) 0x44);
        NUCLEOTIDE_DIPLOID_HASH.put("+-", (byte) 0x45);
        NUCLEOTIDE_DIPLOID_HASH.put("+N", (byte) 0x4F);
        NUCLEOTIDE_DIPLOID_HASH.put("+X", (byte) 0x4F);
        NUCLEOTIDE_DIPLOID_HASH.put("+Z", (byte) 0x4E);

        NUCLEOTIDE_DIPLOID_HASH.put("-A", (byte) 0x50);
        NUCLEOTIDE_DIPLOID_HASH.put("-C", (byte) 0x51);
        NUCLEOTIDE_DIPLOID_HASH.put("-G", (byte) 0x52);
        NUCLEOTIDE_DIPLOID_HASH.put("-T", (byte) 0x53);
        NUCLEOTIDE_DIPLOID_HASH.put("-+", (byte) 0x54);
        NUCLEOTIDE_DIPLOID_HASH.put("--", (byte) 0x55);
        NUCLEOTIDE_DIPLOID_HASH.put("-N", (byte) 0x5F);
        NUCLEOTIDE_DIPLOID_HASH.put("-X", (byte) 0x5F);
        NUCLEOTIDE_DIPLOID_HASH.put("-Z", (byte) 0x5E);

        NUCLEOTIDE_DIPLOID_HASH.put("NA", (byte) 0xF0);
        NUCLEOTIDE_DIPLOID_HASH.put("NC", (byte) 0xF1);
        NUCLEOTIDE_DIPLOID_HASH.put("NG", (byte) 0xF2);
        NUCLEOTIDE_DIPLOID_HASH.put("NT", (byte) 0xF3);
        NUCLEOTIDE_DIPLOID_HASH.put("N+", (byte) 0xF4);
        NUCLEOTIDE_DIPLOID_HASH.put("N-", (byte) 0xF5);
        NUCLEOTIDE_DIPLOID_HASH.put("NN", (byte) 0xFF);
        NUCLEOTIDE_DIPLOID_HASH.put("NX", (byte) 0xFF);
        NUCLEOTIDE_DIPLOID_HASH.put("NZ", (byte) 0xFE);

        NUCLEOTIDE_DIPLOID_HASH.put("XA", (byte) 0xF0);
        NUCLEOTIDE_DIPLOID_HASH.put("XC", (byte) 0xF1);
        NUCLEOTIDE_DIPLOID_HASH.put("XG", (byte) 0xF2);
        NUCLEOTIDE_DIPLOID_HASH.put("XT", (byte) 0xF3);
        NUCLEOTIDE_DIPLOID_HASH.put("X+", (byte) 0xF4);
        NUCLEOTIDE_DIPLOID_HASH.put("X-", (byte) 0xF5);
        NUCLEOTIDE_DIPLOID_HASH.put("XN", (byte) 0xFF);
        NUCLEOTIDE_DIPLOID_HASH.put("XX", (byte) 0xFF);
        NUCLEOTIDE_DIPLOID_HASH.put("XZ", (byte) 0xFE);

        NUCLEOTIDE_DIPLOID_HASH.put("ZA", (byte) 0xE0);
        NUCLEOTIDE_DIPLOID_HASH.put("ZC", (byte) 0xE1);
        NUCLEOTIDE_DIPLOID_HASH.put("ZG", (byte) 0xE2);
        NUCLEOTIDE_DIPLOID_HASH.put("ZT", (byte) 0xE3);
        NUCLEOTIDE_DIPLOID_HASH.put("Z+", (byte) 0xE4);
        NUCLEOTIDE_DIPLOID_HASH.put("Z-", (byte) 0xE5);
        NUCLEOTIDE_DIPLOID_HASH.put("ZN", (byte) 0xEF);
        NUCLEOTIDE_DIPLOID_HASH.put("ZX", (byte) 0xEF);
        NUCLEOTIDE_DIPLOID_HASH.put("ZZ", (byte) 0xEE);

        NUCLEOTIDE_DIPLOID_HASH.put("A", (byte) 0x00); // AA
        NUCLEOTIDE_DIPLOID_HASH.put("C", (byte) 0x11); // CC
        NUCLEOTIDE_DIPLOID_HASH.put("G", (byte) 0x22); // GG
        NUCLEOTIDE_DIPLOID_HASH.put("T", (byte) 0x33); // TT
        NUCLEOTIDE_DIPLOID_HASH.put("+", (byte) 0x44); // ++
        NUCLEOTIDE_DIPLOID_HASH.put("-", (byte) 0x55); // --
        NUCLEOTIDE_DIPLOID_HASH.put("Z", (byte) 0xEE); // ZZ
        NUCLEOTIDE_DIPLOID_HASH.put("N", (byte) 0xFF); // NN
        NUCLEOTIDE_DIPLOID_HASH.put("X", (byte) 0xFF); // NN

        NUCLEOTIDE_DIPLOID_HASH.put("R", (byte) 0x02); // AG
        NUCLEOTIDE_DIPLOID_HASH.put("Y", (byte) 0x13); // CT
        NUCLEOTIDE_DIPLOID_HASH.put("S", (byte) 0x21); // GC
        NUCLEOTIDE_DIPLOID_HASH.put("W", (byte) 0x03); // AT
        NUCLEOTIDE_DIPLOID_HASH.put("K", (byte) 0x23); // GT
        NUCLEOTIDE_DIPLOID_HASH.put("M", (byte) 0x01); // AC
        NUCLEOTIDE_DIPLOID_HASH.put("0", (byte) 0x54); // -+
    }
    private static final byte[] NUCLEOTIDE_DIPLOID_ARRAY = new byte[256];

    static {
        Arrays.fill(NUCLEOTIDE_DIPLOID_ARRAY, UNDEFINED_HOMOZYGOUS);
        for (String temp : NUCLEOTIDE_DIPLOID_HASH.keySet()) {
            NUCLEOTIDE_DIPLOID_ARRAY[getNucleotideDiploidArrayIndex(temp)] = NUCLEOTIDE_DIPLOID_HASH.get(temp);
        }
    }
    private static final int mask = 0x2F;
    private static final int mask2 = 0x81;
    private static final int shift = 2;

    private static int getNucleotideDiploidArrayIndex(String str) {

        if (str.length() == 1) {
            return str.charAt(0);
        } else if (str.length() == 2) {
            return ((((str.charAt(1) << shift) ^ (byte) mask2)) ^ (str.charAt(0) & (byte) mask)) & 0xFF;
        } else {
            throw new IllegalStateException("NucleotideAlignmentConstants: getNucleotideDiploidArrayIndex: illegal str: " + str);
        }

    }
    public static final Map<Byte, String> NUCLEOTIDE_IUPAC_HASH = new HashMap<Byte, String>();

    static {
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x00, "A"); // AA
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x01, "M"); // AC
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x02, "R"); // AG
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x03, "W"); // AT
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x04, "0"); // A+
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x05, "0"); // A-
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x0E, "A"); // AZ
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x0F, "A"); // AN

        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x10, "M"); // CA
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x11, "C"); // CC
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x12, "S"); // CG
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x13, "Y"); // CT
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x14, "0"); // C+
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x15, "0"); // C-
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x1E, "C"); // CZ
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x1F, "C"); // CN

        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x20, "R"); // GA
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x21, "S"); // GC
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x22, "G"); // GG
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x23, "K"); // GT
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x24, "0"); // G+
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x25, "0"); // G-
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x2E, "G"); // GZ
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x2F, "G"); // GN

        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x30, "W"); // TA
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x31, "Y"); // TC
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x32, "K"); // TG
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x33, "T"); // TT
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x34, "0"); // T+
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x35, "0"); // T-
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x3E, "T"); // TZ
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x3F, "T"); // TN

        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x40, "0"); // +A
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x41, "0"); // +C
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x42, "0"); // +G
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x43, "0"); // +T
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x44, "+"); // ++
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x45, "0"); // +-
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x4E, "+"); // +Z
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x4F, "+"); // +N

        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x50, "0"); // -A
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x51, "0"); // -C
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x52, "0"); // -G
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x53, "0"); // -T
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x54, "0"); // -+
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x55, "-"); // --
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x5E, "-"); // -Z
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x5F, "-"); // -N

        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE0, "A"); // ZA
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE1, "C"); // ZC
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE2, "G"); // ZG
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE3, "T"); // ZT
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE4, "+"); // Z+
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE5, "-"); // Z-
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xEE, "Z"); // ZZ
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xEF, "N"); // ZN

        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF0, "A"); // NA
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF1, "C"); // NC
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF2, "G"); // NG
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF3, "T"); // NT
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF4, "+"); // N+
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF5, "-"); // N-
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xFE, "N"); // NZ
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xFF, "N"); // NN

    }
    private static final String[] NUCLEOTIDE_IUPAC_ARRAY = new String[256];
    private static final char[] NUCLEOTIDE_IUPAC_CHARARRAY = new char[256];
    private static final char[] NUCLEOTIDE_IUPAC_CHARCOMPLEMENT_ARRAY = new char[256];

    static {
        Arrays.fill(NUCLEOTIDE_IUPAC_ARRAY, UNDEFINED_ALLELE_STR);
        Arrays.fill(NUCLEOTIDE_IUPAC_CHARARRAY, UNDEFINED_ALLELE_STR.charAt(0));
        Arrays.fill(NUCLEOTIDE_IUPAC_CHARCOMPLEMENT_ARRAY, UNDEFINED_ALLELE_STR.charAt(0));
        for (Byte temp : NUCLEOTIDE_IUPAC_HASH.keySet()) {
            NUCLEOTIDE_IUPAC_ARRAY[temp & 0xFF] = NUCLEOTIDE_IUPAC_HASH.get(temp);
            NUCLEOTIDE_IUPAC_CHARARRAY[temp & 0xFF] = NUCLEOTIDE_IUPAC_HASH.get(temp).charAt(0);
            byte compByte = getNucleotideDiploidComplement(temp);
            NUCLEOTIDE_IUPAC_CHARCOMPLEMENT_ARRAY[NUCLEOTIDE_IUPAC_HASH.get(temp).charAt(0)] = NUCLEOTIDE_IUPAC_HASH.get(compByte).charAt(0);
        }
    }
    private static final Map<String, Byte> NUCLEOTIDE_ALLELE_HASH = new HashMap<String, Byte>();
    private static final byte[] NUCLEOTIDE_ALLELE_ARRAY = new byte[256];

    static {
        NUCLEOTIDE_ALLELE_HASH.put("A", A_ALLELE); // A
        NUCLEOTIDE_ALLELE_HASH.put("C", C_ALLELE); // C
        NUCLEOTIDE_ALLELE_HASH.put("G", G_ALLELE); // G
        NUCLEOTIDE_ALLELE_HASH.put("T", T_ALLELE); // T
        NUCLEOTIDE_ALLELE_HASH.put("a", A_ALLELE); // A
        NUCLEOTIDE_ALLELE_HASH.put("c", C_ALLELE); // C
        NUCLEOTIDE_ALLELE_HASH.put("g", G_ALLELE); // G
        NUCLEOTIDE_ALLELE_HASH.put("t", T_ALLELE); // T
        NUCLEOTIDE_ALLELE_HASH.put("+", INSERT_ALLELE); // +
        NUCLEOTIDE_ALLELE_HASH.put("-", GAP_ALLELE); // -
        NUCLEOTIDE_ALLELE_HASH.put("N", GenotypeTable.UNKNOWN_ALLELE); // N
        NUCLEOTIDE_ALLELE_HASH.put("n", GenotypeTable.UNKNOWN_ALLELE); // N
        Arrays.fill(NUCLEOTIDE_ALLELE_ARRAY, UNDEFINED_ALLELE);
        for (Map.Entry<String, Byte> en : NUCLEOTIDE_ALLELE_HASH.entrySet()) {
            NUCLEOTIDE_ALLELE_ARRAY[en.getKey().charAt(0)] = en.getValue();
        }
    }

    private NucleotideAlignmentConstants() {
        // do not instantiate
    }

    /**
     * Returns diploid byte value for given nucleotide value. First four bits
     * contain first allele value. And second four bits contain second allele
     * value.
     *
     * @param value diploid allele value
     *
     * @return nucleotide diploid allele byte value
     */
    public static byte getNucleotideDiploidByte(String value) {
        try {
            return NUCLEOTIDE_DIPLOID_ARRAY[getNucleotideDiploidArrayIndex(value)];
            // return NUCLEOTIDE_DIPLOID_HASH.get(value).byteValue();
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("NucleotideAlignmentConstants: getNucleotideDiploidByte: unknown allele value: " + value);
        }
    }

    /**
     * Returns Optional diploid byte value for given nucleotide value. First four bits
     * contain first allele value. And second four bits contain second allele
     * value. Optional empty if not present in the Hash.
     *
     * @param value diploid allele value
     *
     * @return nucleotide diploid allele byte value
     */
    public static Optional<Byte> parseNucleotideDiploidByte(String value) {
        try {
            return Optional.ofNullable(NUCLEOTIDE_DIPLOID_ARRAY[getNucleotideDiploidArrayIndex(value)]);
        } catch (IllegalStateException e) {
            return Optional.empty();
        }
    }

    /**
     * Returns char allele for allele byte encoding. THis is called from junit
     * tests in ReferenceGenomeSequenceTest.java
     *
     * @param value haploid allele byte value
     *
     * @return nucleotide haploid allele byte value
     */
    public static Object getNucleotideAlleleValue(byte value) {
        try {
            for (Object alleleObject : NUCLEOTIDE_ALLELE_HASH.keySet()) {
                if (NUCLEOTIDE_ALLELE_HASH.get(alleleObject).equals(value)) {
                    return alleleObject;
                }
            }
            return null;
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("NucleotideAlignmentConstants: getNucleotideAlleleByte: unknown allele value: " + value);
        }
    }

    /**
     * Returns haploid byte value for given nucleotide value. Only right-most
     * four bits used.
     *
     * @param value haploid allele value
     *
     * @return nucleotide haploid allele byte value
     */
    public static byte getNucleotideAlleleByte(String value) {
        try {
            return NUCLEOTIDE_ALLELE_HASH.get(value).byteValue();
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("NucleotideAlignmentConstants: getNucleotideAlleleByte: unknown allele value: " + value);
        }
    }

    /**
     * Returns haploid byte value for given nucleotide value. Only right-most
     * four bits used.
     *
     * @param value haploid allele value
     *
     * @return nucleotide haploid allele byte value
     */
    public static byte getNucleotideAlleleByte(char value) {
        try {
            return NUCLEOTIDE_ALLELE_ARRAY[value];
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("NucleotideAlignmentConstants: getNucleotideAlleleByte: unknown allele value: " + value);
        }
    }

    /**
     * Returns diploid byte value for given nucleotide value. First four bits
     * contain first allele value. And second four bits contain second allele
     * value.
     *
     * @param value diploid allele value
     *
     * @return nucleotide diploid allele byte value
     */
    public static byte getNucleotideDiploidByte(char value) {
        try {
            return NUCLEOTIDE_DIPLOID_ARRAY[value];
            // return NUCLEOTIDE_DIPLOID_HASH.get(String.valueOf(value)).byteValue();
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("NucleotideAlignmentConstants: getNucleotideDiploidByte: unknown allele value: " + value);
        }
    }

    /**
     * Returns the IUPAC String for the given diploid allele value.
     *
     * @param value diploid allele value
     *
     * @return IUPAC String
     */
    public static String getNucleotideIUPAC(byte value) {
        return NUCLEOTIDE_IUPAC_ARRAY[value & 0xFF];
    }

    /**
     * Returns the IUPAC char for the given diploid allele value.
     *
     * @param value diploid allele value
     *
     * @return IUPAC char
     */
    public static char getNucleotideIUPACChar(byte value) {
        return NUCLEOTIDE_IUPAC_CHARARRAY[value & 0xFF];
    }

    /**
     * Returns the Nucleotide String for the given haploid allele value.
     *
     * @param value haploid value
     *
     * @return Nucleotide String
     */
    public static String getHaplotypeNucleotide(byte value) {
        return NUCLEOTIDE_ALLELES[0][value];
    }

    /**
     * Returns the Nucleotide Complement of given byte encoded nucleotide. A
     * returns T. T returns A. C returns G. G returns C. Otherwise given
     * nucleotide is returned.
     *
     * @param nucleotide nucleotide byte value
     *
     * @return Nucleotide Complement
     */
    public static byte getNucleotideComplement(byte nucleotide) {

        if (nucleotide == A_ALLELE) {
            return T_ALLELE;
        } else if (nucleotide == T_ALLELE) {
            return A_ALLELE;
        } else if (nucleotide == C_ALLELE) {
            return G_ALLELE;
        } else if (nucleotide == G_ALLELE) {
            return C_ALLELE;
        } else {
            return nucleotide;
        }

    }

    /**
     * Returns the Nucleotide Complement of the given diploid byte encoded
     * alleles.
     *
     * @param diploidAllele diploid allele value
     *
     * @return Nucleotide Complement
     */
    public static byte getNucleotideDiploidComplement(byte diploidAllele) {

        byte first = (byte) ((diploidAllele >>> 4) & 0xf);
        byte second = (byte) (diploidAllele & 0xf);
        first = getNucleotideComplement(first);
        second = getNucleotideComplement(second);
        return (byte) ((first << 4) | second);
    }

    /**
     * Returns the Nucleotide Complement of the given diploid IUPAC encoded
     * alleles.
     *
     * @param diploidAllele diploid IUPAC
     *
     * @return Nucleotide Complement in IUPAC
     */
    public static char getNucleotideDiploidIUPACComplement(char diploidAllele) {
        if (diploidAllele > 256) {
            return UNDEFINED_ALLELE_STR.charAt(0);
        }
        return NUCLEOTIDE_IUPAC_CHARCOMPLEMENT_ARRAY[diploidAllele];
    }

    /**
     * Returns whether given allele encodings are for Nucleotide Data.
     *
     * @param alleleStates allele encodings
     *
     * @return true if nucleotide encodings
     */
    public static boolean isNucleotideEncodings(String[][] alleleStates) {

        boolean isNucleotide = false;
        if (alleleStates.length == 1) {
            isNucleotide = true;
            if (alleleStates[0].length == NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0].length) {
                for (int i = 0; i < alleleStates.length; i++) {
                    if (!alleleStates[0][i].equals(NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0][i])) {
                        isNucleotide = false;
                    }
                }
            }

        }

        return isNucleotide;

    }

    /**
     * Convert a haploid (allele) string (e.g. ACGTA) to arrays of allele bytes
     * (e.g. {0,1,2,3,0})
     *
     * @param haploString haploid allele string
     * @return encoded array of bytes
     */
    public static byte[] convertHaplotypeStringToAlleleByteArray(String haploString) {
        byte[] haplosBytes = new byte[haploString.length()];
        for (int i = 0; i < haplosBytes.length; i++) {
            haplosBytes[i] = getNucleotideAlleleByte(haploString.charAt(i));
        }
        return haplosBytes;
    }

    /**
     * Convert a genotype (diploid) string (e.g. ACGTA) to arrays of diploid bytes
     * (e.g. {0,17,34,51,0})
     *
     * @param diploidString haploid allele string
     * @return encoded array of bytes
     */
    public static byte[] convertGenotypeStringToDiploidByteArray(String diploidString) {
        byte[] genoBytes = new byte[diploidString.length()];
        for (int i = 0; i < genoBytes.length; i++) {
            genoBytes[i] = getNucleotideDiploidByte(diploidString.charAt(i));
        }
        return genoBytes;
    }

    /**
     * Convert a haploid (allele) string (e.g. {0,1,2,3,0}) to its reverse
     * complement (e.g. {3,0,1,2,3})
     *
     * @param alleles haploid allele byte
     * @return encoded array of bytes
     */
    public static byte[] reverseComplementAlleleByteArray(byte[] alleles) {
        byte[] reverseCompBytes = new byte[alleles.length];
        int allelesIndex = alleles.length - 1;
        for (int i = 0; i < reverseCompBytes.length; i++) {
            reverseCompBytes[i] = getNucleotideComplement(alleles[allelesIndex--]);
        }
        return reverseCompBytes;
    }

    /**
     * Convert a nucleotide byte array to a string
     *
     * @param b the array of nucleotide bytes to convert
     * @return a string representing the nucleotide bytes
     */
    public static String nucleotideBytetoString(byte[] b) {
        StringBuilder result = new StringBuilder();
        for (int index = 0; index < b.length; index++) {
            result.append(getNucleotideAlleleValue(b[index]));
        }
        return result.toString();
    }

    /**
     * Return whether a simple homozgyous DNA state (i.e. no indels or N).
     * @param diploidAllele diploid
     * @return
     */
    public static boolean isHomozygousACGT(Byte diploidAllele) {
        if (!GenotypeTableUtils.isHomozygous(diploidAllele)) return false;
        if((diploidAllele & 0xf) < INSERT_ALLELE) return true;
        return false;
    }

    /**
     * Return whether a simple DNA state (i.e. no indels or N).  Heterozygous resolves to true also.
     * @param diploidAllele diploid
     * @return
     */
    public static boolean isACGT(Byte diploidAllele) {
        if ((((diploidAllele >>> 4) & 0xf) < INSERT_ALLELE) && ((diploidAllele & 0xf) < INSERT_ALLELE)) return true;
        return false;
    }
}
