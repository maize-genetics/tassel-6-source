package net.maizegenetics.dna.map;

import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Defines the chromosome structure and length. The name and length recorded for
 * each chromosome.
 *
 * @author Terry Casstevens and Ed Buckler
 */
public class Chromosome implements Comparable<Chromosome> {

    private static final Pattern DIGITS = Pattern.compile("^\\d+");

    public static Chromosome UNKNOWN = new Chromosome("Unknown");
    private final String myName;
    private final int myChromosomeNumber;
    private final String myCompareString;
    private final int myLength;
    private final GeneralAnnotation myGA;
    private final int hashCode;

    // since there are numerous redundant chromosome, this class use a hash, so that
    // only the pointers are stored.
    private static final ConcurrentMap<Chromosome, Chromosome> CHR_HASH = new ConcurrentHashMap<>(50);

    // this is chromosome cache for when only name is specified
    private static final ConcurrentHashMap<String, Chromosome> CHROMOSOME_NAME_ONLY = new ConcurrentHashMap<>(25);

    public static Chromosome getCanonicalChromosome(Chromosome chr) {
        if (CHR_HASH.size() > 1000) {
            CHR_HASH.clear();
        }
        Chromosome canon = CHR_HASH.putIfAbsent(chr, chr);
        return (canon == null) ? chr : canon;
    }

    /**
     * Creates Chromosome instance with specified name. Returns single instance given same name multiple times.
     *
     * @param name chromosome name
     *
     * @return Chromosome
     */
    public static Chromosome instance(String name) {
        return CHROMOSOME_NAME_ONLY.computeIfAbsent(name, s -> new Chromosome(name));
    }

    /**
     * Creates Chromosome instance with specified name. Returns single instance given same name multiple times.
     *
     * @param name chromosome name
     *
     * @return Chromosome
     */
    public static Chromosome instance(int name) {
        String chr = String.valueOf(name);
        return instance(chr);
    }

    /**
     * @param name Name of the chromosome
     * @param length Length of chromosome in base pairs
     * @param features Map of features about the chromosome
     */
    public Chromosome(String name, int length, GeneralAnnotation features) {
        if (name == null || name.isEmpty()) {
            throw new IllegalArgumentException("Chromosome: name can't be null or empty.");
        }
        myName = parseName(name);
        Matcher matcher = DIGITS.matcher(myName);
        if (matcher.find()) {
            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < 5 - matcher.end(); i++) {
                builder.append("0");
            }
            builder.append(myName);
            myCompareString = builder.toString();
        } else {
            myCompareString = myName;
        }
        myLength = length;
        int convChr = Integer.MAX_VALUE;
        try {
            convChr = Integer.parseInt(myName);
        } catch (NumberFormatException ne) {
            // Use Integer.MAX_VALUE
        }
        myChromosomeNumber = convChr;
        myGA = features;
        hashCode = calcHashCode();
    }

    /**
     * @deprecated use {@link #instance(String)}
     */
    @Deprecated
    public Chromosome(String name) {
        this(name, -1, parseAnnotationFromName(name));
    }

    public String getName() {
        return myName;
    }

    /**
     * Returns the integer value of the chromosome (if name is not a number then
     * Integer.MAX_VALUE is returned)
     */
    public int getChromosomeNumber() {
        return myChromosomeNumber;
    }

    public int getLength() {
        return myLength;
    }

    public GeneralAnnotation getAnnotation() {
        return myGA;
    }

    @Override
    public String toString() {
        return getName();
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    private int calcHashCode() {
        return 79 * 7 + myName.hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) {
            return true;
        }
        if (!(obj instanceof Chromosome)) {
            return false;
        }
        if (hashCode != obj.hashCode()) return false;
        return (compareTo((Chromosome) obj) == 0);
    }

    /**
     * Compares chromosomes numerically if both are numbers.  Otherwise compares as strings.
     *
     * @param o other chromosome
     *
     * @return a negative integer, zero, or a positive integer as this object is less than, equal to, or greater than
     * the specified object.
     */
    @Override
    public int compareTo(Chromosome o) {
        if (this == o) {
            return 0;
        }
        if (myChromosomeNumber != Integer.MAX_VALUE && o.myChromosomeNumber != Integer.MAX_VALUE) {
            return Integer.compare(myChromosomeNumber, o.myChromosomeNumber);
        }
        return myCompareString.compareTo(o.myCompareString);
    }

    /**
     * Takes a string, makes all upper case, removes leading CHROMOSOME/CHR,
     * returns the resulting string
     *
     * @param name name of chromosome
     *
     * @return the input string minus a leading "chr" or "chromsome"
     */
    private static String parseName(String name) {
        String parsedName = name.trim();
        parsedName = parsedName.toUpperCase();
        if (parsedName.startsWith("CHROMOSOME")) {
            parsedName = parsedName.replaceFirst("CHROMOSOME", "");
        }
        if (parsedName.startsWith("CHR")) {
            parsedName = parsedName.replaceFirst("CHR", "");
        }
        int spaceIndex = parsedName.indexOf(" ");
        if (spaceIndex > 0) {
            parsedName = parsedName.substring(0, parsedName.indexOf(" "));
        }
        return parsedName;
    }

    /**
     * Takes a chromosome name, looks for the first space, returns
     * the data beyond as an annotation.  This takes care of lines in
     * a fasta file that look like this:
     * >3 This is a description
     *
     * @param name - the string chromosome passed in
     *
     * @return Annotations built from the string beyond the name
     */
    private static GeneralAnnotation parseAnnotationFromName(String name) {
        GeneralAnnotation annotations = null;
        int spaceIndex = name.indexOf(" ");
        if (spaceIndex > 0) {
            String currChrDesc = name.substring(name.indexOf(" ") + 1);
            annotations = GeneralAnnotationStorage.getBuilder().addAnnotation("Description", currChrDesc).build();
        }

        return annotations;
    }

}
