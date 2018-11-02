/*
 * GenotypeTableUtils
 */
package net.maizegenetics.dna.snp;

import java.io.BufferedReader;

import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

import java.nio.ByteBuffer;
import java.util.*;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;

import static net.maizegenetics.dna.snp.GenotypeTable.*;
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 * Utility methods for comparing, sorting, and counting genotypes.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class GenotypeTableUtils {

    private static final Logger myLogger = Logger.getLogger(GenotypeTableUtils.class);

    private static final Integer ONE = 1;
    private static final byte HIGHMASK = (byte) 0x0F;

    private GenotypeTableUtils() {
        // utility class
    }

    /**
     * This sorts alleles by frequency. Each cell in the given array contains a
     * diploid value which is separated and counted individually. Resulting
     * double dimension array holds alleles (bytes) in result[0]. And the counts
     * are in result[1]. Counts haploid values twice and diploid values once.
     * Higher ploids are not supported.
     *
     * @param data data
     *
     * @return alleles and counts
     */
    public static int[][] getAllelesSortedByFrequency(byte[] data) {

        int[] stateCnt = new int[16];
        for (int i = 0; i < data.length; i++) {
            byte first = (byte) ((data[i] >>> 4) & 0xf);
            byte second = (byte) (data[i] & 0xf);
            if (first < RARE_ALLELE) {
                stateCnt[first]++;
            }
            if (second < RARE_ALLELE) {
                stateCnt[second]++;
            }
        }

        int count = 0;
        for (int j = 0; j < 16; j++) {
            if (stateCnt[j] != 0) {
                count++;
            }
        }

        int result[][] = new int[2][count];
        int index = 0;
        for (int k = 0; k < 16; k++) {
            if (stateCnt[k] != 0) {
                result[0][index] = k;
                result[1][index] = stateCnt[k];
                index++;
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

                if (result[1][k] < result[1][k + 1]) {

                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }

    /**
     * This sorts alleles in a given site by frequency. Each cell in the given
     * array contains a diploid value which is separated and counted
     * individually. Resulting double dimension array holds alleles (bytes) in
     * result[0]. And the counts are in result[1]. Counts haploid values twice
     * and diploid values once. Higher ploids are not supported.
     *
     * @param data data
     * @param site site
     *
     * @return alleles and counts
     */
    public static int[][] getAllelesSortedByFrequency(byte[][] data, int site) {

        int[] stateCnt = new int[16];
        for (int i = 0; i < data.length; i++) {
            byte first = (byte) ((data[i][site] >>> 4) & 0xf);
            byte second = (byte) (data[i][site] & 0xf);
            if (first < RARE_ALLELE) {
                stateCnt[first]++;
            }
            if (second < RARE_ALLELE) {
                stateCnt[second]++;
            }
        }

        int count = 0;
        for (int j = 0; j < 16; j++) {
            if (stateCnt[j] != 0) {
                count++;
            }
        }

        int result[][] = new int[2][count];
        int index = 0;
        for (int k = 0; k < 16; k++) {
            if (stateCnt[k] != 0) {
                result[0][index] = k;
                result[1][index] = stateCnt[k];
                index++;
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

                if (result[1][k] < result[1][k + 1]) {

                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;
    }

    public static byte[] getAlleles(byte[][] data, int site) {
        int[][] alleles = getAllelesSortedByFrequency(data, site);
        int resultSize = alleles[0].length;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i];
        }
        return result;
    }

    public static Object[][] getAllelesSortedByFrequency(String[][] data, int site) {

        Map<String, Integer> stateCnt = new HashMap<>();
        for (int i = 0; i < data.length; i++) {
            String[] temp = data[i][site].split(":");
            String first;
            String second;
            if ((temp == null) || (temp.length == 0)) {
                first = second = UNKNOWN_ALLELE_STR;
            } else if (temp.length == 1) {
                first = second = temp[0].trim();
            } else {
                first = temp[0].trim();
                second = temp[1].trim();
            }
            if (!first.equalsIgnoreCase(UNKNOWN_ALLELE_STR)) {
                Integer count = stateCnt.get(first);
                if (count == null) {
                    stateCnt.put(first, ONE);
                } else {
                    stateCnt.put(first, count + 1);
                }
            }
            if (!second.equalsIgnoreCase(UNKNOWN_ALLELE_STR)) {
                Integer count = stateCnt.get(second);
                if (count == null) {
                    stateCnt.put(second, ONE);
                } else {
                    stateCnt.put(second, count + 1);
                }
            }
        }

        int count = stateCnt.size();

        Object[][] result = new Object[2][count];
        Iterator<String> itr = stateCnt.keySet().iterator();
        int index = 0;
        while (itr.hasNext()) {
            String key = itr.next();
            result[0][index] = key;
            result[1][index] = (Integer) stateCnt.get(key);
            index++;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

                if ((Integer) result[1][k] < (Integer) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;
    }

    /**
     * Converts the byte representation of genotypes to string list of genotypes
     *
     * @param data
     * @return list of the genotypes.
     */
    public static List<String> convertNucleotideGenotypesToStringList(byte[] data) {
        List<String> result = new ArrayList<>();
        for (byte b : data) {
            result.add(NucleotideAlignmentConstants.getHaplotypeNucleotide(b));
        }
        return result;
    }

    public static List<String> getAlleles(String[][] data, int site) {
        Object[][] alleles = getAllelesSortedByFrequency(data, site);
        if ((alleles == null) || (alleles.length == 0)) {
            return null;
        }
        int resultSize = alleles[0].length;
        List<String> result = new ArrayList<>();
        for (int i = 0; i < resultSize; i++) {
            result.add((String) alleles[0][i]);
        }
        return result;
    }

    public static int[][] getAllelesSortedByFrequency(GenotypeTable alignment, int site) {

        int[] stateCnt = new int[16];
        for (int i = 0; i < alignment.numberOfTaxa(); i++) {
            byte[] dipB = alignment.genotypeArray(i, site);
            if (dipB[0] != GenotypeTable.UNKNOWN_ALLELE) {
                stateCnt[dipB[0]]++;
            }
            if (dipB[1] != GenotypeTable.UNKNOWN_ALLELE) {
                stateCnt[dipB[1]]++;
            }
        }

        int count = 0;
        for (int j = 0; j < 16; j++) {
            if (stateCnt[j] != 0) {
                count++;
            }
        }

        int result[][] = new int[2][count];
        int index = 0;
        for (int k = 0; k < 16; k++) {
            if (stateCnt[k] != 0) {
                result[0][index] = k;
                result[1][index] = stateCnt[k];
                index++;
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

                if (result[1][k] < result[1][k + 1]) {

                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }

    public static Object[][] getDiploidsSortedByFrequency(GenotypeTable alignment, int site) {

        Integer ONE_INTEGER = 1;
        int numTaxa = alignment.numberOfTaxa();

        Map<String, Integer> diploidValueCounts = new HashMap<String, Integer>();
        for (int r = 0; r < numTaxa; r++) {
            String current = alignment.genotypeAsString(r, site);
            Integer num = diploidValueCounts.get(current);
            if (num == null) {
                diploidValueCounts.put(current, ONE_INTEGER);
            } else {
                diploidValueCounts.put(current, ++num);
            }
        }

        Object[][] result = new Object[2][diploidValueCounts.size()];

        int i = 0;
        Iterator<String> itr = diploidValueCounts.keySet().iterator();
        while (itr.hasNext()) {
            String key = itr.next();
            Integer count = diploidValueCounts.get(key);
            result[0][i] = key;
            result[1][i++] = count;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0, n = diploidValueCounts.size() - 1; k < n; k++) {

                if ((Integer) result[1][k] < (Integer) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }

    public static String[][] getAlleleStates(String[][] data, int maxNumAlleles) {

        int numSites = data[0].length;

        String[][] alleleStates = new String[numSites][16];
        for (int i = 0; i < numSites; i++) {
            for (int j = 0; j < 16; j++) {
                if (j == RARE_ALLELE) {
                    alleleStates[i][j] = GenotypeTable.RARE_ALLELE_STR;
                } else {
                    alleleStates[i][j] = UNKNOWN_ALLELE_STR;
                }
            }
        }

        for (int site = 0; site < numSites; site++) {
            List<String> alleles = GenotypeTableUtils.getAlleles(data, site);
            if (alleles != null) {
                int numAlleles = Math.min(alleles.size(), maxNumAlleles);
                for (int k = 0; k < numAlleles; k++) {
                    alleleStates[site][k] = (String) alleles.get(k);
                }
            }
        }

        return alleleStates;

    }

    /**
     * remove sites based on minimum frequency (the count of good bases,
     * INCLUDING GAPS) and based on the proportion of good alleles (including
     * gaps) different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the
     * consensus
     * @param minimumCount minimum number of sequences with a good bases (not N
     * or ?), where GAP IS CONSIDERED A GOOD BASE
     */
    public static GenotypeTable removeSitesBasedOnFreqIgnoreMissing(GenotypeTable aa, double minimumProportion, double maximumProportion, int minimumCount) {
        if (!aa.hasGenotype()) {
            return aa;
        }
        int[] includeSites = getIncludedSitesBasedOnFreqIgnoreMissing(aa, minimumProportion, maximumProportion, minimumCount);
        return FilterGenotypeTable.getInstance(aa, includeSites);
    }

    /**
     * This returns the subset of genotypes specified by the given BED file.
     * Start positions are inclusive and end positions are exclusive. If
     * includeSites is false, this returns everything except the subset
     * specified by the BED file.
     *
     * @param input original genotype
     * @param bedFile BED file specifying subset
     * @param includeSites whether to include sites
     *
     * @return subset of genotypes
     */
    public static GenotypeTable filterSitesByBedFile(GenotypeTable input, String bedFile, boolean includeSites) {

        int numSites = input.numberOfSites();
        BitSet sitesToInclude = new OpenBitSet(numSites);
        String line = null;
        try (BufferedReader reader = Utils.getBufferedReader(bedFile)) {
            int lineNum = 1;
            line = reader.readLine();
            while (line != null) {
                String[] tokens = line.trim().split("\t");
                if (tokens.length < 3) {
                    throw new IllegalStateException("filterSitesByBedFile: Expecting at least 3 columns on line: " + lineNum);
                }
                // tokens[1] is start postion from bed file.
                // plus one because bed files are 0-base
                int startSite = input.siteOfPhysicalPosition(Integer.parseInt(tokens[1]) + 1, new Chromosome(tokens[0]));
                if (startSite < 0) {
                    startSite = -startSite - 1;
                }
                // tokens[2] is start postion from bed file.
                // plus one because bed files are 0-base
                int endSite = input.siteOfPhysicalPosition(Integer.parseInt(tokens[2]) + 1, new Chromosome(tokens[0]));
                if (endSite < 0) { // end position doesn't exist, so already excluded
                    endSite = -endSite - 2;
                } else { // end position is exclusive
                    endSite--;
                }
                for (int i = startSite; i <= endSite; i++) {
                    sitesToInclude.fastSet(i);
                }
                line = reader.readLine();
                lineNum++;
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("filterSitesByBedFile: problem reading: " + bedFile + " line: " + line);
        }

        if (!includeSites) {
            sitesToInclude.flip(0, numSites);
        }

        int numNewSites = (int) sitesToInclude.cardinality();
        int[] result = new int[numNewSites];
        int count = 0;
        for (int s = 0; s < numSites; s++) {
            if (sitesToInclude.fastGet(s)) {
                result[count++] = s;
            }
        }

        return FilterGenotypeTable.getInstance(input, result);

    }

    public static GenotypeTable filterSitesByChrPos(GenotypeTable input, PositionList positionList, boolean includeSites) {

        TreeSet<Integer> temp = new TreeSet<>();

        PositionList origPositionList = input.positions();

        Chromosome[] chromosomes = positionList.chromosomes();
        for (Chromosome chromosome : chromosomes) {
            if (origPositionList.chromosomeSiteCount(chromosome) != 0) {
                int[] startEndSubSet = positionList.startAndEndOfChromosome(chromosome);
                int[] startEnd = origPositionList.startAndEndOfChromosome(chromosome);
                int posIndex = startEndSubSet[0];
                for (int site = startEnd[0]; site <= startEnd[1]; site++) {

                    Position current = origPositionList.get(site);
                    int pos = current.getPosition();

                    for (; posIndex <= startEndSubSet[1]; posIndex++) {
                        int currentPos = positionList.get(posIndex).getPosition();
                        if (currentPos == pos) {
                            temp.add(site);
                        } else if (currentPos > pos) {
                            break;
                        }
                    }

                }
            }
        }

        if (temp.size() == origPositionList.size()) {
            return input;
        }

        if (includeSites) {
            return FilterGenotypeTable.getInstance(input, toPrimitive(temp));
        } else {
            int numSites = input.numberOfSites();
            int[] result = new int[numSites - temp.size()];
            int count = 0;
            for (int i = 0; i < numSites; i++) {
                if (!temp.contains(i)) {
                    result[count++] = i;
                }
            }
            return FilterGenotypeTable.getInstance(input, result);
        }

    }

    public static GenotypeTable filterSitesByChrPos(GenotypeTable input, String filename, boolean includeSites) {

        Map<String, List<Integer>> positionsByChr = new HashMap<>();
        int lineNum = 1;
        String[] tokens = null;
        try (BufferedReader reader = Utils.getBufferedReader(filename)) {
            String line = reader.readLine();
            if ((line != null) && (line.toLowerCase().startsWith("chr"))) {
                line = reader.readLine();
                lineNum++;
            }
            while (line != null) {
                tokens = line.split("\t");
                if (tokens.length != 2) {
                    throw new IllegalStateException("GenotypeTableUtils: filterSitesByChrPos: Each line in file must have 2 columns (chr, position). line: " + lineNum + " has: " + tokens.length);
                }
                List<Integer> positions = positionsByChr.get(tokens[0]);
                if (positions == null) {
                    positions = new ArrayList<>();
                    positionsByChr.put(tokens[0], positions);
                }
                positions.add(Integer.valueOf(tokens[1]));
                line = reader.readLine();
                lineNum++;
            }
        } catch (NumberFormatException nfe) {
            throw new IllegalArgumentException("GenotypeTableUtils: filterSitesByChrPos: value: " + tokens[1] + " is not a number on line: " + lineNum);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("GenotypeTableUtils: filterSitesByChrPos: Problem reading file: " + filename + "\n" + e.getMessage());
        }

        return GenotypeTableUtils.filterSitesByChrPos(input, positionsByChr, includeSites);

    }

    private static GenotypeTable filterSitesByChrPos(GenotypeTable input, Map<String, List<Integer>> positionsByChr, boolean includeSites) {
        if (positionsByChr.isEmpty()) {
            return null;
        } else if (positionsByChr.size() == 1) {
            for (Map.Entry<String, List<Integer>> entry : positionsByChr.entrySet()) {
                if (includeSites) {
                    return FilterGenotypeTable.getInstance(input, toPrimitive(getSitesToKeepChrPos(input, entry.getKey(), entry.getValue())));
                } else {
                    throw new UnsupportedOperationException();
                }
            }
            return null; // should never reach this.
        } else {
            List<Integer> allSites = new ArrayList<>();
            for (Map.Entry<String, List<Integer>> entry : positionsByChr.entrySet()) {
                allSites.addAll(getSitesToKeepChrPos(input, entry.getKey(), entry.getValue()));
            }
            if (includeSites) {
                return FilterGenotypeTable.getInstance(input, toPrimitive(allSites));
            } else {
                throw new UnsupportedOperationException();
            }
        }
    }

    public static GenotypeTable keepSitesChrPos(GenotypeTable input, String chromosome, List<Integer> position) {
        return FilterGenotypeTable.getInstance(input, toPrimitive(getSitesToKeepChrPos(input, chromosome, position)));
    }

    private static List<Integer> getSitesToKeepChrPos(GenotypeTable input, String chromosome, List<Integer> position) {

        if ((chromosome == null) || (position == null)) {
            throw new IllegalArgumentException("GenotypeTableUtils: keepSitesChrPos: must specify chromosome and positions.");
        }

        Collections.sort(position);
        PositionList origPositions = input.positions();
        int site = 0;
        for (Position current : origPositions) {
            if (current.getChromosome().getName().equals(chromosome)) {
                break;
            }
            site++;
        }

        int numSites = origPositions.numberOfSites();

        if (site == numSites) {
            return null;
        }

        List<Integer> temp = new ArrayList<>();

        int posIndex = 0;
        for (; site < numSites; site++) {
            Position current = origPositions.get(site);
            if (!current.getChromosome().getName().equals(chromosome)) {
                break;
            }

            int pos = current.getPosition();
            while (position.get(posIndex) < pos) {
                posIndex++;
                if (posIndex == position.size()) {
                    return temp;
                }
            }

            if (position.get(posIndex) == pos) {
                temp.add(site);
            }
        }

        return temp;
    }

    private static int[] toPrimitive(Collection<Integer> list) {
        int[] result = new int[list.size()];
        int index = 0;
        for (int current : list) {
            result[index++] = current;
        }
        return result;
    }

    /**
     * get sites to be included based on minimum frequency (the count of good
     * bases, INCLUDING GAPS) and based on the proportion of good sites
     * (INCLUDING GAPS) different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the
     * consensus
     * @param maximumProportion maximum proportion of sites different from the
     * consensus
     * @param minimumCount minimum number of sequences with a good base or a gap
     * (but not N or ?)
     */
    public static int[] getIncludedSitesBasedOnFreqIgnoreMissing(GenotypeTable aa, double minimumProportion, double maximumProportion, int minimumCount) {

        ArrayList<Integer> includeAL = new ArrayList<>();
        int numSites = aa.numberOfSites();

        if (minimumCount > 0) {
            for (int i = 0; i < numSites; i++) {

                int[][] alleles = aa.allelesSortedByFrequency(i);
                int totalNonMissing = AlleleFreqCache.totalGametesNonMissingForSite(alleles);

                if (totalNonMissing >= (minimumCount * 2)) {

                    double obsMinProp = AlleleFreqCache.minorAlleleFrequency(alleles);

                    if ((obsMinProp >= minimumProportion) && (obsMinProp <= maximumProportion)) {
                        includeAL.add(i);
                    }

                }
            }
        } else {
            for (int i = 0; i < numSites; i++) {

                int[][] alleles = aa.allelesSortedByFrequency(i);
                double obsMinProp = AlleleFreqCache.minorAlleleFrequency(alleles);
                if ((obsMinProp >= minimumProportion) && (obsMinProp <= maximumProportion)) {
                    includeAL.add(i);
                }

            }
        }

        int[] includeSites = new int[includeAL.size()];
        for (int i = 0; i < includeAL.size(); i++) {
            includeSites[i] = includeAL.get(i);
        }
        return includeSites;
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
     * @return true if allele values are the same; false if unknown or not equal
     */
    public static boolean isHomozygous(byte diploidAllele) {
        if (diploidAllele == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
            return false;
        }
        return ((diploidAllele >>> 4) & 0xf) == (diploidAllele & 0xf);
    }

    /**
     * Returns whether two diploid allele values are equal ignoring order.
     *
     * @param alleles1 diploid alleles 1
     * @param alleles2 diploid alleles 2
     *
     * @return true if equal
     */
    public static boolean isEqual(byte[] alleles1, byte[] alleles2) {

        return ((alleles1[0] == alleles2[0]) && (alleles1[1] == alleles2[1]))
                || ((alleles1[0] == alleles2[1]) && (alleles1[1] == alleles2[0]));

    }

    /**
     * Returns whether two diploid allele values are equal ignoring order.
     *
     * @param diploidAllele1 diploid alleles 1
     * @param diploidAllele2 diploid alleles 2
     *
     * @return true if equal
     */
    public static boolean isEqual(byte diploidAllele1, byte diploidAllele2) {

        if (diploidAllele1 != diploidAllele2) {
            byte reversed = (byte) ((diploidAllele1 << 4) | (diploidAllele1 >>> 4));
            if (reversed != diploidAllele2) {
                return false;
            }
        }

        return true;

    }

    /**
     * Returns whether two diploid allele values are equal ignoring order where
     * unknown values equal anything.
     *
     * @param alleles1 diploid alleles 1
     * @param alleles2 diploid alleles 2
     *
     * @return true if equal
     */
    public static boolean isEqualOrUnknown(byte[] alleles1, byte[] alleles2) {

        if (((alleles1[0] == GenotypeTable.UNKNOWN_ALLELE) && (alleles1[1] == GenotypeTable.UNKNOWN_ALLELE))
                || ((alleles2[0] == GenotypeTable.UNKNOWN_ALLELE) && (alleles2[1] == GenotypeTable.UNKNOWN_ALLELE))) {
            return true;
        }

        if (((alleles1[0] == alleles2[0]) && (alleles1[1] == alleles2[1]))
                || ((alleles1[0] == alleles2[1]) && (alleles1[1] == alleles2[0]))) {
            return true;
        } else {
            return false;
        }

    }

    /**
     * Returns whether two diploid allele values are equal ignoring order where
     * unknown values equal anything.
     *
     * @param diploidAllele1 diploid alleles 1
     * @param diploidAllele2 diploid alleles 2
     *
     * @return true if equal
     */
    public static boolean isEqualOrUnknown(byte diploidAllele1, byte diploidAllele2) {

        if ((diploidAllele1 == UNKNOWN_DIPLOID_ALLELE) || (diploidAllele2 == UNKNOWN_DIPLOID_ALLELE)) {
            return true;
        }

        if (diploidAllele1 != diploidAllele2) {
            byte reversed = (byte) ((diploidAllele1 << 4) | (diploidAllele1 >>> 4));
            if (reversed != diploidAllele2) {
                return false;
            }
        }

        return true;

    }

    /**
     * Return true if either at least one allele agree
     *
     * @param genotype1
     * @param genotype2
     * @return true if at least one allele is equal
     */
    public static boolean isPartiallyEqual(byte genotype1, byte genotype2) {
        int low1 = 0xF & genotype1;
        int low2 = 0xF & genotype2;
        if (low1 == low2) {
            return true;
        }
        int high1 = genotype1 >>> 4;
        if (high1 == low2) {
            return true;
        }
        int high2 = genotype2 >>> 4;
        if (low1 == high2) {
            return true;
        }
        if (high1 == high2) {
            return true;
        }
        return false;
    }

    public static boolean areEncodingsEqual(String[][][] encodings) {
        int numEncodings = encodings.length;
        for (int i = 1; i < numEncodings; i++) {
            int numSites = encodings[0].length;
            if (numSites != encodings[i].length) {
                return false;
            }
            for (int s = 0; s < numSites; s++) {
                int numCodes = encodings[0][s].length;
                if (numCodes != encodings[i][s].length) {
                    return false;
                }
                for (int c = 0; c < numCodes; c++) {
                    if (!encodings[0][s][c].equals(encodings[i][s][c])) {
                        return false;
                    }
                }
            }
        }
        return true;
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
     * Ensures diploid value in alphabetical order
     *
     * @param genotype diploid genotype
     *
     * @return diploid value sorted by order A < C < G < T
     */
    public static byte getUnphasedSortedDiploidValue(byte genotype) {
        byte a = (byte) ((genotype >>> 4) & 0xf);
        byte b = (byte) (genotype & 0xf);
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
     * @return diploid value
     */
    public static byte getUnphasedDiploidValueNoHets(byte g1, byte g2) {
        if ((g2 == g1) && (!isHeterozygous(g1))) {
            return g1;
        }
        if (g1 == UNKNOWN_DIPLOID_ALLELE) {
            return UNKNOWN_DIPLOID_ALLELE;
        }
        if (g2 == UNKNOWN_DIPLOID_ALLELE) {
            return UNKNOWN_DIPLOID_ALLELE;
        }
        if (isHeterozygous(g1)) {
            return UNKNOWN_DIPLOID_ALLELE;
        }
        if (isHeterozygous(g2)) {
            return UNKNOWN_DIPLOID_ALLELE;
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
        return result;
    }

    /**
     * Method for getting TBits rapidly from major and minor allele arrays
     *
     * @param genotype
     * @param mjA
     * @param mnA
     * @return
     */
    public static BitSet[] calcBitPresenceFromGenotype(byte[] genotype, byte[] mjA, byte[] mnA) {
        int sites = genotype.length;
        if ((genotype.length != mjA.length) || (genotype.length != mnA.length)) {
            throw new ArrayIndexOutOfBoundsException("Input genotypes unequal in length");
        }
        OpenBitSet rMj = new OpenBitSet(genotype.length);
        OpenBitSet rMn = new OpenBitSet(genotype.length);
        for (int i = 0; i < sites; i++) {
            byte g = genotype[i];
            byte mj = mjA[i];
            byte mn = mnA[i];
            //           System.out.printf("inc:%d g:%d mj:%d mn:%d %n", i, g, mj, mn);
            if (mj == UNKNOWN_ALLELE) {
                continue;
            }
            if (g == getDiploidValuePhased(mj, mj)) {
                rMj.fastSet(i);
                continue;
            }
            if (mn == UNKNOWN_ALLELE) {
                continue;
            }
            if (g == getDiploidValuePhased(mn, mn)) {
                rMn.fastSet(i);
                continue;
            }
            byte het = getUnphasedDiploidValue(mj, mn);
            if (isEqual(g, het)) {
                rMj.fastSet(i);
                rMn.fastSet(i);
            }
        }
        return new BitSet[]{rMj, rMn};
    }

    /**
     * Returns BitSet indicating presence of given diploid value in genotype.
     * This does unphased comparisons, so the order of the two allele values do
     * not matter. If given diploid value is UNKNOWN, then it doesn't match
     * anything. Bits set to 1 indicate match.
     *
     * @param genotype genotype
     * @param diploidValue diploid value
     *
     * @return BitSet indicating presence
     */
    public static BitSet calcBitPresenceOfDiploidValueFromGenotype(byte[] genotype, byte diploidValue) {
        int sites = genotype.length;
        OpenBitSet result = new OpenBitSet(genotype.length);
        for (int i = 0; i < sites; i++) {

            if (diploidValue == UNKNOWN_DIPLOID_ALLELE) {
                // do nothing
            } else if (isEqual(diploidValue, genotype[i])) {
                result.fastSet(i);
            }

        }
        return result;
    }

    public static BitSet calcBitPresenceFromGenotype(byte[] genotype, byte[] referenceValues) {
        int sites = genotype.length;
        if (genotype.length != referenceValues.length) {
            throw new ArrayIndexOutOfBoundsException("GenotypeTableUtils: calcBitPresenceFromGenotype: Input genotypes unequal in length");
        }
        OpenBitSet result = new OpenBitSet(genotype.length);
        for (int i = 0; i < sites; i++) {

            if (referenceValues[i] == UNKNOWN_ALLELE) {
                // do nothing
            } else if (referenceValues[i] == (byte) (genotype[i] & 0xf)) {
                result.fastSet(i);
            } else if (referenceValues[i] == (byte) ((genotype[i] >>> 4) & 0xf)) {
                result.fastSet(i);
            }

        }
        return result;
    }

    public static BitSet calcBitUnknownPresenceFromGenotype(byte[] genotype) {
        int length = genotype.length;
        OpenBitSet result = new OpenBitSet(genotype.length);
        for (int i = 0; i < length; i++) {

            if (genotype[i] == UNKNOWN_DIPLOID_ALLELE) {
                result.fastSet(i);
            }

        }
        return result;
    }

    /**
     * Method for getting TBits rapidly from major and minor allele arrays
     *
     * @param genotype
     * @param mjA
     * @param mnA
     * @return
     */
    public static BitSet[] calcBitPresenceFromGenotype15(byte[] genotype, byte[] mjA, byte[] mnA) {
        int sites = genotype.length;
        if ((genotype.length != mjA.length) || (genotype.length != mnA.length)) {
            throw new ArrayIndexOutOfBoundsException("Input genotypes unequal in length");
        }
        ByteBuffer gBB = ByteBuffer.wrap(genotype);  //byte buffer make the code 20% faster for short sequences
        ByteBuffer mjBB = ByteBuffer.wrap(mjA);
        ByteBuffer mnBB = ByteBuffer.wrap(mnA);
        OpenBitSet rMj = new OpenBitSet(genotype.length);
        OpenBitSet rMn = new OpenBitSet(genotype.length);
        for (int i = 0; i < sites; i++) {
            byte g = gBB.get();
            byte mj = mjBB.get();
            byte mn = mnBB.get();
            // System.out.printf("inc:%d g:%d mj:%d mn:%d %n", i, g, mj, mn);
            if (mj == GenotypeTable.UNKNOWN_ALLELE) {
                continue;
            }
            if (g == GenotypeTableUtils.getDiploidValuePhased(mj, mj)) {
                rMj.fastSet(i);
                continue;
            }
            if (mn == GenotypeTable.UNKNOWN_ALLELE) {
                continue;
            }
            if (g == GenotypeTableUtils.getDiploidValuePhased(mn, mn)) {
                rMn.fastSet(i);
                continue;
            }
            byte het = GenotypeTableUtils.getUnphasedDiploidValue(mj, mn);
            if (GenotypeTableUtils.isEqual(g, het)) {
                rMj.fastSet(i);
                rMn.fastSet(i);
            }
        }
        return new BitSet[]{rMj, rMn};
    }

    /**
     * Method for getting Site Bits rapidly from major and minor alleles
     *
     * @param genotype
     * @param mj
     * @param mn
     * @return
     */
    public static BitSet[] calcBitPresenceFromGenotype(byte[] genotype, byte mj, byte mn) {
        int sites = genotype.length;
        OpenBitSet rMj = new OpenBitSet(genotype.length);
        OpenBitSet rMn = new OpenBitSet(genotype.length);
        for (int i = 0; i < sites; i++) {
            byte g = genotype[i];
            if (mj == UNKNOWN_ALLELE) {
                continue;
            }
            if (g == GenotypeTableUtils.getDiploidValuePhased(mj, mj)) {
                rMj.fastSet(i);
                continue;
            }
            if (mn == UNKNOWN_ALLELE) {
                continue;
            }
            if (g == GenotypeTableUtils.getDiploidValuePhased(mn, mn)) {
                rMn.fastSet(i);
                continue;
            }
            byte het = GenotypeTableUtils.getUnphasedDiploidValue(mj, mn);
            if (GenotypeTableUtils.isEqual(g, het)) {
                rMj.fastSet(i);
                rMn.fastSet(i);
            }
        }
        return new BitSet[]{rMj, rMn};
    }

    public static BitSet calcBitPresenceFromGenotype(byte[] genotype, byte referenceValue) {
        int length = genotype.length;
        OpenBitSet result = new OpenBitSet(genotype.length);
        if (referenceValue == UNKNOWN_ALLELE) {
            return result;
        }
        for (int i = 0; i < length; i++) {

            if (referenceValue == (byte) (genotype[i] & 0xf)) {
                result.fastSet(i);
            } else if (referenceValue == (byte) ((genotype[i] >>> 4) & 0xf)) {
                result.fastSet(i);
            }

        }
        return result;
    }

    public static float[][] convertGenotypeToFloatProbability(GenotypeTable genotype, boolean sitesByTaxa) {
        int nsites = genotype.numberOfSites();
        int ntaxa = genotype.numberOfTaxa();
        float[][] out = null;
        if (sitesByTaxa) {
            out = new float[nsites][ntaxa];
            for (int s = 0; s < nsites; s++) {
                byte major = genotype.majorAllele(s);
                for (int t = 0; t < ntaxa; t++) {
                    byte geno = genotype.genotype(t, s);
                    if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                        out[s][t] = Float.NaN;
                    }
                    byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
                    if (alleles[0] == major) {
                        out[s][t] += 0.5;
                    }
                    if (alleles[1] == major) {
                        out[s][t] += 0.5;
                    }
                }
            }
        }
        if (!sitesByTaxa) {
            out = new float[ntaxa][nsites];
            for (int s = 0; s < nsites; s++) {
                byte major = genotype.majorAllele(s);
                for (int t = 0; t < ntaxa; t++) {
                    byte geno = genotype.genotype(t, s);
                    if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                        out[t][s] = Float.NaN;
                    }
                    byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
                    if (alleles[0] == major) {
                        out[t][s] += 0.5;
                    }
                    if (alleles[1] == major) {
                        out[t][s] += 0.5;
                    }
                }
            }
        }
        return out;

    }

    public static double[][] convertGenotypeToDoubleProbability(GenotypeTable genotype, boolean sitesByTaxa) {
        int nsites = genotype.numberOfSites();
        int ntaxa = genotype.numberOfTaxa();
        double[][] out = null;
        if (sitesByTaxa) {
            out = new double[nsites][ntaxa];
            for (int s = 0; s < nsites; s++) {
                byte major = genotype.majorAllele(s);
                for (int t = 0; t < ntaxa; t++) {
                    byte geno = genotype.genotype(t, s);
                    if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                        out[s][t] = Double.NaN;
                    }
                    byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
                    if (alleles[0] == major) {
                        out[s][t] += 0.5;
                    }
                    if (alleles[1] == major) {
                        out[s][t] += 0.5;
                    }
                }
            }
        }
        if (!sitesByTaxa) {
            out = new double[ntaxa][nsites];
            for (int s = 0; s < nsites; s++) {
                byte major = genotype.majorAllele(s);
                for (int t = 0; t < ntaxa; t++) {
                    byte geno = genotype.genotype(t, s);
                    if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                        out[t][s] = Double.NaN;
                    }
                    byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
                    if (alleles[0] == major) {
                        out[t][s] += 0.5;
                    }
                    if (alleles[1] == major) {
                        out[t][s] += 0.5;
                    }
                }
            }
        }
        return out;
    }

}
