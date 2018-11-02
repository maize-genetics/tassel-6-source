/*
 *  AlleleFreqCache
 */
package net.maizegenetics.dna.snp.genotypecall;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import java.util.concurrent.CopyOnWriteArraySet;
import java.util.concurrent.ForkJoinPool;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;

/**
 * Cache for allele frequency statistics. Allele frequency can be expensive to
 * recalculate at large scale, so these class efficiently loops through blocks
 * of sites and caches the statistics on them.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class AlleleFreqCache {

    private static final int DEFAULT_MAX_NUM_ALLELES = 6;
    private static final int SHIFT_AMOUNT = 8;
    private static final int NUM_SITES_TO_CACHE = 1 << SHIFT_AMOUNT;
    public static final int SITE_BLOCK_MASK = ~(NUM_SITES_TO_CACHE - 1);
    private final GenotypeCallTable myGenotype;
    private final int myMaxNumAlleles;
    private final Cache<Integer, int[][][]> myCachedAlleleFreqs;
    private final ForkJoinPool myThreadPool;

    public AlleleFreqCache(GenotypeCallTable genotype, int maxNumAlleles) {
        myGenotype = genotype;
        myMaxNumAlleles = maxNumAlleles;
        myCachedAlleleFreqs = CacheBuilder.newBuilder()
                .initialCapacity(myGenotype.numberOfSites() / NUM_SITES_TO_CACHE)
                .build();
        myThreadPool = ForkJoinPool.commonPool();
    }

    private static int getStartSite(int site) {
        return site & SITE_BLOCK_MASK;
    }

    private int[][] getCachedAlleleFreq(int site) {
        int startSite = getStartSite(site);
        int[][][] result = myCachedAlleleFreqs.getIfPresent(startSite);
        if (result == null) {
            startLookAhead(startSite);
            if (site == startSite) {
                startLookAhead(startSite + NUM_SITES_TO_CACHE);
                startLookAhead(startSite + NUM_SITES_TO_CACHE * 2);
                startLookAhead(startSite + NUM_SITES_TO_CACHE * 3);
            }
            return alleleFreq(site);
        }
        return result[site - startSite];
    }

    public int[][] getAllelesSortedByFrequency(int site) {
        return getCachedAlleleFreq(site);
    }

    private int[][][] calculateAlleleFreq(int site) {

        int startSite = getStartSite(site);
        int numSites = Math.min(NUM_SITES_TO_CACHE, myGenotype.numberOfSites() - startSite);
        int[][][] alleleCounts = new int[numSites][][];
        for (int s = 0; s < numSites; s++) {
            alleleCounts[s] = alleleFreq(startSite + s);
        }

        myCachedAlleleFreqs.put(startSite, alleleCounts);
        return alleleCounts;
    }

    private int[][] alleleFreq(int site) {
        return alleleFreq(myGenotype.genotypeForAllTaxa(site), myMaxNumAlleles);
    }

    public static int[][] allelesSortedByFrequencyNucleotide(byte[] data) {
        return alleleFreq(data, DEFAULT_MAX_NUM_ALLELES);
    }

    private static int[][] alleleFreq(byte[] data, int maxNumAlleles) {
        int numTaxa = data.length;
        int[] alleleFreq = new int[maxNumAlleles];
        for (int taxon = 0; taxon < numTaxa; taxon++) {
            if (((data[taxon] >>> 4) & 0xf) < maxNumAlleles) {
                alleleFreq[((data[taxon] >>> 4) & 0xf)]++;
            }
            if ((data[taxon] & 0xf) < maxNumAlleles) {
                alleleFreq[(data[taxon] & 0xf)]++;
            }
        }

        for (byte i = 0; i < maxNumAlleles; i++) {
            // size | allele (the 5-i is to get the sort right, so if case of ties A is first)
            alleleFreq[i] = (alleleFreq[i] << 4) | (maxNumAlleles - 1 - i);
        }

        int numAlleles = sort(alleleFreq, maxNumAlleles);
        int[][] alleleCounts = new int[2][numAlleles];
        for (int i = 0; i < numAlleles; i++) {
            alleleCounts[0][i] = (byte) (5 - (0xF & alleleFreq[i]));
            alleleCounts[1][i] = alleleFreq[i] >>> 4;
        }

        return alleleCounts;
    }

    public static Stats allelesSortedByFrequencyAndCountsNucleotide(int index, byte[] data) {
        return alleleFreqCounts(index, data, DEFAULT_MAX_NUM_ALLELES);
    }

    public static final int UNKNOWN_COUNT = 0;
    public static final int UNKNOWN_GAMETE_COUNT = 1;
    public static final int HETEROZYGOUS_COUNT = 2;
    public static final int HOMOZYGOUS_COUNT = 3;

    private static final byte UNKNOWN_COUNT_BIT = 0x1;
    private static final byte UNKNOWN_SINGLE_GAMETE_COUNT_BIT = 0x2;
    private static final byte HETEROZYGOUS_COUNT_BIT = 0x4;
    private static final byte HOMOZYGOUS_COUNT_BIT = 0x8;

    private static final byte[] OTHER_COUNTS = new byte[256];

    static {
        for (int i = 0; i < 256; i++) {

            if (((byte) i) == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                OTHER_COUNTS[i] = UNKNOWN_COUNT_BIT;
                continue;
            }

            byte[] alleles = GenotypeTableUtils.getDiploidValues((byte) i);

            if (alleles[0] == GenotypeTable.UNKNOWN_ALLELE) {
                OTHER_COUNTS[i] = UNKNOWN_SINGLE_GAMETE_COUNT_BIT;
            } else if (alleles[1] == GenotypeTable.UNKNOWN_ALLELE) {
                OTHER_COUNTS[i] = UNKNOWN_SINGLE_GAMETE_COUNT_BIT;
            } else if (alleles[0] == alleles[1]) {
                OTHER_COUNTS[i] = HOMOZYGOUS_COUNT_BIT;
            } else {
                OTHER_COUNTS[i] = HETEROZYGOUS_COUNT_BIT;
            }
        }
    }

    private static Stats alleleFreqCounts(int index, byte[] data, int maxNumAlleles) {

        int numGenotypes = data.length;
        int[] alleleFreq = new int[maxNumAlleles];
        int[] otherCounts = new int[4];
        for (int i = 0; i < numGenotypes; i++) {

            int count = Byte.toUnsignedInt(data[i]);
            otherCounts[UNKNOWN_COUNT] += OTHER_COUNTS[count] & UNKNOWN_COUNT_BIT;
            // << 1 is times 2
            otherCounts[UNKNOWN_GAMETE_COUNT] += (OTHER_COUNTS[count] & UNKNOWN_COUNT_BIT) << 1;
            // this is zero if both gametes where Unknown
            otherCounts[UNKNOWN_GAMETE_COUNT] += (OTHER_COUNTS[count] & UNKNOWN_SINGLE_GAMETE_COUNT_BIT) >>> 1;
            otherCounts[HETEROZYGOUS_COUNT] += (OTHER_COUNTS[count] & HETEROZYGOUS_COUNT_BIT) >>> 2;
            otherCounts[HOMOZYGOUS_COUNT] += (OTHER_COUNTS[count] & HOMOZYGOUS_COUNT_BIT) >>> 3;

            if (((data[i] >>> 4) & 0xf) < maxNumAlleles) {
                alleleFreq[((data[i] >>> 4) & 0xf)]++;
            }
            if ((data[i] & 0xf) < maxNumAlleles) {
                alleleFreq[(data[i] & 0xf)]++;
            }

        }

        for (byte i = 0; i < maxNumAlleles; i++) {
            // size | allele (the 5-i is to get the sort right, so if case of ties A is first)
            alleleFreq[i] = (alleleFreq[i] << 4) | (maxNumAlleles - 1 - i);
        }

        int numAlleles = sort(alleleFreq, maxNumAlleles);
        int[][] alleleCounts = new int[2][numAlleles];
        for (int i = 0; i < numAlleles; i++) {
            alleleCounts[0][i] = (byte) (5 - (0xF & alleleFreq[i]));
            alleleCounts[1][i] = alleleFreq[i] >>> 4;
        }

        return Stats.getInstance(alleleCounts, otherCounts, numGenotypes, index);

    }

    private static int sort(int[] data, int maxNumAlleles) {
        int countNotZero = 0;
        for (int j = 0; j < maxNumAlleles - 1; j++) {
            int imax = j;
            for (int i = j + 1; i < maxNumAlleles; i++) {
                if (data[i] > data[imax]) {
                    imax = i;
                }
            }
            if (data[imax] > 0xF) {
                int temp = data[j];
                data[j] = data[imax];
                data[imax] = temp;
                countNotZero++;
            } else {
                return countNotZero;
            }
        }
        if (data[5] > 0xF) {
            countNotZero++;
        }
        return countNotZero;
    }

    public static byte majorAllele(int[][] alleles) {
        if (alleles[0].length >= 1) {
            return (byte) alleles[0][0];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    public static double majorAlleleFrequency(int[][] alleles) {

        int numAlleles = alleles[0].length;
        if (numAlleles >= 1) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing += alleles[1][i];
            }
            return (double) alleles[1][0] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    public static byte minorAllele(int[][] alleles) {
        if (alleles[0].length >= 2) {
            return (byte) alleles[0][1];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    public static double minorAlleleFrequency(int[][] alleles) {

        int numAlleles = alleles[0].length;
        if (numAlleles >= 2) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing += alleles[1][i];
            }
            return (double) alleles[1][1] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    public static int totalGametesNonMissingForSite(int[][] alleles) {
        int numAlleles = alleles[0].length;
        int result = 0;
        for (int i = 0; i < numAlleles; i++) {
            result += alleles[1][i];
        }
        return result;
    }

    public static double proportionHeterozygous(int[] counts, int totalCount) {
        return counts[HETEROZYGOUS_COUNT] / (double) totalCount;
    }

    public static double proportionHeterozygous(byte[] data) {
        int numNotMissing = 0;
        int numHeterozygous = 0;
        for (byte current : data) {
            if (current != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                numNotMissing++;
                if (GenotypeTableUtils.isHeterozygous(current)) {
                    numHeterozygous++;
                }
            }
        }
        if (numNotMissing == 0) {
            return 0.0;
        } else {
            return (double) numHeterozygous / (double) numNotMissing;
        }
    }

    private final CopyOnWriteArraySet<Integer> myCurrentlyProcessingBlocks = new CopyOnWriteArraySet<>();

    private void startLookAhead(int site) {
        int startSite = getStartSite(site);
        if ((startSite < myGenotype.numberOfSites()) && (myCachedAlleleFreqs.getIfPresent(startSite) == null)) {
            if (myCurrentlyProcessingBlocks.add(startSite)) {
                myThreadPool.execute(new LookAheadSiteStats(startSite));
            }
        }
    }

    private class LookAheadSiteStats implements Runnable {

        private final int myStartSite;

        public LookAheadSiteStats(int site) {
            myStartSite = getStartSite(site);
        }

        @Override
        public void run() {
            try {
                calculateAlleleFreq(myStartSite);
            } finally {
                myCurrentlyProcessingBlocks.remove(myStartSite);
            }
        }
    }
}
