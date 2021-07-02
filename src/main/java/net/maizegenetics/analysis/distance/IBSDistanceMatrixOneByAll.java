/*
 *  IBSDistanceMatrixOneByAll
 * 
 *  Created on Dec 17, 2016
 */
package net.maizegenetics.analysis.distance;

import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache;
import net.maizegenetics.util.ProgressListener;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class IBSDistanceMatrixOneByAll {

    private static final Logger myLogger = LogManager.getLogger(IBSDistanceMatrixOneByAll.class);

    private IBSDistanceMatrixOneByAll() {
        // utility
    }

    public static double[] getInstance(GenotypeTable genotype, int taxonIndex) {
        return getInstance(genotype, taxonIndex, 0, null);
    }

    public static double[] getInstance(GenotypeTable genotype, int taxonIndex, ProgressListener listener) {
        return getInstance(genotype, taxonIndex, 0, listener);
    }

    public static double[] getInstance(GenotypeTable genotype, int taxonIndex, int minSiteComp, ProgressListener listener) {
        return computeHetBitDistances(genotype, taxonIndex, listener, minSiteComp);
    }

    private static double[] computeHetBitDistances(GenotypeTable genotype, int taxonIndex, ProgressListener listener, int minSitesComp) {

        int numSeqs = genotype.numberOfTaxa();
        long time = System.currentTimeMillis();

        Counters temp = new Counters(numSeqs);
        stream(genotype, taxonIndex, listener).forEach((long[] t) -> {
            temp.add(t);
        });

        int[] counters = temp.myCounters;

        double[] result = new double[numSeqs];
        int index = 0;
        for (int i = 0; i < numSeqs; i++) {
            int sameCount = counters[index++];
            int diffCount = counters[index++];
            int hetCount = counters[index++];
            long sites = sameCount + diffCount - hetCount;
            double identity = ((double) (sameCount) - 0.5 * hetCount) / (double) (sites);
            double dist = 1 - identity;

            if (sites < minSitesComp) {
                dist = Double.NaN;
            }
            result[i] = dist;
        }

        myLogger.info("IBSDistanceMatrixOneByAll: computeHetBitDistances time = " + (System.currentTimeMillis() - time) / 1000 + " seconds");

        return result;

    }

    protected static void fireProgress(int percent, ProgressListener listener) {
        if (listener != null) {
            listener.progress(percent, null);
        }

    }

    //
    // Each CPU thread (process) creates an instance of this class
    // to acculate counters of the IBS Distance Matrix. The add()
    // method parses out the three counts from each long that's
    // coming from the stream. These are
    // combined with addAll() to result in one instance at the end.
    // Each three consecutive int holds the same, different, and het
    // count for a pair-wise comparison.
    //
    private static class Counters {

        private final int[] myCounters;
        private final int myNumTaxa;

        public Counters(int numTaxa) {
            myNumTaxa = numTaxa;
            myCounters = new int[myNumTaxa * 3];
        }

        public synchronized void add(long[] values) {
            int index = 0;
            for (int i = 0; i < myNumTaxa; i++) {
                myCounters[index++] += (int) (values[i] & 0x1FFFFFl);
                myCounters[index++] += (int) ((values[i] >>> 21) & 0x1FFFFFl);
                myCounters[index++] += (int) ((values[i] >>> 42) & 0x1FFFFFl);
            }
        }

        public void addAll(Counters counters) {
            int[] other = counters.myCounters;
            for (int t = 0; t < myNumTaxa; t++) {
                myCounters[t] += other[t];
            }
        }

    }

    //
    // These constants named whether a pair-wise comparison is
    // SAME_DIFFERENT.  The value has a 1 in the appropriate
    // 3 x 20 bits depending whether same, different, or het
    //
    private static final long TRUE_TRUE_LONG = 0x40000200001l;
    private static final long TRUE_FALSE_LONG = 0x1l;
    private static final long FALSE_TRUE_LONG = 0x200000l;
    private static final long FALSE_FALSE_LONG = 0x0l;

    //
    // This precalculates the counts for every combination
    // of three sites.
    //
    private static long[] PRECALCULATED_COUNTS = null;

    static {

        long[] possibleTerms = new long[32];
        possibleTerms[22] = TRUE_FALSE_LONG;
        possibleTerms[20] = TRUE_TRUE_LONG;
        possibleTerms[18] = TRUE_TRUE_LONG;
        possibleTerms[6] = FALSE_TRUE_LONG;
        possibleTerms[2] = FALSE_TRUE_LONG;
        possibleTerms[0] = FALSE_FALSE_LONG;
        possibleTerms[21] = TRUE_TRUE_LONG;
        possibleTerms[17] = TRUE_TRUE_LONG;
        possibleTerms[4] = TRUE_TRUE_LONG;
        possibleTerms[1] = TRUE_TRUE_LONG;
        possibleTerms[5] = FALSE_TRUE_LONG;
        possibleTerms[19] = TRUE_TRUE_LONG;
        possibleTerms[3] = TRUE_TRUE_LONG;
        possibleTerms[14] = TRUE_FALSE_LONG;
        possibleTerms[10] = TRUE_TRUE_LONG;
        possibleTerms[11] = TRUE_TRUE_LONG;
        possibleTerms[7] = TRUE_FALSE_LONG;

        PRECALCULATED_COUNTS = new long[23255];

        for (int i = 0; i < 23255; i++) {
            int firstCode = i & 0x1f;
            int secondCode = (i >>> 5) & 0x1f;
            int thirdCode = (i >>> 10) & 0x1f;
            PRECALCULATED_COUNTS[i] = possibleTerms[firstCode] + possibleTerms[secondCode] + possibleTerms[thirdCode];
        }

    }

    //
    // This defines the codes for each possible state at a given
    // site and taxon.
    //
    private static final byte[] PRECALCULATED_ENCODINGS = new byte[8];

    static {
        // 22, 21, 19, 14, 11, 7, 0
        PRECALCULATED_ENCODINGS[1] = 0x16; // Major
        PRECALCULATED_ENCODINGS[3] = 0x15; // Major and Minor
        PRECALCULATED_ENCODINGS[5] = 0x13; // Major and Second Minor
        PRECALCULATED_ENCODINGS[2] = 0xE; // Minor
        PRECALCULATED_ENCODINGS[6] = 0xB; // Minor and Second Minor
        PRECALCULATED_ENCODINGS[4] = 0x7; // Second Minor
        PRECALCULATED_ENCODINGS[0] = 0x0; // Unknown
    }

    private static final int NUM_CORES_TO_USE = Runtime.getRuntime().availableProcessors() - 1;

    //
    // Used to report progress.  This is not thread-safe but
    // works well enough for this purpose.
    //
    private static int myNumSitesProcessed = 0;

    private static final int MAX_NUMBER_20_BITS = 0xFFFFF;

    //
    // Creates stream from IBSSiteSpliterator and Genotype Table
    //
    private static Stream<long[]> stream(GenotypeTable genotypes, int taxonIndex, ProgressListener listener) {
        myNumSitesProcessed = 0;
        return StreamSupport.stream(new IBSSiteSpliterator(genotypes, taxonIndex, 0, genotypes.numberOfSites(), listener), true);
    }

    //
    // Spliterator that splits the sites into halves each time for
    // processing.
    //
    static class IBSSiteSpliterator implements Spliterator<long[]> {

        private int myCurrentSite;
        private final int myFence;
        private final GenotypeTable myGenotypes;
        private final int myNumTaxa;
        private final int myNumSites;
        private final ProgressListener myProgressListener;
        private int myMinSitesToProcess;
        private final int myTaxonIndex;

        IBSSiteSpliterator(GenotypeTable genotypes, int taxonIndex, int currentIndex, int fence, ProgressListener listener) {
            myGenotypes = genotypes;
            myNumTaxa = myGenotypes.numberOfTaxa();
            myNumSites = myGenotypes.numberOfSites();
            myCurrentSite = currentIndex;
            myFence = fence;
            myProgressListener = listener;
            myMinSitesToProcess = myNumSites / NUM_CORES_TO_USE;
            if (myMinSitesToProcess == 0) {
                myMinSitesToProcess = myNumSites;
            }
            myTaxonIndex = taxonIndex;
        }

        @Override
        public void forEachRemaining(Consumer<? super long[]> action) {

            int numSitesProcessed = myFence - myCurrentSite;

            //
            // This prevents overrunning the max number that can
            // be held in 20 bits of the long.
            //
            for (; myCurrentSite < myFence;) {

                int currentBlockFence = Math.min(myCurrentSite + MAX_NUMBER_20_BITS, myFence);
                long[] counts = new long[myNumTaxa];

                for (; myCurrentSite < currentBlockFence;) {

                    int[] numSites = new int[1];

                    //
                    // Gets encodings for several blocks of sites.
                    //
                    short[] encodings1 = getBlockOfSites(myCurrentSite, numSites, currentBlockFence);

                    short[] encodings2 = getBlockOfSites(myCurrentSite + numSites[0], numSites, currentBlockFence);

                    short[] encodings3 = getBlockOfSites(myCurrentSite + numSites[0], numSites, currentBlockFence);

                    short[] encodings4 = getBlockOfSites(myCurrentSite + numSites[0], numSites, currentBlockFence);

                    short[] encodings5 = getBlockOfSites(myCurrentSite + numSites[0], numSites, currentBlockFence);

                    short[] encodings6 = getBlockOfSites(myCurrentSite + numSites[0], numSites, currentBlockFence);

                    myCurrentSite += numSites[0];

                    //
                    // Iterates through all pair-wise combinations of taxa
                    //
                    for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                        //
                        // Combine first taxon's encoding with
                        // second taxon's encoding to
                        // create index into pre-calculated counts
                        //
                        counts[firstTaxa] += PRECALCULATED_COUNTS[encodings1[firstTaxa] & encodings1[myTaxonIndex]]
                                + PRECALCULATED_COUNTS[encodings2[firstTaxa] & encodings2[myTaxonIndex]]
                                + PRECALCULATED_COUNTS[encodings3[firstTaxa] & encodings3[myTaxonIndex]]
                                + PRECALCULATED_COUNTS[encodings4[firstTaxa] & encodings4[myTaxonIndex]]
                                + PRECALCULATED_COUNTS[encodings5[firstTaxa] & encodings5[myTaxonIndex]]
                                + PRECALCULATED_COUNTS[encodings6[firstTaxa] & encodings6[myTaxonIndex]];
                    }
                }

                action.accept(counts);
            }
            myNumSitesProcessed += numSitesProcessed;
            fireProgress((int) ((double) myNumSitesProcessed / (double) myNumSites * 100.0), myProgressListener);
        }

        private static final int NUM_SITES_PER_BLOCK = 3;

        private short[] getBlockOfSites(int currentSite, int[] numSites, int currentBlockFence) {

            int currentSiteNum = 0;

            //
            // This holds the encoding for every taxa.  Each
            // short has encodings in each 5 bits.
            //
            short[] encodings = new short[myNumTaxa];

            while ((currentSiteNum < NUM_SITES_PER_BLOCK) && (currentSite < currentBlockFence)) {

                byte[] genotype = myGenotypes.genotypeAllTaxa(currentSite);
                int[][] alleles = AlleleFreqCache.allelesSortedByFrequencyNucleotide(genotype);
                int numAlleles = alleles[0].length;

                //
                // If whole site is Unknown, then skip the site.
                //
                if (numAlleles != 0) {

                    //
                    // Records presence of major, minor, and second minor alleles
                    // for current site in 5 bits.
                    //
                    for (int i = 0; i < myNumTaxa; i++) {
                        byte first = (byte) (genotype[i] & 0xf);
                        byte second = (byte) (genotype[i] >>> 4 & 0xf);
                        int allelePresent = 0;
                        if ((alleles[0][0] == first) || (alleles[0][0] == second)) {
                            allelePresent = 0x1;
                        }
                        if (numAlleles >= 2) {
                            if ((alleles[0][1] == first) || (alleles[0][1] == second)) {
                                allelePresent |= 0x2;
                            }
                            if (numAlleles >= 3) {
                                if ((alleles[0][2] == first) || (alleles[0][2] == second)) {
                                    allelePresent |= 0x4;
                                }
                            }
                        }

                        encodings[i] = (short) (encodings[i] << 5 | PRECALCULATED_ENCODINGS[allelePresent]);
                    }

                    currentSiteNum++;
                }

                currentSite++;
                numSites[0]++;
            }

            return encodings;

        }

        @Override
        public boolean tryAdvance(Consumer<? super long[]> action) {
            if (myCurrentSite < myFence) {
                forEachRemaining(action);
                return true;
            } else {
                return false;
            }
        }

        @Override
        /**
         * Splits sites
         */
        public Spliterator<long[]> trySplit() {
            int lo = myCurrentSite;
            int mid = lo + myMinSitesToProcess;
            if (mid < myFence) {
                myCurrentSite = mid;
                return new IBSSiteSpliterator(myGenotypes, myTaxonIndex, lo, mid, myProgressListener);
            } else {
                return null;
            }
        }

        @Override
        public long estimateSize() {
            return (long) (myFence - myCurrentSite);
        }

        @Override
        public int characteristics() {
            return IMMUTABLE;
        }
    }

}
