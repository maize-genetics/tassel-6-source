/*
 *  EndelmanDistanceMatrix
 *
 *  Created on April 19, 2021
 */
package net.maizegenetics.analysis.distance;

import net.maizegenetics.dna.factor.FactorTable;
import net.maizegenetics.dna.factor.site.AlleleFreq;
import net.maizegenetics.dna.factor.site.FactorSite;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;
import net.maizegenetics.util.GeneralAnnotationStorage;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

import java.util.Arrays;
import java.util.Optional;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * @author Terry Casstevens
 */
public class EndelmanDistanceMatrix {

    private static final Logger myLogger = Logger.getLogger(EndelmanDistanceMatrix.class);

    private static final int DEFAULT_MAX_ALLELES = 255;

    /**
     * Compute Endelman kinship for all pairs of taxa. Missing sites are
     * ignored. http://www.g3journal.org/content/2/11/1405.full.pdf Equation-13
     *
     * @return Endelman Kinship Matrix
     */
    private EndelmanDistanceMatrix() {
        // utility
    }

    /**
     * Compute Endelman Kinship Matrix. Maximum alleles per site to evaluate
     * defaults to DEFAULT_MAX_ALLELES.
     *
     * @param genotype Factor Table used to compute kinship
     *
     * @return Endelman Kinship Matrix
     */
    public static DistanceMatrix getInstance(FactorTable genotype) {
        return getInstance(genotype, DEFAULT_MAX_ALLELES, null);
    }

    /**
     * Compute Endelman Kinship Matrix
     *
     * @param genotype Factor Table used to compute kinship
     * @param maxAlleles maximum alleles per site to evaluate. i.e. Set to 3 to
     * evaluate the three most frequent allele states.
     *
     * @return Endelman Kinship Matrix
     */
    public static DistanceMatrix getInstance(FactorTable genotype, int maxAlleles) {
        return getInstance(genotype, maxAlleles, null);
    }

    /**
     * Compute Endelman Kinship Matrix. Maximum alleles per site to evaluate
     * defaults to 2.
     *
     * @param genotype Factor Table used to compute kinship
     * @param listener Progress listener
     *
     * @return Endelman Kinship Matrix
     */
    public static DistanceMatrix getInstance(FactorTable genotype, ProgressListener listener) {
        return computeEndelmanDistances(genotype, DEFAULT_MAX_ALLELES, listener);
    }

    /**
     * Compute Endelman Kinship Matrix
     *
     * @param genotype Factor Table used to compute kinship
     * @param maxAlleles maximum alleles per site to evaluate. i.e. Set to 3 to
     * evaluate the three most frequent allele states.
     * @param listener Progress listener
     *
     * @return Endelman Kinship Matrix
     */
    public static DistanceMatrix getInstance(FactorTable genotype, int maxAlleles, ProgressListener listener) {
        return computeEndelmanDistances(genotype, maxAlleles, listener);
    }

    private static DistanceMatrix computeEndelmanDistances(FactorTable genotype, int maxAlleles, ProgressListener listener) {

        if ((maxAlleles < 2) || (maxAlleles > 255)) {
            throw new IllegalArgumentException("EndelmanDistanceMatrix: computeEndelmanDistances: max alleles must be between 2 and 6 inclusive.");
        }

        int numSeqs = genotype.numTaxa();

        long estimatedNumMinutesToRun = Math.round((double) numSeqs * ((double) numSeqs + 1.0) / 2.0 * (double) genotype.numFactors() / (double) NUM_CORES_TO_USE / 85000000000.0);
        if (estimatedNumMinutesToRun < 60l) {
            myLogger.info("EndelmanDistanceMatrix: estimated time: " + estimatedNumMinutesToRun + " minutes");
        } else {
            myLogger.info("EndelmanDistanceMatrix: estimated time: " + (estimatedNumMinutesToRun / 60) + " hours " + (estimatedNumMinutesToRun % 60) + " minutes");
        }

        long time = System.currentTimeMillis();

        //
        // Sets up parellel stream to divide up sites for processing.
        // Also reduces the distance sums and sum of frequencies into one instance.
        //
        Optional<CountersDistances> optional = stream(genotype, maxAlleles, listener).reduce((CountersDistances t, CountersDistances u) -> {
            t.addAll(u);
            return t;
        });

        if (!optional.isPresent()) {
            return null;
        }
        CountersDistances counters = optional.get();
        double sumpk = counters.mySumPi;
        float[] distances = counters.myDistances;

        //
        // This does the final division of the frequency sum into
        // the distance sums.
        //
        sumpk *= 2.0;

        GeneralAnnotationStorage.Builder annotations = GeneralAnnotationStorage.getBuilder();
        annotations.addAnnotation(DistanceMatrixBuilder.MATRIX_TYPE, KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString());
        annotations.addAnnotation(DistanceMatrixBuilder.CENTERED_IBS_SUMPK, sumpk);

        DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(genotype.taxa());
        builder.annotation(annotations.build());
        int index = 0;
        for (int t = 0; t < numSeqs; t++) {
            for (int i = 0, n = numSeqs - t; i < n; i++) {
                if (i == 0 && t == 32) {
                    System.out.println("taxa 32: distances: " + distances[index] + "  sumpk: " + sumpk);
                }
                builder.set(t, t + i, distances[index] / sumpk);
                index++;
            }
        }

        long actualNumMinutesToRun = Math.round((System.currentTimeMillis() - time) / 60000l);
        if (actualNumMinutesToRun < 60l) {
            myLogger.info("EndelmanDistanceMatrix: actual time: " + actualNumMinutesToRun + " minutes");
        } else {
            myLogger.info("EndelmanDistanceMatrix: actual time: " + (actualNumMinutesToRun / 60) + " hours " + (actualNumMinutesToRun % 60) + " minutes");
        }

        return builder.build();

    }

    public static DistanceMatrix subtractEndelmanDistance(DistanceMatrix[] matrices, DistanceMatrix superMatrix, ProgressListener listener) {

        int numTaxa = superMatrix.numberOfTaxa();
        String matrixType = superMatrix.annotations().getTextAnnotation(DistanceMatrixBuilder.MATRIX_TYPE)[0];
        if (!matrixType.equals(KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString())) {
            throw new IllegalArgumentException("subtractEndelmanDistance: superset matrix must be matrix type: " + KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString());
        }
        for (DistanceMatrix current : matrices) {
            int currentNumTaxa = current.numberOfTaxa();
            if (currentNumTaxa != numTaxa) {
                throw new IllegalArgumentException("subtractEndelmanDistance: subset and superset must have same number of taxa.");
            }
            String[] currentMatrixType = current.annotations().getTextAnnotation(DistanceMatrixBuilder.MATRIX_TYPE);
            if (currentMatrixType.length == 0) {
                throw new IllegalArgumentException("subtractEndelmanDistance: subset matrix must be created with a more recent build of Tassel that adds neccessary annotations to the matrix");
            }
            if (!matrixType.equals(currentMatrixType[0])) {
                throw new IllegalArgumentException("subtractEndelmanDistance: subset matrix must be matrix type: " + KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString());
            }
        }

        TaxaList superTaxaList = superMatrix.getTaxaList();
        for (DistanceMatrix current : matrices) {
            TaxaList subsetTaxaList = current.getTaxaList();
            for (int t = 0; t < numTaxa; t++) {
                if (!superTaxaList.get(t).equals(subsetTaxaList.get(t))) {
                    throw new IllegalArgumentException("subtractEndelmanDistance: superset taxon: " + superTaxaList.get(t).getName() + " doesn't match subset taxon: " + subsetTaxaList.taxaName(t));
                }
            }
        }

        DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(superTaxaList);

        //
        // This does the final division of the frequency sum into
        // the distance sums.
        //
        int numMatrices = matrices.length;
        double superSumpk = superMatrix.annotations().getQuantAnnotation(DistanceMatrixBuilder.CENTERED_IBS_SUMPK)[0];
        double resultSumpk = superSumpk;
        double[] matricesSumpk = new double[numMatrices];
        for (int i = 0; i < numMatrices; i++) {
            matricesSumpk[i] = matrices[i].annotations().getQuantAnnotation(DistanceMatrixBuilder.CENTERED_IBS_SUMPK)[0];
            resultSumpk -= matricesSumpk[i];
        }

        GeneralAnnotationStorage.Builder resultAnnotations = GeneralAnnotationStorage.getBuilder();
        resultAnnotations.addAnnotation(DistanceMatrixBuilder.MATRIX_TYPE, KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString());
        resultAnnotations.addAnnotation(DistanceMatrixBuilder.CENTERED_IBS_SUMPK, resultSumpk);
        builder.annotation(resultAnnotations.build());

        for (int t = 0; t < numTaxa; t++) {
            for (int i = 0, n = numTaxa - t; i < n; i++) {
                double resultValue = superMatrix.getDistance(t, t + i) * superSumpk;
                for (int j = 0; j < numMatrices; j++) {
                    resultValue -= (matrices[j].getDistance(t, t + i) * matricesSumpk[j]);
                }
                builder.set(t, t + i, resultValue / resultSumpk);
            }
        }

        return builder.build();

    }

    public static DistanceMatrix addEndelmanDistance(DistanceMatrix[] matrices, ProgressListener listener) {

        int numTaxa = matrices[0].numberOfTaxa();
        String matrixType = matrices[0].annotations().getTextAnnotation(DistanceMatrixBuilder.MATRIX_TYPE)[0];
        if (!matrixType.equals(KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString())) {
            throw new IllegalArgumentException("addEndelmanDistance: superset matrix must be matrix type: " + KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString());
        }
        for (int i = 1; i < matrices.length; i++) {
            DistanceMatrix current = matrices[i];
            int currentNumTaxa = current.numberOfTaxa();
            if (currentNumTaxa != numTaxa) {
                throw new IllegalArgumentException("addEndelmanDistance: all matrices must have same number of taxa.");
            }
            String[] currentMatrixType = current.annotations().getTextAnnotation(DistanceMatrixBuilder.MATRIX_TYPE);
            if (currentMatrixType.length == 0) {
                throw new IllegalArgumentException("addEndelmanDistance: matrix must be created with a more recent build of Tassel that adds neccessary annotations to the matrix");
            }
            if (!matrixType.equals(currentMatrixType[0])) {
                throw new IllegalArgumentException("addEndelmanDistance: matrix must be matrix type: " + KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString());
            }
        }

        TaxaList superTaxaList = matrices[0].getTaxaList();
        for (int i = 1; i < matrices.length; i++) {
            DistanceMatrix current = matrices[i];
            TaxaList subsetTaxaList = current.getTaxaList();
            for (int t = 0; t < numTaxa; t++) {
                if (!superTaxaList.get(t).equals(subsetTaxaList.get(t))) {
                    throw new IllegalArgumentException("addEndelmanDistance: superset taxon: " + superTaxaList.get(t).getName() + " doesn't match subset taxon: " + subsetTaxaList.taxaName(t));
                }
            }
        }

        DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(superTaxaList);

        //
        // This does the final division of the frequency sum into
        // the distance sums.
        //
        int numMatrices = matrices.length;
        double resultSumpk = 0.0;
        double[] matricesSumpk = new double[numMatrices];
        for (int i = 0; i < numMatrices; i++) {
            matricesSumpk[i] = matrices[i].annotations().getQuantAnnotation(DistanceMatrixBuilder.CENTERED_IBS_SUMPK)[0];
            resultSumpk += matricesSumpk[i];
        }

        GeneralAnnotationStorage.Builder resultAnnotations = GeneralAnnotationStorage.getBuilder();
        resultAnnotations.addAnnotation(DistanceMatrixBuilder.MATRIX_TYPE, KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString());
        resultAnnotations.addAnnotation(DistanceMatrixBuilder.CENTERED_IBS_SUMPK, resultSumpk);
        builder.annotation(resultAnnotations.build());

        for (int t = 0; t < numTaxa; t++) {
            for (int i = 0, n = numTaxa - t; i < n; i++) {
                double resultValue = 0.0;
                for (int j = 0; j < numMatrices; j++) {
                    resultValue += (matrices[j].getDistance(t, t + i) * matricesSumpk[j]);
                }
                builder.set(t, t + i, resultValue / resultSumpk);
            }
        }

        return builder.build();

    }

    protected static void fireProgress(int percent, ProgressListener listener) {
        if (listener != null) {
            if (percent > 100) {
                percent = 100;
            }
            listener.progress(percent, null);
        }
    }

    //
    // Each CPU thread (process) creates an instance of this class
    // to acculate terms of the Endelman equation. These are
    // combined with addAll() to result in one instance at the end.
    //
    private static class CountersDistances {

        private double mySumPi = 0.0;
        private final float[] myDistances;
        private final int myNumTaxa;

        public CountersDistances(int numTaxa) {
            myNumTaxa = numTaxa;
            myDistances = new float[myNumTaxa * (myNumTaxa + 1) / 2];
        }

        public void addAll(CountersDistances counters) {
            float[] otherDistances = counters.myDistances;
            for (int t = 0, n = myDistances.length; t < n; t++) {
                myDistances[t] += otherDistances[t];
            }
            mySumPi += counters.mySumPi;
        }

    }

    private static final byte calculateCount(byte allele, byte value1, byte value2) {
        if (allele == (byte) 0xFF || (value1 == (byte) 0xFF && value2 == (byte) 0xFF)) return 7;
        byte result = 0;
        if (value1 == allele) {
            result = 2;
        }
        if (value2 == allele) {
            result += 2;
        }
        if (result == 0) {
            result = 1;
        }
        return result;
    }

    private static final int NUM_CORES_TO_USE = TasselPrefs.getMaxThreads();

    //
    // Used to report progress.  This is not thread-safe but
    // works well enough for this purpose.
    //
    private static int myNumSitesProcessed = 0;

    //
    // Creates stream from EndelmanSiteSpliterator and Genotype Table
    //
    private static Stream<CountersDistances> stream(FactorTable genotypes, int maxAlleles, ProgressListener listener) {
        myNumSitesProcessed = 0;
        return StreamSupport.stream(new EndelmanSiteSpliterator(genotypes, 0, genotypes.numFactors(), maxAlleles, listener), true);
    }

    //
    // Spliterator that splits the sites into halves each time for
    // processing.
    //
    static class EndelmanSiteSpliterator implements Spliterator<CountersDistances> {

        private int myCurrentSite;
        private final int myFence;
        private final FactorTable myGenotypes;
        private final int myNumTaxa;
        private final int myNumSites;
        private final int myMaxAlleles;
        private final ProgressListener myProgressListener;
        private final int myMinSitesToProcess;
        private final int myNumSitesPerBlockForProgressReporting;

        EndelmanSiteSpliterator(FactorTable genotypes, int currentIndex, int fence, int maxAlleles, ProgressListener listener) {
            myGenotypes = genotypes;
            myNumTaxa = myGenotypes.numTaxa();
            myNumSites = myGenotypes.numFactors();
            myCurrentSite = currentIndex;
            myFence = fence;
            myMaxAlleles = maxAlleles;
            myProgressListener = listener;
            myMinSitesToProcess = Math.max(myNumSites / NUM_CORES_TO_USE, 1000);
            myNumSitesPerBlockForProgressReporting = Math.max((myFence - myCurrentSite) / 10, 100);
        }

        @Override
        public void forEachRemaining(Consumer<? super CountersDistances> action) {

            CountersDistances result = new CountersDistances(myNumTaxa);
            float[] distances = result.myDistances;
            double[] sumpi = new double[1];

            float[] answer1 = new float[32768];
            float[] answer2 = new float[32768];
            float[] answer3 = new float[32768];

            for (; myCurrentSite < myFence; ) {

                int currentBlockFence = Math.min(myCurrentSite + myNumSitesPerBlockForProgressReporting, myFence);

                int numSitesProcessed = currentBlockFence - myCurrentSite;

                for (; myCurrentSite < currentBlockFence; ) {

                    FactorSite factorSite = myGenotypes.site(myCurrentSite);
                    int[][] alleles = AlleleFreq.alleleFreq(factorSite, myMaxAlleles);
                    int totalNumAlleles = Math.min(alleles[0].length - 1, myMaxAlleles - 1);

                    int totalAlleleCount = 0;
                    for (int i = 0; i < alleles[1].length; i++) {
                        //System.out.println(i + ": " + alleles[0][i]);
                        totalAlleleCount += alleles[1][i];
                    }

                    for (int alleleBlock = 0; alleleBlock < totalNumAlleles; alleleBlock += 15) {

                        //
                        // Pre-calculates possible terms and gets counts for
                        // three blocks for five (pseudo-)sites.
                        //
                        //getBlockOfSites(FactorSite factorSite, int[][] alleles, int startAllele, int numAlleles, int totalAlleleCount, double[] sumpi)
                        int startAllele = alleleBlock;
                        int numAlleles = Math.min(5, totalNumAlleles - startAllele);
                        Tuple<short[], float[]> firstBlock = getBlockOfSites(factorSite, alleles, startAllele, numAlleles, totalAlleleCount, sumpi);
                        float[] possibleTerms = firstBlock.y;
                        short[] alleleCount1 = firstBlock.x;

                        //System.out.println("myCurrentSite: " + myCurrentSite);
                        if (myCurrentSite == 3092) {
                            System.out.println("myCurrentSite: " + myCurrentSite + " possibleTerms");
                            for (float term : possibleTerms) {
                                System.out.print(term + " ");
                            }
                            System.out.println();
                            System.out.println("myCurrentSite: " + myCurrentSite + " alleleCount1");
                            for (short count : alleleCount1) {
                                System.out.print(count + " ");
                            }
                            System.out.println();
                        }

                        startAllele = Math.min(startAllele + 5, totalNumAlleles);
                        numAlleles = Math.min(5, totalNumAlleles - startAllele);
                        Tuple<short[], float[]> secondBlock = getBlockOfSites(factorSite, alleles, startAllele, numAlleles, totalAlleleCount, sumpi);
                        float[] possibleTerms2 = secondBlock.y;
                        short[] alleleCount2 = secondBlock.x;

                        if (myCurrentSite == 3092) {
                            System.out.println("myCurrentSite: " + myCurrentSite + " possibleTerms2");
                            for (float term : possibleTerms2) {
                                System.out.print(term + " ");
                            }
                            System.out.println();
                            System.out.println("myCurrentSite: " + myCurrentSite + " alleleCount2");
                            for (short count : alleleCount2) {
                                System.out.print(count + " ");
                            }
                            System.out.println();
                        }

                        startAllele = Math.min(startAllele + 10, totalNumAlleles);
                        numAlleles = Math.min(5, totalNumAlleles - startAllele);
                        Tuple<short[], float[]> thirdBlock = getBlockOfSites(factorSite, alleles, startAllele, numAlleles, totalAlleleCount, sumpi);
                        float[] possibleTerms3 = thirdBlock.y;
                        short[] alleleCount3 = thirdBlock.x;

                        if (myCurrentSite == 3092) {
                            System.out.println("myCurrentSite: " + myCurrentSite + " possibleTerms3");
                            for (float term : possibleTerms3) {
                                System.out.print(term + " ");
                            }
                            System.out.println();
                            System.out.println("myCurrentSite: " + myCurrentSite + " alleleCount3");
                            for (short count : alleleCount3) {
                                System.out.print(count + " ");
                            }
                            System.out.println();
                        }

                        //
                        // Using possible terms, calculates all possible answers
                        // for each site block.
                        //
                        for (int i = 0; i < 32768; i++) {
                            answer1[i] = possibleTerms[(i & 0x7000) >>> 12] + possibleTerms[((i & 0xE00) >>> 9) | 0x8] + possibleTerms[((i & 0x1C0) >>> 6) | 0x10] + possibleTerms[((i & 0x38) >>> 3) | 0x18] + possibleTerms[(i & 0x7) | 0x20];
                            answer2[i] = possibleTerms2[(i & 0x7000) >>> 12] + possibleTerms2[((i & 0xE00) >>> 9) | 0x8] + possibleTerms2[((i & 0x1C0) >>> 6) | 0x10] + possibleTerms2[((i & 0x38) >>> 3) | 0x18] + possibleTerms2[(i & 0x7) | 0x20];
                            answer3[i] = possibleTerms3[(i & 0x7000) >>> 12] + possibleTerms3[((i & 0xE00) >>> 9) | 0x8] + possibleTerms3[((i & 0x1C0) >>> 6) | 0x10] + possibleTerms3[((i & 0x38) >>> 3) | 0x18] + possibleTerms3[(i & 0x7) | 0x20];
                        }

                        //
                        // Iterates through all pair-wise combinations of taxa adding
                        // distance comparisons and site counts.
                        //
                        int index = 0;
                        for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                            //
                            // Can skip inter-loop if all fifteen sites for first
                            // taxon is Unknown diploid allele values
                            //
                            if ((alleleCount1[firstTaxa] != 0x7FFF) || (alleleCount2[firstTaxa] != 0x7FFF) || (alleleCount3[firstTaxa] != 0x7FFF)) {
                                for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                                    //
                                    // Combine first taxon's allele counts with
                                    // second taxon's major allele counts to
                                    // create index into pre-calculated answers
                                    //
                                    distances[index] += answer1[alleleCount1[firstTaxa] | alleleCount1[secondTaxa]] + answer2[alleleCount2[firstTaxa] | alleleCount2[secondTaxa]] + answer3[alleleCount3[firstTaxa] | alleleCount3[secondTaxa]];
                                    index++;
                                }
                            } else {
                                index += myNumTaxa - firstTaxa;
                            }
                        }

                    }

                    myCurrentSite++;

                }

                myNumSitesProcessed += numSitesProcessed;
                fireProgress((int) ((double) myNumSitesProcessed / (double) myNumSites * 100.0), myProgressListener);

            }

            result.mySumPi = sumpi[0];
            action.accept(result);

        }

        private static final int NUM_SITES_PER_BLOCK = 5;

        private Tuple<short[], float[]> getBlockOfSites(FactorSite factorSite, int[][] alleles, int startAllele, int numAlleles, int totalAlleleCount, double[] sumpi) {

            int currentSiteNum = 0;

            //
            // This hold possible terms for the Endelman summation given
            // site's allele frequency.  First three bits
            // identifies relative site (0, 1, 2, 3, 4).  Remaining three bits
            // the allele counts encoding.
            //
            float[] possibleTerms = new float[40];

            //
            // This holds count of allele for each taxa.
            // Each short holds count (0, 1, 2, 3) for all four sites
            // at given taxon.  The count encodings are stored in three
            // bits each.
            //
            short[] alleleCount = new short[myNumTaxa];

            //
            // This initializes the counts to 0x7FFF.  That means
            // diploid allele values for the five pseudo-sites are Unknown.
            //
            Arrays.fill(alleleCount, (short) 0x7FFF);

            for (int a = 0; a < numAlleles; a++) {

                byte allele = (byte) alleles[0][a + startAllele];
                float alleleFreq = (float) alleles[1][a + startAllele] / (float) totalAlleleCount;
                float alleleFreqTimes2 = alleleFreq * 2.0f;
                sumpi[0] += alleleFreq * (1.0 - alleleFreq);

                //
                // Temporarily stores component terms of equation for
                // individual allele counts (0, 1, 2)
                //
                float[] term = new float[3];

                //System.out.println((a + startAllele) + ": allele: " + allele);

                //
                // If allele is Unknown, the entire
                // site is skipped.
                //
                if (allele != (byte) 0xFF) {

                    term[0] = 0.0f - alleleFreqTimes2;
                    term[1] = 1.0f - alleleFreqTimes2;
                    term[2] = 2.0f - alleleFreqTimes2;

                    //
                    // Pre-calculates all possible terms of the summation
                    // for this current (pseudo-) site.
                    // Counts (0,0; 0,1; 0,2; 1,1; 1,2; 2,2)
                    //
                    int siteNumIncrement = currentSiteNum * 8;
                    possibleTerms[siteNumIncrement + 1] = term[0] * term[0];
                    possibleTerms[siteNumIncrement + 3] = term[0] * term[1];
                    possibleTerms[siteNumIncrement + 5] = term[0] * term[2];
                    possibleTerms[siteNumIncrement + 2] = term[1] * term[1];
                    possibleTerms[siteNumIncrement + 6] = term[1] * term[2];
                    possibleTerms[siteNumIncrement + 4] = term[2] * term[2];

                    //
                    // Records allele counts for current site in
                    // three bits.
                    //
                    //int temp = (allele & 0x7) << 6;
                    int shift = (NUM_SITES_PER_BLOCK - currentSiteNum - 1) * 3;
                    int mask = ~(0x7 << shift) & 0x7FFF;
                    for (int i = 0; i < myNumTaxa; i++) {
                        byte[] taxonAlleles = factorSite.genotype(i);
                        alleleCount[i] = (short) (alleleCount[i] & (mask | calculateCount(allele, taxonAlleles[0], taxonAlleles[1]) << shift));
                        //alleleCount[i] = (short) (alleleCount[i] & (mask | PRECALCULATED_COUNTS[temp | ((genotypes[i] & 0x70) >>> 1) | (genotypes[i] & 0x7)] << shift));
                    }
                }

                currentSiteNum++;
            }

            return new Tuple<>(alleleCount, possibleTerms);

        }

        @Override
        public boolean tryAdvance(Consumer<? super CountersDistances> action) {
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
        public Spliterator<CountersDistances> trySplit() {
            int lo = myCurrentSite;
            int mid = lo + myMinSitesToProcess;
            if (mid < myFence) {
                myCurrentSite = mid;
                return new EndelmanSiteSpliterator(myGenotypes, lo, mid, myMaxAlleles, myProgressListener);
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
