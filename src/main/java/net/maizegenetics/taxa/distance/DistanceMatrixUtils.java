// DistanceMatrixUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.taxa.distance;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Utility functions for distance matrices
 *
 * @author Alexei Drummond
 * @author Terry Casstevens
 */
public class DistanceMatrixUtils {

    private DistanceMatrixUtils() {
        // utility class
    }

    /**
     * Get Genetic Relationship Matrix (grm) file names.
     *
     * @param base filename base
     *
     * @return array of file names. Index 0 is the id filename (.grm.id); Index 1 is the binary matrix (.grm.bin); Index
     * 2 is the binary counts (.grm.N.bin); Index 3 is the raw (text) matrix (.grm.raw); Index 4 is the text matrix
     * (.txt)
     */
    public static String[] getGRMFilenames(String base) {

        String[] result = new String[5];
        String filename = base.toLowerCase();

        if (filename.endsWith(".txt")) {
            int txtIndex = filename.lastIndexOf(".txt");
            String temp = base.substring(0, txtIndex);
            result[0] = temp + ".grm.id";
            result[1] = temp + ".grm.bin";
            result[2] = temp + ".grm.N.bin";
            result[3] = temp + ".grm.raw";
            result[4] = filename;
            return result;
        }

        int grmIndex = filename.lastIndexOf(".grm");
        if (grmIndex == -1) {
            result[0] = base + ".grm.id";
            result[1] = base + ".grm.bin";
            result[2] = base + ".grm.N.bin";
            result[3] = base + ".grm.raw";
            result[4] = base + ".txt";
        } else {
            String temp = base.substring(0, grmIndex);
            result[0] = temp + ".grm.id";
            result[1] = temp + ".grm.bin";
            result[2] = temp + ".grm.N.bin";
            result[3] = temp + ".grm.raw";
            result[4] = temp + ".txt";
        }

        return result;

    }

    /**
     * Get DARwin file names.
     *
     * @param base filename base
     *
     * @return array of file names. Index 0 is the id filename (.don); Index 1 is the dissimilarity matrix (.dis)
     */
    public static String[] getDARwinFilenames(String base) {

        String[] result = new String[2];
        int index = base.lastIndexOf(".");
        if (index == -1) {
            result[0] = base + ".don";
            result[1] = base + ".dis";
            return result;
        } else {
            String temp = base.substring(0, index);
            result[0] = temp + ".don";
            result[1] = temp + ".dis";
            return result;
        }

    }

    /**
     * compute squared distance to second distance matrix. If both matrices have
     * the same size it is assumed that the order of the taxa is identical.
     */
    public static double squaredDistance(DistanceMatrix mat1, DistanceMatrix mat2, boolean weighted) {

        boolean aliasNeeded = false;
        if (mat1.getSize() != mat2.getSize()) {
            aliasNeeded = true;
        }

        int[] alias = null;

        if (aliasNeeded) {
            if (mat1.getSize() > mat2.getSize()) {
                //swap so mat1 is the smaller of the two
                DistanceMatrix temp = mat2;
                mat2 = mat1;
                mat1 = temp;
            }
            alias = new int[mat1.getSize()];
            for (int i = 0; i < alias.length; i++) {
                alias[i] = mat2.whichIdNumber(mat1.getTaxon(i).getName());
            }
        } else {
            alias = new int[mat1.getSize()];
            for (int i = 0; i < alias.length; i++) {
                alias[i] = i;
            }
        }

        double sum = 0;
        int ai;
        final double[][] mat1Distance = mat1.getDistances();
        final double[][] mat2Distance = mat2.getDistances();
        for (int i = 0; i < mat1.getSize() - 1; i++) {
            ai = alias[i];

            for (int j = i + 1; j < mat1.getSize(); j++) {
                double diff = mat1Distance[i][j] - mat2Distance[ai][alias[j]];
                double weight;
                if (weighted) {
                    // Fitch-Margoliash weight
                    // (variances proportional to distances)
                    weight = 1.0 / (mat1Distance[i][j] * mat2Distance[ai][alias[j]]);
                } else {
                    // Cavalli-Sforza-Edwards weight
                    // (homogeneity of variances)
                    weight = 1.0;
                }
                sum += weight * diff * diff;
            }
        }

        return 2.0 * sum; // we counted only half the matrix
    }

    /**
     * Returns a distance matrix with the specified taxa removed.
     */
    public static DistanceMatrix minus(DistanceMatrix parent, int taxaToRemove) {

        int size = parent.numberOfTaxa() - 1;

        double[][] distances = new double[size][size];
        Taxon[] ids = new Taxon[size];
        int counti = 0, countj = 0;
        for (int i = 0; i < size; i++) {
            if (counti == taxaToRemove) {
                counti += 1;
            }
            ids[i] = parent.getTaxon(counti);

            countj = 0;
            final double[][] parentDistance = parent.getDistances();
            for (int j = 0; j < size; j++) {
                if (countj == taxaToRemove) {
                    countj += 1;
                }
                distances[i][j] = parentDistance[counti][countj];
                countj += 1;
            }
            counti += 1;
        }
        TaxaList tl = new TaxaListBuilder().addAll(ids).build();
        DistanceMatrix smaller = new DistanceMatrix(distances, tl);

        return smaller;
    }

    /**
     * @param parent the DistanceMatrix from which to extract a subset
     * @param taxaToKeep an index of the taxa to keep
     *
     * @return A DistanceMatrix with all the taxa that are in both parent and taxaToKeep in the same order as taxaToKeep
     */
    public static DistanceMatrix keepTaxa(DistanceMatrix parent, int[] taxaToKeep) {
        int ntaxa = taxaToKeep.length;
        double[][] newDistances = new double[ntaxa][ntaxa];
        for (int r = 0; r < ntaxa; r++) {
            for (int c = 0; c < ntaxa; c++) {
                newDistances[r][c] = parent.getDistance(taxaToKeep[r], taxaToKeep[c]);
            }
        }

        TaxaListBuilder taxaBuilder = new TaxaListBuilder();
        for (int ndx : taxaToKeep) {
            taxaBuilder.add(parent.getTaxon(ndx));
        }
        TaxaList taxaListToKeep = taxaBuilder.build();

        return new DistanceMatrix(newDistances, taxaListToKeep);
    }

    /**
     * @param parent the DistanceMatrix from which to extract a subset
     * @param taxaToKeep a TaxaList of taxa to be kept
     *
     * @return a DistanceMatrix that contains only the taxa that are in both taxaToKeep and parent. The taxa will be in
     * the same order as taxaToKeep.
     */
    public static DistanceMatrix keepTaxa(DistanceMatrix parent, TaxaList taxaToKeep) {
        int[] keepIndex = taxaToKeep.stream()
                .mapToInt(t -> parent.whichIdNumber(t))
                .filter(i -> i > -1)
                .toArray();
        return keepTaxa(parent, keepIndex);
    }

    public static DistanceMatrix clusterBySmallestDistance(DistanceMatrix orig) {

        TaxaList taxa = orig.getTaxaList();
        int numTaxa = taxa.numberOfTaxa();

        TaxaPairLowestDistance[] lowValues = new TaxaPairLowestDistance[numTaxa];
        for (int t = 0; t < numTaxa; t++) {
            lowValues[t] = new TaxaPairLowestDistance(t);
        }

        for (int x = 0; x < numTaxa; x++) {
            for (int y = x + 1; y < numTaxa; y++) {

                float value = orig.getDistance(x, y);

                if (!Float.isNaN(value)) {

                    if (lowValues[x].myLowValue > value) {
                        lowValues[x].myLowValue = value;
                        lowValues[x].myTaxon2 = y;
                    }

                    if (lowValues[y].myLowValue > value) {
                        lowValues[y].myLowValue = value;
                        lowValues[y].myTaxon2 = x;
                    }

                }

            }
        }

        Arrays.sort(lowValues);

        List<List<Integer>> clusters = new ArrayList<>();
        List<Integer>[] whichCluster = new ArrayList[numTaxa];
        List<Integer> unknownList = new ArrayList<>();
        for (int t = 0; t < numTaxa; t++) {

            int taxon1 = lowValues[t].myTaxon1;
            int taxon2 = lowValues[t].myTaxon2;

            if (taxon2 == -1) {
                unknownList.add(taxon1);
            } else if (whichCluster[taxon1] == null && whichCluster[taxon2] == null) {
                List<Integer> temp = new ArrayList<>();
                temp.add(taxon1);
                temp.add(taxon2);
                clusters.add(temp);
                whichCluster[taxon1] = temp;
                whichCluster[taxon2] = temp;
            } else if (whichCluster[taxon1] == null) {
                whichCluster[taxon1] = whichCluster[taxon2];
                whichCluster[taxon1].add(taxon1);
            } else if (whichCluster[taxon2] == null) {
                whichCluster[taxon2] = whichCluster[taxon1];
                whichCluster[taxon2].add(taxon2);
            }

        }

        clusters.add(unknownList);

        DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(numTaxa);
        int count = 0;
        for (List<Integer> current : clusters) {
            int currentNumTaxa = current.size();
            for (int taxon = 0; taxon < currentNumTaxa; taxon++) {
                builder.addTaxon(taxa.get(current.get(taxon)));
                for (int x = taxon; x < currentNumTaxa; x++) {
                    builder.set(x + count, taxon + count, orig.getDistance(current.get(taxon), current.get(x)));
                }
            }
            count += currentNumTaxa;
        }

        return builder.build();

    }

    private static class TaxaPairLowestDistance implements Comparable<TaxaPairLowestDistance> {

        private final int myTaxon1;
        private int myTaxon2;
        private float myLowValue = Float.POSITIVE_INFINITY;

        public TaxaPairLowestDistance(int taxon1) {
            myTaxon1 = taxon1;
            myTaxon2 = -1;
            myLowValue = Float.POSITIVE_INFINITY;
        }

        @Override
        public int compareTo(TaxaPairLowestDistance o) {
            if (myLowValue < o.myLowValue) {
                return -1;
            } else if (myLowValue > o.myLowValue) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    /**
     * Calculates the IBS distance between two taxa with bitsets for for major
     * and minor allele
     *
     * @param iMajor
     * @param iMinor
     * @param jMajor
     * @param jMinor
     *
     * @return
     */
    public static double getIBSDistance(long[] iMajor, long[] iMinor, long[] jMajor, long[] jMinor) {
        int sameCnt = 0, diffCnt = 0, hetCnt = 0;
        for (int x = 0; x < iMajor.length; x++) {
            long same = (iMajor[x] & jMajor[x]) | (iMinor[x] & jMinor[x]);
            long diff = (iMajor[x] & jMinor[x]) | (iMinor[x] & jMajor[x]);
            long hets = same & diff;
            sameCnt += BitUtil.pop(same);
            diffCnt += BitUtil.pop(diff);
            hetCnt += BitUtil.pop(hets);
        }
        double identity = (double) (sameCnt + (hetCnt / 2)) / (double) (sameCnt + diffCnt + hetCnt);
        double dist = 1 - identity;
        return dist;
    }

    public static double getIBSDistance(BitSet iMajor, BitSet iMinor, BitSet jMajor, BitSet jMinor) {
        return getIBSDistance(iMajor.getBits(), iMinor.getBits(), jMajor.getBits(), jMinor.getBits());
    }
}
