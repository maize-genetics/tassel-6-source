/*
 *  AbstractGenotypeCallTable
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Spliterator;
import static java.util.Spliterator.IMMUTABLE;
import static java.util.Spliterator.ORDERED;
import static java.util.Spliterator.SIZED;
import static java.util.Spliterator.SUBSIZED;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Abstract implementation of methods of GenotypeCallTable.
 *
 * @author Terry Casstevens
 */
abstract class AbstractGenotypeCallTable implements GenotypeCallTable {

    private static final Logger myLogger = Logger.getLogger(AbstractGenotypeCallTable.class);
    private static final int DEFAULT_MAX_NUM_ALLELES = NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES;
    protected final int myTaxaCount;
    protected final int mySiteCount;
    private final String[][] myAlleleEncodings;
    private final boolean myIsPhased;
    private final AlleleFreqCache myAlleleFreqCache;

    AbstractGenotypeCallTable(int numTaxa, int numSites, boolean phased, String[][] alleleEncodings, int maxNumAlleles) {
        myTaxaCount = numTaxa;
        mySiteCount = numSites;
        myIsPhased = phased;
        myAlleleEncodings = alleleEncodings;
        myAlleleFreqCache = new AlleleFreqCache(this, maxNumAlleles);
    }

    AbstractGenotypeCallTable(int numTaxa, int numSites, boolean phased, String[][] alleleEncodings) {
        this(numTaxa, numSites, phased, alleleEncodings, DEFAULT_MAX_NUM_ALLELES);
    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        return GenotypeTableUtils.getDiploidValues(genotype(taxon, site));
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {
        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i - startSite] = genotype(taxon, i);
        }
        return result;
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = genotype(taxon, i);
        }
        return result;
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        String[][] alleleStates = alleleDefinitions();
        byte[] temp = genotypeArray(taxon, site);
        return alleleStates[0][temp[0]] + ":" + alleleStates[0][temp[1]];
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            if (i != startSite) {
                builder.append(";");
            }
            builder.append(genotypeAsString(taxon, i));
        }
        return builder.toString();
    }

    @Override
    public String genotypeAsStringRow(int taxon) {
        return genotypeAsStringRange(taxon, 0, mySiteCount);
    }

    @Override
    public String[] genotypeAsStringArray(int taxon, int site) {
        String[][] alleleStates = alleleDefinitions();
        byte[] temp = genotypeArray(taxon, site);
        return new String[]{alleleStates[0][temp[0]], alleleStates[0][temp[1]]};
    }

    @Override
    public String[] genotypeAsStringArray(int site, byte value) {
        String[][] alleleStates = alleleDefinitions();
        byte[] temp = GenotypeTableUtils.getDiploidValues(value);
        return new String[]{alleleStates[0][temp[0]], alleleStates[0][temp[1]]};
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        return myAlleleFreqCache.getAllelesSortedByFrequency(site);
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        byte[] values = genotypeArray(taxon, site);
        return values[0] != values[1];
    }

    @Override
    public int heterozygousCount(int site) {
        int result = 0;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            if (isHeterozygous(i, site)) {
                result++;
            }
        }
        return result;
    }

    @Override
    public boolean isPolymorphic(int site) {

        byte first = GenotypeTable.UNKNOWN_ALLELE;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            byte[] current = genotypeArray(i, site);
            if (current[0] != GenotypeTable.UNKNOWN_ALLELE) {
                if (first == GenotypeTable.UNKNOWN_ALLELE) {
                    first = current[0];
                } else if (first != current[0]) {
                    return true;
                }
            }
            if (current[1] != GenotypeTable.UNKNOWN_ALLELE) {
                if (first == GenotypeTable.UNKNOWN_ALLELE) {
                    first = current[1];
                } else if (first != current[1]) {
                    return true;
                }
            }
        }

        return false;

    }

    @Override
    public boolean isAllPolymorphic() {

        for (int i = 0, n = mySiteCount; i < n; i++) {
            if (!isPolymorphic(i)) {
                return false;
            }
        }

        return true;

    }

    @Override
    public boolean isPhased() {
        return myIsPhased;
    }

    @Override
    public boolean retainsRareAlleles() {
        return true;
    }

    @Override
    public String[][] alleleDefinitions() {
        return myAlleleEncodings;
    }

    @Override
    public String[] alleleDefinitions(int site) {
        if (myAlleleEncodings.length == 1) {
            return myAlleleEncodings[0];
        } else {
            return myAlleleEncodings[site];
        }
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        return alleleDefinitions(site)[value];
    }

    @Override
    public String diploidAsString(int site, byte value) {
        String[] alleleStates = alleleDefinitions(site);
        return alleleStates[(value >>> 4) & 0xf] + ":" + alleleStates[value & 0xf];
    }

    @Override
    public int maxNumAlleles() {
        return DEFAULT_MAX_NUM_ALLELES;
    }

    @Override
    public int totalGametesNonMissingForSite(int site) {

        int result = 0;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            byte[] current = genotypeArray(i, site);
            if (current[0] != GenotypeTable.UNKNOWN_ALLELE) {
                result++;
            }
            if (current[1] != GenotypeTable.UNKNOWN_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public int totalNonMissingForSite(int site) {

        int result = 0;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            byte current = genotype(i, site);
            if (current != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public byte[] majorAlleleForAllSites() {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = majorAllele(i);
        }
        return result;
    }

    @Override
    public byte[] minorAlleleForAllSites() {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = minorAllele(i);
        }
        return result;
    }

    @Override
    public byte[] thirdAlleleForAllSites() {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = thirdAllele(i);
        }
        return result;
    }

    @Override
    public byte thirdAllele(int site) {
        int[][] alleles = allelesSortedByFrequency(site);

        if (alleles[0].length > 2) {
            return (byte) alleles[0][2];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    @Override
    public int minorAlleleCount(int site) {

        int[][] alleles = allelesSortedByFrequency(site);

        if (alleles[0].length >= 2) {
            return alleles[1][1];
        } else {
            return 0;
        }

    }

    @Override
    public byte minorAllele(int site) {
        int[][] alleles = allelesSortedByFrequency(site);

        if (alleles[0].length >= 2) {
            return (byte) alleles[0][1];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    @Override
    public String minorAlleleAsString(int site) {
        return genotypeAsString(site, minorAllele(site));
    }

    @Override
    public byte[] minorAlleles(int site) {
        int[][] alleles = allelesSortedByFrequency(site);
        int resultSize = Math.max(0, alleles[0].length - 1);
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i + 1];
        }
        return result;
    }

    @Override
    public int majorAlleleCount(int site) {

        int[][] alleles = allelesSortedByFrequency(site);

        if (alleles[0].length >= 1) {
            return alleles[1][0];
        } else {
            return 0;
        }

    }

    @Override
    public byte majorAllele(int site) {
        int[][] alleles = allelesSortedByFrequency(site);

        if (alleles[0].length >= 1) {
            return (byte) alleles[0][0];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    @Override
    public String majorAlleleAsString(int site) {
        return genotypeAsString(site, majorAllele(site));
    }

    @Override
    public double majorAlleleFrequency(int site) {

        int[][] alleles = allelesSortedByFrequency(site);

        int numAlleles = alleles[0].length;
        if (numAlleles >= 1) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing = totalNonMissing + alleles[1][i];
            }
            return (double) alleles[1][0] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    @Override
    public double minorAlleleFrequency(int site) {

        int[][] alleles = allelesSortedByFrequency(site);

        int numAlleles = alleles[0].length;
        if (numAlleles >= 2) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing = totalNonMissing + alleles[1][i];
            }
            return (double) alleles[1][1] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    @Override
    public Object[][] genosSortedByFrequency(int site) {

        Integer ONE_INTEGER = 1;

        Map<String, Integer> diploidValueCounts = new HashMap<String, Integer>();
        for (int r = 0; r < myTaxaCount; r++) {
            String current = genotypeAsString(r, site);
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

    @Override
    public Object[][] genoCounts() {

        Map<String, Long> diploidValueCounts = new HashMap<String, Long>();
        for (int c = 0; c < mySiteCount; c++) {
            Object[][] diploids = genosSortedByFrequency(c);
            for (int i = 0; i < diploids[0].length; i++) {
                String current = (String) diploids[0][i];
                Long count = (long) ((Integer) diploids[1][i]).intValue();
                Long num = diploidValueCounts.get(current);
                if (num == null) {
                    diploidValueCounts.put(current, count);
                } else {
                    diploidValueCounts.put(current, (num + count));
                }
            }
        }

        Object[][] result = new Object[2][diploidValueCounts.size()];

        int i = 0;
        Iterator itr = diploidValueCounts.keySet().iterator();
        while (itr.hasNext()) {
            String key = (String) itr.next();
            Long count = diploidValueCounts.get(key);
            result[0][i] = key;
            result[1][i++] = count;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0, n = diploidValueCounts.size() - 1; k < n; k++) {

                if ((Long) result[1][k] < (Long) result[1][k + 1]) {

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

    @Override
    public Object[][] majorMinorCounts() {

        String[][] alleleStates = alleleDefinitions();

        if (alleleStates.length != 1) {
            return new Object[0][0];
        }

        long[][] counts = new long[16][16];

        for (int site = 0; site < mySiteCount; site++) {
            byte[] alleles = alleles(site);
            if ((alleles == null) || alleles.length == 0) {
                // do nothing
            } else if (alleles.length == 1) {
                counts[alleles[0]][alleles[0]]++;
            } else {
                counts[alleles[0]][alleles[1]]++;
            }
        }

        int numAlleles = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    numAlleles++;
                }
            }
        }

        Object[][] result = new Object[2][numAlleles];
        int nextResult = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    result[0][nextResult] = genotypeAsString(0, x) + ":" + genotypeAsString(0, y);
                    result[1][nextResult++] = counts[x][y];
                }
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if ((Long) result[1][k] < (Long) result[1][k + 1]) {

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

    @Override
    public int totalGametesNonMissingForTaxon(int taxon) {

        int result = 0;
        for (int i = 0, n = mySiteCount; i < n; i++) {
            byte[] current = genotypeArray(taxon, i);
            if (current[0] != GenotypeTable.UNKNOWN_ALLELE) {
                result++;
            }
            if (current[1] != GenotypeTable.UNKNOWN_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public int heterozygousCountForTaxon(int taxon) {
        int result = 0;
        for (int i = 0, n = mySiteCount; i < n; i++) {
            if (isHeterozygous(taxon, i)) {
                result++;
            }
        }
        return result;
    }

    @Override
    public int totalNonMissingForTaxon(int taxon) {

        int result = 0;
        for (int i = 0, n = mySiteCount; i < n; i++) {
            byte current = genotype(taxon, i);
            if (current != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public byte[] alleles(int site) {
        int[][] alleles = allelesSortedByFrequency(site);
        int resultSize = alleles[0].length;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i];
        }
        return result;
    }

    @Override
    public int numberOfSites() {
        return mySiteCount;
    }

    @Override
    public int numberOfTaxa() {
        return myTaxaCount;
    }

    @Override
    public byte[] genotypeForAllSites(int taxon) {
        int numSites = numberOfSites();
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = genotype(taxon, i);
        }
        return result;
    }

    @Override
    public byte[] genotypeForSiteRange(int taxon, int start, int end) {
        int numSites = end - start;
        byte[] result = new byte[numSites];
        for (int i = start; i < end; i++) {
            result[i] = genotype(taxon, i);
        }
        return result;
    }

    @Override
    public byte[] genotypeForAllTaxa(int site) {
        int numTaxa = numberOfTaxa();
        byte[] result = new byte[numTaxa];
        for (int i = 0; i < numTaxa; i++) {
            result[i] = genotype(i, site);
        }
        return result;
    }

    @Override
    public Stats siteStats(int site) {
        return AlleleFreqCache.allelesSortedByFrequencyAndCountsNucleotide(site, genotypeForAllTaxa(site));
    }

    @Override
    public Stats taxonStats(int taxon) {
        return AlleleFreqCache.allelesSortedByFrequencyAndCountsNucleotide(taxon, genotypeForAllSites(taxon));
    }

    @Override
    public Stream<Byte> stream() {
        return StreamSupport.stream(spliterator(), true);
    }

    @Override
    public Stream<Byte> stream(int taxon) {
        return StreamSupport.stream(new AbstractGenotypeCallTableSpliterator<>(taxon, 0, numberOfSites(), taxon, numberOfSites()), true);
    }

    public Spliterator<Byte> spliterator() {
        return new AbstractGenotypeCallTableSpliterator<>(0, 0, numberOfSites(), numberOfTaxa() - 1, numberOfSites());
    }

    class AbstractGenotypeCallTableSpliterator<T extends Byte> implements Spliterator<Byte> {

        protected int myTaxaOrigin;
        protected int mySiteOrigin;
        protected final int myNumSites;
        protected final int myTaxaFence;
        protected final int mySiteFence;

        AbstractGenotypeCallTableSpliterator(int taxaOrigin, int siteOrigin, int numSites, int taxaFence, int siteFence) {
            myTaxaOrigin = taxaOrigin;
            mySiteOrigin = siteOrigin;
            myNumSites = numSites;
            myTaxaFence = taxaFence;
            mySiteFence = siteFence;
        }

        @Override
        public void forEachRemaining(Consumer<? super Byte> action) {
            for (; myTaxaOrigin < myTaxaFence; myTaxaOrigin++) {
                for (; mySiteOrigin < myNumSites; mySiteOrigin++) {
                    action.accept(genotype(myTaxaOrigin, mySiteOrigin));
                }
                mySiteOrigin = 0;
            }
            for (; mySiteOrigin < mySiteFence; mySiteOrigin++) {
                action.accept(genotype(myTaxaOrigin, mySiteOrigin));
            }
        }

        @Override
        public boolean tryAdvance(Consumer<? super Byte> action) {
            if (((myTaxaOrigin < myTaxaFence) && (mySiteOrigin < myNumSites))
                    || ((myTaxaOrigin == myTaxaFence) && (mySiteOrigin < mySiteFence))) {
                action.accept(genotype(myTaxaOrigin, mySiteOrigin));
                mySiteOrigin++;
                if (mySiteOrigin >= myNumSites) {
                    mySiteOrigin = 0;
                    myTaxaOrigin++;
                }
                return true;
            } else {
                return false;
            }
        }

        @Override
        public Spliterator<Byte> trySplit() {
            long size = estimateSize();
            if (size > 1) {
                size >>>= 1;
                int loTaxa = myTaxaOrigin;
                int loSite = mySiteOrigin;
                int midTaxa = myTaxaOrigin;
                int midSite = mySiteOrigin;
                midTaxa += size / myNumSites;
                midSite += size % myNumSites;
                if (midSite > myNumSites) {
                    midTaxa++;
                    midSite -= myNumSites;
                }
                myTaxaOrigin = midTaxa;
                mySiteOrigin = midSite;
                return new AbstractGenotypeCallTableSpliterator<>(loTaxa, loSite, myNumSites, midTaxa, midSite);
            } else {
                return null;
            }
        }

        @Override
        public long estimateSize() {
            return (long) (myTaxaFence - myTaxaOrigin) * (long) myNumSites - (long) mySiteOrigin + (long) mySiteFence;
        }

        @Override
        public int characteristics() {
            return ORDERED | SIZED | IMMUTABLE | SUBSIZED;
        }
    }

}
