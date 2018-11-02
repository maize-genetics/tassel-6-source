/*
 *  HDF5ByteGenotype
 */
package net.maizegenetics.dna.snp.genotypecall;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.Tassel5HDF5Constants;
import net.maizegenetics.util.HDF5Utils;

import java.util.concurrent.ExecutionException;
import java.util.function.Consumer;
import java.util.Spliterator;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.apache.log4j.Logger;

/**
 * HDF5 implementation of GenotypeTable. Uses caching of GenotypeTable,
 * alleleCounts, MAF, and siteCoverage
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
class HDF5ByteGenotypeCallTable extends AbstractGenotypeCallTable {

    private static final Logger myLogger = Logger.getLogger(HDF5ByteGenotypeCallTable.class);

    private static final int SHIFT_AMOUNT = 16;

    private final String[] genotypePaths;
    /**
     * Byte representations of DNA sequences are stored in blocks of 65536 sites
     */
    private static final int HDF5_GENOTYPE_BLOCK_SIZE = 1 << SHIFT_AMOUNT;
    public static final int SITE_BLOCK_MASK = ~(HDF5_GENOTYPE_BLOCK_SIZE - 1);
    private final IHDF5Reader myHDF5Reader;

    private final LoadingCache<Long, byte[]> myGenoCache;
    private final CacheLoader<Long, byte[]> myGenoLoader = new CacheLoader<Long, byte[]>() {
        @Override
        public byte[] load(Long key) {
            long offset = getSiteStartFromKey(key) << SHIFT_AMOUNT;
            byte[] data;
            synchronized (myHDF5Reader) {
                data = myHDF5Reader.int8().readArrayBlockWithOffset(getTaxaGenoPath(getTaxonFromKey(key)), HDF5_GENOTYPE_BLOCK_SIZE, offset);
            }
            return data;
        }
    };

    private final LoadingCache<Integer, SiteBlockAttr> mySiteAnnoCache; //key = site
    private final CacheLoader<Integer, SiteBlockAttr> siteAnnotLoader = new CacheLoader<Integer, SiteBlockAttr>() {
        int lastCachedStartSite = Integer.MIN_VALUE;
        int[][] af;
        byte[][] afOrder;
        float[] maf;
        float[] paf;

        @Override
        public SiteBlockAttr load(Integer key) {
            int startSite = getStartSite(key);
            int length = Math.min(HDF5_GENOTYPE_BLOCK_SIZE, numberOfSites() - startSite);
            System.out.println("Reading from HDF5 site anno:" + startSite);
            System.out.println("");
            synchronized (myHDF5Reader) {
                af = myHDF5Reader.int32().readMatrixBlockWithOffset(Tassel5HDF5Constants.ALLELE_CNT, 6, length, 0l, startSite);
                afOrder = myHDF5Reader.int8().readMatrixBlockWithOffset(Tassel5HDF5Constants.ALLELE_FREQ_ORD, 6, length, 0l, startSite);
                maf = myHDF5Reader.float32().readArrayBlockWithOffset(Tassel5HDF5Constants.MAF, length, startSite);
                paf = myHDF5Reader.float32().readArrayBlockWithOffset(Tassel5HDF5Constants.SITECOV, length, startSite);
                lastCachedStartSite = startSite;
            }
            return new SiteBlockAttr(startSite, afOrder, af, maf, paf);
        }
    };

    private class SiteBlockAttr {

        private final int startSite;  //4
        private final byte[][] myAlleleFreqOrder;  //[
        private final int[][] myAlleleCnt;  //2-6*4=24,  sorted by allele frequency of myAlleleFreqOrder
        private final float[] maf; //4
        private final float[] siteCov;  //4

        public SiteBlockAttr(int startSite, byte[][] myAlleleFreqOrder, int[][] myAlleleCnt,
                float[] maf, float[] siteCov) {
            this.startSite = startSite;
            this.myAlleleFreqOrder = myAlleleFreqOrder;
            this.myAlleleCnt = myAlleleCnt;
            this.maf = maf;
            this.siteCov = siteCov;
        }

        public int[][] getAllelesSortedByFrequency(int site) {
            int offset = site - startSite;
            int alleleCnt = 0;
            while (myAlleleFreqOrder[alleleCnt][offset] != GenotypeTable.UNKNOWN_ALLELE) {
                alleleCnt++;
            }
            int result[][] = new int[2][alleleCnt];
            for (int i = 0; i < alleleCnt; i++) {
                result[0][i] = myAlleleFreqOrder[i][offset];
                result[1][i] = myAlleleCnt[result[0][i]][offset];
            }
            return result;
        }

        public float getMAF(int site) {
            return maf[site - startSite];
        }

        public float getSiteCoverage(int site) {
            return siteCov[site - startSite];
        }

    }

    private static long getCacheKey(int taxon, int site) {
        return ((long) taxon << 33) + (site / HDF5_GENOTYPE_BLOCK_SIZE);
    }

    private static int getTaxonFromKey(long key) {
        return (int) (key >>> 33);
    }

    private static int getSiteStartFromKey(long key) {
        return (int) ((key << 33) >>> 33);
    }

    private static int getStartSite(int site) {
        return site & SITE_BLOCK_MASK;
    }

    private String getTaxaGenoPath(int taxon) {
        return genotypePaths[taxon];
    }

    private HDF5ByteGenotypeCallTable(IHDF5Reader reader, int numTaxa, int numSites, boolean phased, String[][] alleleEncodings) {
        super(numTaxa, numSites, phased, alleleEncodings);
        genotypePaths = new String[numTaxa];
        TaxaList tL = new TaxaListBuilder().buildFromHDF5Genotypes(reader);  //not the most efficient thing to do, but ensures sort is the same.
        for (int i = 0; i < numTaxa; i++) {
            genotypePaths[i] = Tassel5HDF5Constants.getGenotypesCallsPath(tL.taxaName(i));
        }
        myHDF5Reader = reader;
        long oneThirdMemory = Runtime.getRuntime().maxMemory() / 196608l;
        long oneColumnBlockForEachProcess = numTaxa * Runtime.getRuntime().availableProcessors();
        myGenoCache = CacheBuilder.newBuilder()
                .maximumSize(Math.min(oneThirdMemory, oneColumnBlockForEachProcess))
                .build(myGenoLoader);
        mySiteAnnoCache = CacheBuilder.newBuilder()
                .maximumSize(150)
                .build(siteAnnotLoader);
    }

    static HDF5ByteGenotypeCallTable getInstance(IHDF5Reader reader) {
        if (!HDF5Utils.isHDF5GenotypeLocked(reader)) {
            throw new IllegalStateException("The Genotype module of this HDF5 file hasn't been locked, and therefore can't be opened for reading. This could occur if the file was created using the -ko (keep open) option when running the plugin ProductionSNPCallerPluginV2. Please check your file, close if appropriate, and try again.");
        }
        int numTaxa = reader.int32().getAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA);
        int numSites = reader.int32().getAttr(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES);
        String[][] alleleEncodings = NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;
        return new HDF5ByteGenotypeCallTable(reader, numTaxa, numSites, false, alleleEncodings);
    }

    @Override
    public byte genotype(int taxon, int site) {
        try {
            byte[] data = myGenoCache.get(getCacheKey(taxon, site));
            return data[site % HDF5_GENOTYPE_BLOCK_SIZE];
        } catch (ExecutionException ex) {
            myLogger.error(ex.getMessage(), ex);
            throw new IllegalStateException("HDF5ByteGenotyeCallTable: getBase: Error getting base from cache: " + ex.getMessage());
        }
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(genotype(taxon, site));
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        try {
            SiteBlockAttr sa = mySiteAnnoCache.get(getStartSite(site));
            return sa.getAllelesSortedByFrequency(site);
        } catch (ExecutionException e) {
            e.printStackTrace();
            throw new UnsupportedOperationException("Error in getMinorAlleleFrequency from cache");
        }
    }

    @Override
    public double minorAlleleFrequency(int site) {
        try {
            SiteBlockAttr sa = mySiteAnnoCache.get(getStartSite(site));
            return sa.getMAF(site);
        } catch (ExecutionException e) {
            e.printStackTrace();
            throw new UnsupportedOperationException("Error in getMinorAlleleFrequency from cache");
        }
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isSiteOptimized() {
        return false;
    }

    @Override
    public Stream<Byte> stream() {
        return StreamSupport.stream(spliterator(), true);
    }

    @Override
    public Stream<Byte> stream(int taxon) {
        return StreamSupport.stream(new HDF5ByteGenotypeCallTableSpliterator<>(taxon, 0, numberOfSites(), taxon, numberOfSites()), true);
    }

    @Override
    public Spliterator<Byte> spliterator() {
        return new HDF5ByteGenotypeCallTableSpliterator<>(0, 0, numberOfSites(), numberOfTaxa() - 1, numberOfSites());
    }

    class HDF5ByteGenotypeCallTableSpliterator<T extends Byte> extends AbstractGenotypeCallTableSpliterator<Byte> {

        HDF5ByteGenotypeCallTableSpliterator(int taxaOrigin, int siteOrigin, int numSites, int taxaFence, int siteFence) {
            super(taxaOrigin, siteOrigin, numSites, taxaFence, siteFence);
        }

        @Override
        public void forEachRemaining(Consumer<? super Byte> action) {

            for (; myTaxaOrigin < myTaxaFence; myTaxaOrigin++) {
                while (mySiteOrigin < myNumSites) {
                    try {
                        byte[] data = myGenoCache.get(getCacheKey(myTaxaOrigin, mySiteOrigin));
                        int startIndex = mySiteOrigin % HDF5_GENOTYPE_BLOCK_SIZE;
                        int endIndex = data.length;
                        for (int i = startIndex; i < endIndex; i++) {
                            action.accept(data[i]);
                        }
                        mySiteOrigin += endIndex - startIndex;
                    } catch (ExecutionException ex) {
                        myLogger.error(ex.getMessage(), ex);
                        throw new IllegalStateException("HDF5ByteGenotyeCallTable: HDF5ByteGenotypeCallTableSpliterator: forEachRemaining: Error getting base from cache: " + ex.getMessage());
                    }
                }
                mySiteOrigin = 0;
            }
            while (mySiteOrigin < mySiteFence) {
                try {
                    byte[] data = myGenoCache.get(getCacheKey(myTaxaOrigin, mySiteOrigin));
                    int startIndex = mySiteOrigin % HDF5_GENOTYPE_BLOCK_SIZE;
                    int endIndex = Math.min(data.length, mySiteFence - mySiteOrigin);
                    for (int i = startIndex; i < endIndex; i++) {
                        action.accept(data[i]);
                    }
                    mySiteOrigin += endIndex - startIndex;
                } catch (ExecutionException ex) {
                    myLogger.error(ex.getMessage(), ex);
                    throw new IllegalStateException("HDF5ByteGenotyeCallTable: HDF5ByteGenotypeCallTableSpliterator: forEachRemaining: Error getting base from cache: " + ex.getMessage());
                }
            }
        }

        @Override
        public Spliterator<Byte> trySplit() {
            long size = estimateSize();
            if (size > HDF5_GENOTYPE_BLOCK_SIZE) {
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
                midSite -= midSite % HDF5_GENOTYPE_BLOCK_SIZE;
                myTaxaOrigin = midTaxa;
                mySiteOrigin = midSite;
                return new HDF5ByteGenotypeCallTableSpliterator<>(loTaxa, loSite, myNumSites, midTaxa, midSite);
            } else {
                return null;
            }
        }

    }
}
