/*
 *  LineIndexHapmapGenotypeCallTable
 * 
 *  Created on Aug 23, 2015
 */
package net.maizegenetics.dna.snp.genotypecall;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import htsjdk.samtools.util.BlockCompressedInputStream;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.io.LineIndex;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CopyOnWriteArraySet;
import java.util.concurrent.ForkJoinPool;

/**
 * @author Terry Casstevens
 */
public class LineIndexHapmapGenotypeCallTable extends AbstractGenotypeCallTable {

    private static final Logger myLogger = Logger.getLogger(LineIndexHapmapGenotypeCallTable.class);
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;
    private static final int NUM_LOOK_AHEAD_BLOCKS = 103;

    private final String myFilename;
    private final LineIndex myIndex;
    private final boolean myIsOneLetter;
    private final int myNumLinesPerInterval;
    private final ConcurrentLinkedQueue<BlockCompressedInputStream> myReaders = new ConcurrentLinkedQueue<>();
    private final CopyOnWriteArraySet<Integer> myCurrentlyProcessingBlocks = new CopyOnWriteArraySet<>();

    private final Cache<Integer, byte[][]> myGenoCache;

    private final Cache<Integer, byte[]> mySmallGenoCache;

    private final ForkJoinPool myThreadPool;

    private LineIndexHapmapGenotypeCallTable(int numTaxa, int numSites, boolean phased, boolean isOneLetter, LineIndex index, String filename) {
        super(numTaxa, numSites, phased, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        myIsOneLetter = isOneLetter;
        myIndex = index;
        myNumLinesPerInterval = index.numLinesPerInterval();
        myFilename = filename;

        long oneThirdMemory = Runtime.getRuntime().maxMemory() / (numTaxa * myNumLinesPerInterval * 3);
        int maxCacheSize = (int) Math.min((long) (110 * Runtime.getRuntime().availableProcessors()), oneThirdMemory);

        myGenoCache = CacheBuilder.newBuilder()
                .initialCapacity(maxCacheSize)
                .maximumSize(maxCacheSize)
                .build();

        mySmallGenoCache = CacheBuilder.newBuilder()
                .initialCapacity(1000)
                .maximumSize(1000)
                .build();

        myThreadPool = ForkJoinPool.commonPool();

    }

    public static LineIndexHapmapGenotypeCallTable getInstance(int numTaxa, int numSites, boolean phased, boolean isOneLetter, LineIndex index, String filename) {
        return new LineIndexHapmapGenotypeCallTable(numTaxa, numSites, phased, isOneLetter, index, filename);
    }

    private byte[] getFromCache(int site) {

        int blockNumber = site / myNumLinesPerInterval;

        byte[][] result = myGenoCache.getIfPresent(blockNumber);

        if (result == null) {

            if (myCurrentlyProcessingBlocks.add(blockNumber + 1)) {
                myThreadPool.submit(new ProcessLines(blockNumber + 1));
            }

            result = parseBlock(site);

        }

        return result[site % myNumLinesPerInterval];

    }

    private BlockCompressedInputStream getReader() {
        BlockCompressedInputStream reader = myReaders.poll();
        if (reader == null) {
            try {
                reader = new BlockCompressedInputStream(new File(myFilename));
            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
            }
        }
        return reader;
    }

    @Override
    public byte genotype(int taxon, int site) {
        try {
            byte[] result = mySmallGenoCache.getIfPresent(site);
            if (result != null) {
                return result[taxon];
            } else {
                result = getFromCache(site);
                mySmallGenoCache.put(site, result);
                return result[taxon];
            }
        } catch (Exception ex) {
            myLogger.error(ex.getMessage(), ex);
            throw new IllegalStateException("LineIndexHapmapGenotypeCallTable: genotype: Error getting genotype from cache: " + ex.getMessage());
        }
    }

    @Override
    public byte[] genotypeForAllTaxa(int site) {
        byte[] result = new byte[myTaxaCount];
        System.arraycopy(getFromCache(site), 0, result, 0, myTaxaCount);
        return result;
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
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isSiteOptimized() {
        return true;
    }

    /**
     * Parse line from Hapmap file to genotypes for a site.
     *
     * @param input input line
     * @param numTaxa number of taxa
     * @param site site
     * @param isOneLetter is genotypes code as one letter or two
     *
     * @return genotypes
     */
    private static byte[] parseLine(String input, int numTaxa, int site, boolean isOneLetter) {

        int len = input.length();
        int tabIndex = 0;
        int offset = 0;
        for (int i = 0; (tabIndex < NUM_HAPMAP_NON_TAXA_HEADERS) && (i < len); i++) {
            if (input.charAt(i) == '\t') {
                tabIndex++;
                offset = i + 1;
            }
        }

        byte[] data = new byte[numTaxa];
        int taxon = 0;
        if (isOneLetter) {
            for (int i = offset; i < len; i += 2) {
                if (taxon >= numTaxa) {
                    throw new IllegalStateException("LineIndexHapmapGenotypeCallTable: Site: " + site + " has too many values.");
                }
                byte value = NucleotideAlignmentConstants.getNucleotideDiploidByte(input.charAt(i));
                if (value == NucleotideAlignmentConstants.UNDEFINED_HOMOZYGOUS) {
                    throw new IllegalStateException("LineIndexHapmapGenotypeCallTable: Site: " + site + " has illegal value: " + input.charAt(i));
                }
                data[taxon++] = value;
            }
        } else {
            for (int i = offset; i < len; i += 3) {
                if (taxon >= numTaxa) {
                    throw new IllegalStateException("LineIndexHapmapGenotypeCallTable: Site: " + site + " has too many values.");
                }
                // there is a phasing conflict with the existing import approach
                byte value = GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.getNucleotideDiploidByte(input.charAt(i + 1)),
                        NucleotideAlignmentConstants.getNucleotideDiploidByte(input.charAt(i)));
                if (value == NucleotideAlignmentConstants.UNDEFINED_HOMOZYGOUS) {
                    throw new IllegalStateException("LineIndexHapmapGenotypeCallTable: Site: " + site + " has illegal value: " + input.charAt(i) + input.charAt(i + 1));
                }
                data[taxon++] = value;
            }
        }

        return data;
    }

    private byte[] parseLine(int site) {

        BlockCompressedInputStream reader = null;
        try {

            int seekIndex = site / myNumLinesPerInterval;

            reader = getReader();
            reader.seek(myIndex.virtualOffset(seekIndex));

            int siteIndex = site % myNumLinesPerInterval;
            for (int i = 0; i < siteIndex; i++) {
                reader.readLine();
            }
            return parseLine(reader.readLine(), myTaxaCount, site, myIsOneLetter);

        } catch (Exception e) {
            myLogger.error(e.getMessage(), e);
            throw new IllegalStateException("LineIndexHapmapGenotypeCallTable: parseBlock: problem parsing site: " + site);
        } finally {
            myReaders.add(reader);
        }

    }

    private byte[][] parseBlock(int site) {

        BlockCompressedInputStream reader = null;
        try {

            int processBlock = site / myNumLinesPerInterval;
            int startSite = processBlock * myNumLinesPerInterval;
            int seekIndex = startSite / myNumLinesPerInterval;

            if (startSite >= mySiteCount) {
                return null;
            }

            reader = getReader();
            reader.seek(myIndex.virtualOffset(seekIndex));

            return parseBlock(processBlock, startSite, reader);

        } catch (Exception e) {
            myLogger.error(e.getMessage(), e);
            throw new IllegalStateException("LineIndexHapmapGenotypeCallTable: parseBlock: problem parsing site: " + site);
        } finally {
            myReaders.add(reader);
        }

    }

    private byte[][] parseBlock(int processBlock, int startSite, BlockCompressedInputStream reader) {

        try {

            int numSites = Math.min(myNumLinesPerInterval, mySiteCount - startSite);
            byte[][] result = new byte[numSites][];
            for (int i = 0; i < numSites; i++) {
                result[i] = parseLine(reader.readLine(), myTaxaCount, startSite + i, myIsOneLetter);
            }
            myGenoCache.put(processBlock, result);
            // This get to prevent early eviction from cache
            myGenoCache.getIfPresent(processBlock);
            myCurrentlyProcessingBlocks.remove(processBlock);

            return result;

        } catch (Exception e) {
            myLogger.error(e.getMessage(), e);
            throw new IllegalStateException("LineIndexHapmapGenotypeCallTable: parseBlock: problem parsing site: " + processBlock);
        }

    }

    private class ProcessLines implements Runnable {

        private int myStartSite;
        private final int mySeekIndex;
        private final int myProcessBlock;

        public ProcessLines(int block) {
            myProcessBlock = block;
            myStartSite = myProcessBlock * myNumLinesPerInterval;
            mySeekIndex = myStartSite / myNumLinesPerInterval;
        }

        @Override
        public void run() {

            if (myStartSite >= mySiteCount) {
                return;
            }

            BlockCompressedInputStream reader = getReader();
            try {

                reader.seek(myIndex.virtualOffset(mySeekIndex));

                parseBlock(myProcessBlock, myStartSite, reader);

                myStartSite += myNumLinesPerInterval;
                if (myStartSite >= mySiteCount) {
                    return;
                }

                for (int b = 1; b < NUM_LOOK_AHEAD_BLOCKS; b++) {

                    if (myGenoCache.getIfPresent(myProcessBlock + b) != null) {
                        return;
                    }
                    if (!myCurrentlyProcessingBlocks.add(myProcessBlock + b)) {
                        return;
                    }

                    parseBlock(myProcessBlock + b, myStartSite, reader);

                    myStartSite += myNumLinesPerInterval;
                    if (myStartSite >= mySiteCount) {
                        return;
                    }

                }

            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
                throw new IllegalStateException("LineIndexHapmapGenotypeCallTable: ProcessLines: problem parsing site: " + myStartSite);
            } finally {
                myReaders.add(reader);
            }

        }

    }

}
