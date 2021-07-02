/*
 *  GenotypeBuilder
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;

/**
 * Builder to construct a GenotypeCallTable. This builder is generally only used
 * in complex situations, where the GenotypeTableBuilder does not suffice.
 *
 * @author Terry Casstevens
 */
public class GenotypeCallTableBuilder {

    private static Logger myLogger = LogManager.getLogger(GenotypeCallTableBuilder.class);

    private SuperByteMatrix myGenotype;
    private boolean myIsPhased = false;
    private String[][] myAlleleEncodings = NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;

    private GenotypeCallTableBuilder(SuperByteMatrix genotype) {
        myGenotype = genotype;
    }

    /**
     * Get Genotype Builder given number of taxa and sites. Performance
     * optimized for site loop inside taxon loop. Default is unphased and
     * NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES encoding.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites.
     *
     * @return Genotype Builder
     */
    public static GenotypeCallTableBuilder getInstance(int numTaxa, int numSites) {
        return getUnphasedNucleotideGenotypeBuilder(numTaxa, numSites);
    }

    /**
     * Get Genotype Builder given number of taxa and sites. Performance
     * optimized for site loop inside taxon loop. Default is unphased and
     * NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES encoding.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites.
     *
     * @return Genotype Builder
     */
    public static GenotypeCallTableBuilder getUnphasedNucleotideGenotypeBuilder(int numTaxa, int numSites) {
        SuperByteMatrix matrix = SuperByteMatrixBuilder.getInstance(numTaxa, numSites);
        matrix.setAll(GenotypeTable.UNKNOWN_GENOTYPE);
        return new GenotypeCallTableBuilder(matrix);
    }

    public static GenotypeCallTableBuilder getInstanceCopy(GenotypeCallTable genotype) {
        if (genotype instanceof ByteGenotypeCallTable) {
            SuperByteMatrix matrix = SuperByteMatrixBuilder.getInstanceCopy(((ByteGenotypeCallTable) genotype).myGenotype);
            return new GenotypeCallTableBuilder(matrix).isPhased(genotype.isPhased()).alleleEncodings(genotype.alleleDefinitions());
        } else {
            final int NUM_TAXA_TO_COPY = 10;
            int numTaxa = genotype.numberOfTaxa();
            int numSites = genotype.numberOfSites();
            GenotypeCallTableBuilder builder = GenotypeCallTableBuilder.getInstance(numTaxa, numSites).isPhased(genotype.isPhased()).alleleEncodings(genotype.alleleDefinitions());
            int numThreads = Runtime.getRuntime().availableProcessors();
            ExecutorService pool = Executors.newFixedThreadPool(numThreads);
            List<Future<?>> futures = new ArrayList<>();
            for (int t = 0; t < numTaxa; t += NUM_TAXA_TO_COPY) {
                int numTaxaToCopy = Math.min(NUM_TAXA_TO_COPY, numTaxa - t);
                Future<?> future = pool.submit(new CopyAllSitesFromTaxa(genotype, builder, t, numTaxaToCopy));
                futures.add(future);
            }

            for (Future<?> future : futures) {
                try {
                    future.get();
                } catch (InterruptedException | ExecutionException e) {
                    myLogger.debug(e.getMessage(), e);
                    throw new IllegalStateException("GenotypeCallTableBuilder: getInstanceCopy: " + e.getMessage());
                }
            }
            pool.shutdown();

            return builder;
        }
    }

    private static class CopyAllSitesFromTaxa implements Runnable {

        private final GenotypeCallTable mySrc;
        private final GenotypeCallTableBuilder myDest;
        private final int myStartTaxon;
        private final int myNumTaxaToCopy;
        private final int myNumSites;

        public CopyAllSitesFromTaxa(GenotypeCallTable src, GenotypeCallTableBuilder dest, int startTaxon, int numTaxaToCopy) {
            mySrc = src;
            myDest = dest;
            myStartTaxon = startTaxon;
            myNumTaxaToCopy = numTaxaToCopy;
            myNumSites = src.numberOfSites();
        }

        @Override
        public void run() {
            for (int t = myStartTaxon, n = myStartTaxon + myNumTaxaToCopy; t < n; t++) {
                for (int s = 0; s < myNumSites; s++) {
                    myDest.setBase(t, s, mySrc.genotype(t, s));
                }
            }
        }
    }

    public GenotypeCallTableBuilder setBase(int taxon, int site, byte value) {
        myGenotype.set(taxon, site, value);
        return this;
    }

    public GenotypeCallTableBuilder setBaseRangeForTaxon(int taxon, int startSite, byte[] value) {
        myGenotype.arraycopy(taxon, value, startSite);
        return this;
    }

    public GenotypeCallTableBuilder setBases(String[] data) {

        int numTaxa = data.length;

        int numSites = data[0].length();

        for (int site = 0; site < numSites; site++) {
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                setBase(taxon, site, NucleotideAlignmentConstants.getNucleotideDiploidByte(data[taxon].charAt(site)));
            }
        }

        return this;

    }

    public GenotypeCallTableBuilder setBases(String[][] data) {

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("BitAlignment: getInstance: data can not be empty.");
        }

        int numTaxa = data.length;
        int numSites = data[0].length;

        for (int site = 0; site < numSites; site++) {
            if (data[0][0].contains(":")) {
                Pattern colon = Pattern.compile(":");
                for (int taxon = 0; taxon < numTaxa; taxon++) {
                    if (data[taxon][site].equalsIgnoreCase(GenotypeTable.UNKNOWN_GENOTYPE_STR)) {
                        setBase(taxon, site, GenotypeTable.UNKNOWN_GENOTYPE);
                    } else if (data[taxon][site].equals("?") || data[taxon][site].equals("?:?")) {
                        setBase(taxon, site, GenotypeTable.UNKNOWN_GENOTYPE);
                    } else {
                        String[] siteval = colon.split(data[taxon][site]);
                        byte first = NucleotideAlignmentConstants.getNucleotideAlleleByte(siteval[0]);
                        byte second = NucleotideAlignmentConstants.getNucleotideAlleleByte(siteval[1]);
                        setBase(taxon, site, (byte) ((first << 4) | second));
                    }
                }
            } else {
                for (int taxon = 0; taxon < numTaxa; taxon++) {
                    if (data[taxon][site].equalsIgnoreCase(GenotypeTable.UNKNOWN_ALLELE_STR)) {
                        setBase(taxon, site, GenotypeTable.UNKNOWN_GENOTYPE);
                    } else if (data[taxon][site].equals("?")) {
                        setBase(taxon, site, GenotypeTable.UNKNOWN_GENOTYPE);
                    } else {
                        setBase(taxon, site, NucleotideAlignmentConstants.getNucleotideAlleleByte(data[taxon][site]));
                    }
                }
            }
        }
        return this;
    }

    public GenotypeCallTableBuilder isPhased(boolean isPhased) {
        myIsPhased = isPhased;
        return this;
    }

    public GenotypeCallTableBuilder alleleEncodings(String[][] alleleEncodings) {
        myAlleleEncodings = alleleEncodings;
        return this;
    }

    public int getTaxaCount() {
        return myGenotype.getNumRows();
    }

    public int getSiteCount() {
        return myGenotype.getNumColumns();
    }

    public void reorderTaxa(int[] newIndices) {
        myGenotype.reorderRows(newIndices);
    }

    public void reorderPositions(int[] newIndices) {
        myGenotype.reorderColumns(newIndices);
    }

    public GenotypeCallTable build() {
        SuperByteMatrix temp = myGenotype;
        myGenotype = null;
        if (NucleotideAlignmentConstants.isNucleotideEncodings(myAlleleEncodings)) {
            return new NucleotideGenotypeCallTable(temp, myIsPhased);
        } else {
            return new ByteGenotypeCallTable(temp, myIsPhased, myAlleleEncodings);
        }
    }
}
