package net.maizegenetics.dna.snp.io;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;

import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.*;

/**
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class BuilderFromPLINK {

    private static final Logger myLogger = Logger.getLogger(BuilderFromPLINK.class);
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    private static final int NUM_PLINK_NON_SITE_HEADERS = 6;

    private final String myPedFile;
    private final String myMapFile;
    private final ProgressListener myProgressListener;
    private boolean mySortTaxaAlphabetically = false;
    private boolean mySortPositions = false;

    private BuilderFromPLINK(String pedfile, String mapfile, ProgressListener listener) {
        myPedFile = pedfile;
        myMapFile = mapfile;
        myProgressListener = listener;
    }

    public static BuilderFromPLINK getBuilder(String pedfile, String mapfile, ProgressListener listener) {
        return new BuilderFromPLINK(pedfile, mapfile, listener);
    }

    public GenotypeTable build() {

        GenotypeTable result = null;
        try {
            int numThreads = Runtime.getRuntime().availableProcessors();
            ExecutorService pool = Executors.newFixedThreadPool(numThreads);

            myLogger.info("Reading: " + myPedFile + " and " + myMapFile);

            PositionListBuilder posBuild = processSites(myMapFile);
            myLogger.info("Number of sites: " + posBuild.size());

            int numOfTaxa = Utils.getNumberLinesNotHashOrBlank(myPedFile);
            myLogger.info("Number of taxa: " + numOfTaxa);

            GenotypeCallTableBuilder genotypeCallTableBuilder = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(numOfTaxa, posBuild.size());

            int linesAtTime = (int) Math.ceil((1 << 24) / posBuild.size());
            ArrayList<String> textLines = new ArrayList<>(linesAtTime);

            List<Future<ProcessPLINKBlock>> futures = new ArrayList<>();
            int numLines = 0;
            try (BufferedReader reader = Utils.getBufferedReader(myPedFile)) {
                String currLine = reader.readLine();
                while (currLine != null) {
                    textLines.add(currLine);
                    numLines++;
                    if (numLines % linesAtTime == 0) {
                        ProcessPLINKBlock processBlock = new ProcessPLINKBlock(textLines, genotypeCallTableBuilder, numLines - linesAtTime, numLines * 100 / numOfTaxa);
                        futures.add(pool.submit(processBlock));
                        textLines = new ArrayList<>(linesAtTime);
                    }
                    currLine = reader.readLine();
                }
            }

            if (textLines.size() > 0) {
                ProcessPLINKBlock processBlock = new ProcessPLINKBlock(textLines, genotypeCallTableBuilder, numLines - textLines.size(), 100);
                futures.add(pool.submit(processBlock));
            }

            TaxaListBuilder taxaBuild = new TaxaListBuilder();
            for (Future<ProcessPLINKBlock> future : futures) {
                ProcessPLINKBlock pb = future.get();
                taxaBuild.addAll(pb.getBlockTaxa());
            }
            TaxaList taxaList = taxaBuild.build();

            pool.shutdown();

            // Check that result is in correct order. If not, either try to sort
            // or just throw an error (determined by what was passed to fullSort)
            if (posBuild.validateOrdering() == false) {
                if (mySortPositions) {
                    posBuild.sortPositions(genotypeCallTableBuilder);
                    // Double-check post-sort ordering. Should never happen, but just to be safe
                    if (posBuild.validateOrdering() == false) {
                        throw new IllegalStateException("BuilderFromPLINK: Ordering of PLINK failed.");
                    }
                } else {
                    throw new IllegalStateException("BuilderFromPLINK: Ordering incorrect. PLINK must be ordered by position.");
                }
            }
            GenotypeCallTable g = genotypeCallTableBuilder.build();
            result = GenotypeTableBuilder.getInstance(g, posBuild.build(), taxaList);

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("BuilderFromPLINK: build: problem processing: " + e.getMessage());
        }

        return result;
    }

    /**
     * Set the builder so that when built it will sort the taxa
     */
    public BuilderFromPLINK sortTaxa() {
        mySortTaxaAlphabetically = true;
        return this;
    }

    /**
     * Set the builder so that when built it will sort the positions.
     */
    public BuilderFromPLINK sortPositions() {
        mySortPositions = true;
        return this;
    }

    // chromosome (1-22, X, Y or 0 if unplaced)
    // rs# or snp identifier
    // Genetic distance (morgans)
    // Base-pair position (bp units)
    private static final int PLINK_MAP_CHROMOSOME_INDEX = 0;
    private static final int PLINK_MAP_SND_ID_INDEX = 1;
    private static final int PLINK_MAP_GENETIC_DISTANCE_INDEX = 2;
    private static final int PLINK_MAP_POSITION_INDEX = 3;
    private static final int NUM_PLINK_MAP_COLUMNS = 4;

    private static PositionListBuilder processSites(String mapfile) {

        Map<String, Chromosome> chromosomes = new HashMap<>();
        List<Position> positions = new ArrayList<>();

        int lineNum = 1;
        try (BufferedReader reader = Utils.getBufferedReader(mapfile)) {
            String line = reader.readLine();
            while (line != null) {
                String[] tokens = WHITESPACE_PATTERN.split(line);
                if (tokens.length < NUM_PLINK_MAP_COLUMNS) {
                    throw new IllegalStateException("BuilderFromPLINK: processSites: Not all columns defined line : \"" + line + "\" of file: " + mapfile);
                }
                Chromosome chr = chromosomes.get(tokens[PLINK_MAP_CHROMOSOME_INDEX]);
                if (chr == null) {
                    chr = Chromosome.instance(new String(tokens[PLINK_MAP_CHROMOSOME_INDEX]));
                    chromosomes.put(tokens[PLINK_MAP_CHROMOSOME_INDEX], chr);
                }
                GeneralPosition current = new GeneralPosition.Builder(chr, Integer.parseInt(tokens[PLINK_MAP_POSITION_INDEX]))
                        .snpName(new String(tokens[PLINK_MAP_SND_ID_INDEX])).build();
                positions.add(current);
                line = reader.readLine();
                lineNum++;
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("BuilderFromPLINK: processSites: problem with: " + mapfile + " line: " + lineNum);
        }

        PositionListBuilder result = new PositionListBuilder();
        result.addAll(positions);
        return result;

    }

    // Family ID
    // Individual ID
    // Paternal ID
    // Maternal ID
    // Sex (1=male; 2=female; other=unknown)
    // Phenotype
    private static final int PLINK_PED_FAMILY_ID_INDEX = 0;
    private static final int PLINK_PED_INDIVIDUAL_ID_INDEX = 1;
    private static final int PLINK_PED_PATERNAL_ID_INDEX = 2;
    private static final int PLINK_PED_MATERNAL_ID_INDEX = 3;
    private static final int PLINK_PED_SEX_INDEX = 4;
    private static final int PLINK_PED_PHENOTYPE_INDEX = 5;

    private class ProcessPLINKBlock implements Callable<ProcessPLINKBlock> {

        private final int myNumTaxaToProcess;
        private ArrayList<String> myTextLines;
        private final ArrayList<Taxon> myBlockTaxaList;
        private final GenotypeCallTableBuilder myBuilder;
        private final int myStartTaxon;
        private final int myProgress;

        private ProcessPLINKBlock(ArrayList<String> textLines, GenotypeCallTableBuilder builder, int startTaxon, int progress) {
            myNumTaxaToProcess = textLines.size();
            myTextLines = textLines;
            myBlockTaxaList = new ArrayList<>(myNumTaxaToProcess);
            myBuilder = builder;
            myStartTaxon = startTaxon;
            myProgress = progress;
        }

        @Override
        public ProcessPLINKBlock call() throws Exception {

            for (int t = 0; t < myNumTaxaToProcess; t++) {
                String input = myTextLines.get(t);
                try {
                    String[] tokens = WHITESPACE_PATTERN.split(input, NUM_PLINK_NON_SITE_HEADERS + 1);

                    String taxonName = new String(tokens[PLINK_PED_INDIVIDUAL_ID_INDEX].trim());
                    Taxon taxon = new Taxon.Builder(taxonName).build();
                    myBlockTaxaList.add(taxon);
                    int taxonIndex = t + myStartTaxon;
                    for (int i = 0, n = tokens[NUM_PLINK_NON_SITE_HEADERS].length(); i < n; i += 4) {
                        myBuilder.setBase(taxonIndex, i / 4, GenotypeTableUtils.getDiploidValue(getPLINKAlleleByte(tokens[NUM_PLINK_NON_SITE_HEADERS].charAt(i)),
                                getPLINKAlleleByte(tokens[NUM_PLINK_NON_SITE_HEADERS].charAt(i + 2))));
                    }
                } catch (ArrayIndexOutOfBoundsException e) {
                    myLogger.debug(e.getMessage(), e);
                    throw new IllegalStateException("BuilderFromPLINK: ProcessPLINKBlock: problem: to many genotypes in file: " + myPedFile + " line: " + input.substring(0, Math.min(input.length(), 50)));
                } catch (Exception e) {
                    myLogger.debug(e.getMessage(), e);
                    throw new IllegalStateException("BuilderFromPLINK: ProcessPLINKBlock: problem: " + e.getMessage() + " file: " + myPedFile + " line: " + input.substring(0, Math.min(input.length(), 50)));
                }
            }
            myTextLines = null;
            if (myProgressListener != null) {
                myProgressListener.progress(myProgress, null);
            }

            return this;

        }

        List<Taxon> getBlockTaxa() {
            return myBlockTaxaList;
        }

    }

    private static final Map<String, Byte> PLINK_ALLELE_HASH = new HashMap<>();
    private static final byte[] PLINK_ALLELE_ARRAY = new byte[256];

    static {
        PLINK_ALLELE_HASH.put("A", A_ALLELE);
        PLINK_ALLELE_HASH.put("C", C_ALLELE);
        PLINK_ALLELE_HASH.put("G", G_ALLELE);
        PLINK_ALLELE_HASH.put("T", T_ALLELE);
        PLINK_ALLELE_HASH.put("+", INSERT_ALLELE);
        PLINK_ALLELE_HASH.put("-", GAP_ALLELE);
        PLINK_ALLELE_HASH.put("N", GenotypeTable.UNKNOWN_ALLELE);
        PLINK_ALLELE_HASH.put("0", GenotypeTable.UNKNOWN_ALLELE);
        PLINK_ALLELE_HASH.put("1", A_ALLELE);
        PLINK_ALLELE_HASH.put("2", C_ALLELE);
        PLINK_ALLELE_HASH.put("3", G_ALLELE);
        PLINK_ALLELE_HASH.put("4", T_ALLELE);
        Arrays.fill(PLINK_ALLELE_ARRAY, UNDEFINED_ALLELE);
        for (Map.Entry<String, Byte> en : PLINK_ALLELE_HASH.entrySet()) {
            PLINK_ALLELE_ARRAY[en.getKey().charAt(0)] = en.getValue();
        }
    }

    /**
     * Returns haploid byte value for given PLINK value. Only right-most four
     * bits used.
     *
     * @param value haploid allele value
     *
     * @return nucleotide haploid allele byte value
     */
    public static byte getPLINKAlleleByte(char value) {
        try {
            return PLINK_ALLELE_ARRAY[value];
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("BuilderFromPLINK: getPLINKAlleleByte: unknown allele value: " + value);
        }
    }

}
