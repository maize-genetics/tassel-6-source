/*
 * ImportUtils
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.io.BuilderFromHapMap;
import net.maizegenetics.dna.snp.io.BuilderFromPLINK;
import net.maizegenetics.dna.snp.io.BuilderFromVCF;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * Methods for importing GenotypeTables from various file formats.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class ImportUtils {

    private static final Logger myLogger = Logger.getLogger(ImportUtils.class);

    private ImportUtils() {
        // Utility Class - do not instantiate.
    }

    /**
     * @deprecated Use ImportUtils.read()
     */
    public static GenotypeTable readGuessFormat(String fileName) {

        if (fileName.endsWith("hmp.txt.gz") || fileName.endsWith("hmp.txt")) {
            return readFromHapmap(fileName, null);
        } else if (fileName.endsWith(".vcf") || fileName.endsWith(".vcf.gz")) {
            return readFromVCF(fileName, null);
        }
        return null;

    }

    public static GenotypeTable readFromVCF(final String filename, ProgressListener listener, boolean keepDepth, boolean sortPositions) {
        BuilderFromVCF builder = BuilderFromVCF.getBuilder(filename, listener);
        if (keepDepth) {
            builder.keepDepth();
        }
        if (sortPositions) {
            return builder.buildAndSortInMemory();
        } else {
            return builder.build();
        }
    }

    public static GenotypeTable readFromVCF(final String filename, ProgressListener listener, boolean keepDepth) {
        return readFromVCF(filename, listener, keepDepth, false);
    }

    public static GenotypeTable readFromVCF(final String filename, ProgressListener listener) {
        return readFromVCF(filename, listener, true);
    }

    public static GenotypeTable readFromVCF(final String filename) {
        return readFromVCF(filename, null);
    }

    /**
     * Read GenotypeTable from HapMap file
     *
     * @param filename input HapMap file name
     *
     * @return a genotype table
     */
    public static GenotypeTable readFromHapmap(final String filename) {
        return readFromHapmap(filename, null);
    }

    /**
     * Read GenotypeTable from HapMap file
     *
     * @param filename input HapMap file name
     * @param listener progress listener to track reading rate
     *
     * @return a genotype table
     */
    public static GenotypeTable readFromHapmap(final String filename, ProgressListener listener) {
        return readFromHapmap(filename, listener, false);
    }

    /**
     * Read GenotypeTable from HapMap file
     *
     * @param filename input HapMap file name
     * @param listener progress listener to track reading rate
     * @param sortPositions whether to sort positions
     *
     * @return
     */
    public static GenotypeTable readFromHapmap(final String filename, ProgressListener listener, boolean sortPositions) {
        BuilderFromHapMap builder = BuilderFromHapMap.getBuilder(filename, listener);
        if (sortPositions) {
            builder.sortPositions();
        }
        return builder.build();
    }

    public static GenotypeTable readFromPLink(final String pedFilename, final String mapFilename, ProgressListener listener) {
        return BuilderFromPLINK.getBuilder(pedFilename, mapFilename, listener).build();
    }

    public static GenotypeTable readFromPLink(final String pedFilename, final String mapFilename, ProgressListener listener, boolean sortPositions) {
        if (sortPositions) {
            return BuilderFromPLINK.getBuilder(pedFilename, mapFilename, listener).sortPositions().build();
        } else {
            return BuilderFromPLINK.getBuilder(pedFilename, mapFilename, listener).build();
        }
    }

    public static GenotypeTable readFasta(String filename) throws FileNotFoundException, IOException {

        BufferedReader reader = Utils.getBufferedReader(filename);

        List<String> taxa = new ArrayList<>();
        List<String> sequences = new ArrayList<>();

        String line = reader.readLine();
        boolean sequence = false;
        int sequenceLength = -1;
        int count = 1;
        while (line != null) {

            line = line.trim();

            if (line.startsWith(";")) {
                line = reader.readLine();
            } else if (line.startsWith(">")) {
                StringTokenizer tokens = new StringTokenizer(line);
                String taxaName = tokens.nextToken();
                if (taxaName.length() == 1) {
                    taxaName = tokens.nextToken();
                } else {
                    taxaName = taxaName.substring(1).trim();
                }
                taxa.add(taxaName);
                sequence = true;
                line = reader.readLine();
            } else if (sequence) {
                StringBuilder builder = new StringBuilder();
                while ((line != null) && (!line.startsWith(">")) && (!line.startsWith(";"))) {
                    line = line.trim().toUpperCase();
                    builder.append(line);
                    line = reader.readLine();
                }
                String temp = builder.toString();
                if (sequenceLength == -1) {
                    sequenceLength = temp.length();
                } else if (sequenceLength != temp.length()) {
                    throw new IllegalStateException("ImportUtils: readFasta: Sequence: " + count + " Differs in Length.");
                }
                sequences.add(temp);
                sequence = false;
                count++;
            } else {
                myLogger.error("readFasta: file: " + filename + " invalid format.");
                throw new IllegalArgumentException("Import: readFasta: invalid format.");
            }

        }

        String[] taxaNames = new String[taxa.size()];
        taxa.toArray(taxaNames);
        TaxaList taxaList = (new TaxaListBuilder()).addAll(taxaNames).build();

        String[] sequenceArray = new String[sequences.size()];
        sequences.toArray(sequenceArray);

        GenotypeCallTable genotype = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(sequences.size(), sequenceLength)
                .setBases(sequenceArray)
                .build();

        return GenotypeTableBuilder.getInstance(genotype, PositionListBuilder.getInstance(sequenceLength), taxaList);

    }
}
