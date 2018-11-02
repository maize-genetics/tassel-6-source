//
// ReadDistanceMatrix.java
//
package net.maizegenetics.taxa.distance;

import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;
import net.maizegenetics.util.GeneralAnnotationStorage;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;

import java.io.IOException;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import java.util.regex.Pattern;

import org.apache.log4j.Logger;

/**
 * @author Terry Casstevens
 */
public class ReadDistanceMatrix {

    private static final Logger myLogger = Logger.getLogger(ReadDistanceMatrix.class);

    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s+");

    private ReadDistanceMatrix() {
    }

    public static DistanceMatrix readDistanceMatrix(String filename) {

        String[] grmFilenames = DistanceMatrixUtils.getGRMFilenames(filename);
        String grmNBinFilename = grmFilenames[2];

        try (BufferedReader reader = Utils.getBufferedReader(filename)) {

            GeneralAnnotationStorage.Builder annotations = GeneralAnnotationStorage.getBuilder();
            boolean notFinished = true;
            String line = null;
            while (notFinished) {
                line = reader.readLine();
                if (line == null) {
                    throw new IllegalArgumentException("ReadDistanceMatrix: There is no matrix data in this file: " + filename);
                }
                line = line.trim();
                if (line.isEmpty()) {
                    // do nothing
                } else if (line.startsWith("##")) {
                    String[] keyValue = line.substring(2).split("=", 2);
                    annotations.addAnnotation(keyValue[0].trim(), keyValue[1].trim());
                } else {
                    notFinished = false;
                }
            }

            int numTaxa = 0;
            try {
                numTaxa = Integer.parseInt(line);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("ReadDistanceMatrix: The number of taxa is not a number: " + line);
            }

            DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(numTaxa);
            builder.annotation(annotations.build());

            String current = reader.readLine();
            int index = 0;
            while (current != null) {

                if (index >= numTaxa) {
                    throw new IllegalArgumentException("ReadDistancMatrix: There are too many lines in this file.  Expected: " + (numTaxa + 1) + " counting the first line (number of taxa)");
                }

                String[] tokens = WHITESPACE_PATTERN.split(current);
                if (tokens.length != numTaxa + 1) {
                    throw new IllegalStateException("ReadDistanceMatrix: Incorrect number of values on line number: " + (index + 2) + " expected: " + (numTaxa + 1) + " counting taxon name. actual: " + tokens.length);
                }
                builder.addTaxon(new Taxon(tokens[0]));

                for (int i = 0; i < numTaxa; i++) {
                    try {
                        builder.set(index, i, Double.parseDouble(tokens[i + 1]));
                    } catch (NumberFormatException nfex) {
                        myLogger.debug(nfex.getMessage(), nfex);
                        throw new IllegalArgumentException("ReadDistanceMatrix: Incorrectly formatted number: " + tokens[i + 1] + " on line number: " + (index + 2));
                    }
                }

                current = reader.readLine();
                index++;
            }

            if (index != numTaxa) {
                throw new IllegalArgumentException("ReadDistanceMatrix: There are too few lines in this file.  Expected: " + (numTaxa + 1) + " counting the first line (number of taxa)");
            }

            if (new File(grmNBinFilename).exists()) {
                readBinMultiBlupCounts(grmNBinFilename, builder);
            }

            return builder.build();

        } catch (IOException ioex) {
            myLogger.debug(ioex.getMessage(), ioex);
            throw new IllegalStateException("ReadDistanceMatrix: Problem reading file: " + filename);
        }

    }

    public static DistanceMatrix readBinMultiBlupMatrix(String filename) {

        String[] grmFilenames = DistanceMatrixUtils.getGRMFilenames(filename);
        String grmIdFilename = grmFilenames[0];
        String grmBinFilename = grmFilenames[1];
        String grmNBinFilename = grmFilenames[2];

        TaxaList taxa = readBinMultiBlupID(grmIdFilename);
        int numTaxa = taxa.numberOfTaxa();

        DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(taxa);
        int taxaOne = 0;
        int taxaTwo = 0;
        try (InputStream input = Utils.getInputStream(grmBinFilename)) {
            byte[] current = new byte[4];
            int numRead = input.read(current);
            while (numRead == 4) {
                float f = ByteBuffer.wrap(current).order(ByteOrder.LITTLE_ENDIAN).getFloat();
                builder.set(taxaTwo, taxaOne++, f);
                if (taxaOne > taxaTwo) {
                    taxaOne = 0;
                    taxaTwo++;
                }
                numRead = input.read(current);
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("ReadDistanceMatrix: readBinMultiBlupMatrix: Problem reading: " + grmBinFilename);
        }

        readBinMultiBlupCounts(grmNBinFilename, builder);

        return builder.build();

    }

    public static TaxaList readBinMultiBlupID(String filename) {

        TaxaListBuilder builder = new TaxaListBuilder();

        try (BufferedReader reader = Utils.getBufferedReader(filename)) {
            String current = reader.readLine();
            while (current != null) {
                String[] tokens = current.split("\t");
                if (tokens.length != 2) {
                    throw new IllegalStateException("ReadDistanceMatrix: readBinMultiBlupID: Rows must have two columns: " + filename);
                }
                builder.add(new Taxon(tokens[1]));
                current = reader.readLine();
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("ReadDistanceMatrix: readBinMultiBlupID: Problem reading: " + filename);
        }

        return builder.build();

    }

    public static void readBinMultiBlupCounts(String filename, DistanceMatrixBuilder builder) {

        try (InputStream input = Utils.getInputStream(filename)) {
            int taxaOne = 0;
            int taxaTwo = 0;
            byte[] current = new byte[4];
            int numRead = input.read(current);
            while (numRead == 4) {
                float f = ByteBuffer.wrap(current).order(ByteOrder.LITTLE_ENDIAN).getFloat();
                builder.setCount(taxaTwo, taxaOne++, Math.round(f));
                if (taxaOne > taxaTwo) {
                    taxaOne = 0;
                    taxaTwo++;
                }
                numRead = input.read(current);
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("ReadDistanceMatrix: readBinMultiBlupCounts: Problem reading: " + filename);
        }

    }

}
