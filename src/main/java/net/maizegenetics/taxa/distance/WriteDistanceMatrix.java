package net.maizegenetics.taxa.distance;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Map;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 * @author Terry Casstevens
 * @author Zack Miller
 */
public class WriteDistanceMatrix {

    private static final Logger myLogger = Logger.getLogger(WriteDistanceMatrix.class);

    private WriteDistanceMatrix() {
        //utility
    }

    public static void saveDelimitedDistanceMatrix(DistanceMatrix matrix, String saveFile) {

        if ((saveFile == null) || (saveFile.isEmpty())) {
            throw new IllegalArgumentException("WriteDistanceMatrix: saveDelimitedDistanceMatrix: No file specified.");
        }

        try (BufferedWriter bw = Utils.getBufferedWriter(saveFile)) {

            GeneralAnnotation annotations = matrix.annotations();
            if (annotations != null) {
                for (Map.Entry<String, String> current : annotations.getAllAnnotationEntries()) {
                    bw.write("##");
                    bw.write(current.getKey());
                    bw.write("=");
                    bw.write(current.getValue());
                    bw.write("\n");
                }
            }

            bw.write(String.valueOf(matrix.getRowCount()));
            bw.write("\n");

            for (long r = 0, n = matrix.getRowCount(); r < n; r++) {
                Object[] theRow = matrix.getRow(r);
                for (int i = 0; i < theRow.length; i++) {
                    if (i != 0) {
                        bw.write("\t");
                    }
                    bw.write(theRow[i].toString());
                }
                bw.write("\n");
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("WriteDistanceMatrix: saveDelimitedDistanceMatrix: problem writing file: " + e.getMessage());
        }

        if (matrix instanceof DistanceMatrixWithCounts) {
            String[] grmFilenames = DistanceMatrixUtils.getGRMFilenames(saveFile);
            String grmNBinFilename = grmFilenames[2];
            saveBinMultiBlupCounts(grmNBinFilename, (DistanceMatrixWithCounts) matrix);
        }

    }

    public static void saveRawMultiBlupMatrix(DistanceMatrix matrix, String filename) {
        String[] grmFilenames = DistanceMatrixUtils.getGRMFilenames(filename);
        saveRawMultiBlupMatrix(matrix, grmFilenames[0], grmFilenames[3]);
    }

    public static void saveRawMultiBlupMatrix(DistanceMatrix matrix, String taxaFile, String matrixFile) {

        if (matrixFile == null || taxaFile == null) {
            return;
        }

        saveMultiBlupIDs(taxaFile, matrix.getTaxaList());

        try (BufferedWriter bw = Utils.getBufferedWriter(matrixFile)) {

            for (long r = 0, n = matrix.getRowCount(); r < n; r++) {
                Object[] theRow = matrix.getRow(r);

                for (int i = 1; i < theRow.length; i++) {
                    if (i != 1) {
                        bw.write("\t");
                    }
                    bw.write(theRow[i].toString());
                }
                bw.write("\n");
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            myLogger.error("WriteDistanceMatrix: saveDelimitedDistanceMatrix: problem writing file: " + matrixFile);
        }

    }

    public static void saveBinMultiBlupMatrix(DistanceMatrix matrix, String filename) {
        String[] grmFilenames = DistanceMatrixUtils.getGRMFilenames(filename);
        saveBinMultiBlupMatrix(matrix, grmFilenames[0], grmFilenames[1], grmFilenames[2]);
    }

    public static void saveBinMultiBlupMatrix(DistanceMatrix matrix, String taxaFile, String matrixFile, String countFile) {

        if (matrixFile == null || taxaFile == null || countFile == null) {
            return;
        }

        saveMultiBlupIDs(taxaFile, matrix.getTaxaList());
        saveBinMultiBlupCounts(countFile, matrix);

        try (BufferedOutputStream bw = Utils.getBufferedOutputStream(matrixFile)) {

            for (long r = 0, n = matrix.getRowCount(); r < n; r++) {
                Object[] theRow = matrix.getRow(r);

                for (int i = 1; i < r + 2; i++) {

                    ByteBuffer kinsBuffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);

                    kinsBuffer.putFloat(Float.parseFloat(theRow[i].toString()));
                    bw.write(kinsBuffer.array());

                }
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            myLogger.error("WriteDistanceMatrix: saveDelimitedDistanceMatrix: problem writing file: " + matrixFile + ".  " + e.getMessage());
        }

    }

    public static void saveMultiBlupIDs(String filename, TaxaList taxa) {

        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {
            for (Taxon taxon : taxa) {
                writer.write(taxon.getName());
                writer.write("\t");
                writer.write(taxon.getName());
                writer.write("\n");
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("WriteDistanceMatrix: saveMultiBlupIDs: Problem writing: " + filename + ".  " + e.getMessage());
        }

    }

    public static void saveBinMultiBlupCounts(String filename, DistanceMatrix matrix) {

        if (matrix instanceof DistanceMatrixWithCounts) {

            DistanceMatrixWithCounts matrixWCounts = (DistanceMatrixWithCounts) matrix;
            int numTaxa = matrixWCounts.numberOfTaxa();

            try (BufferedOutputStream writer = Utils.getBufferedOutputStream(filename)) {
                ByteBuffer countsBuffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
                for (int taxaOne = 0; taxaOne < numTaxa; taxaOne++) {
                    for (int taxaTwo = 0; taxaTwo <= taxaOne; taxaTwo++) {
                        countsBuffer.clear();
                        countsBuffer.putFloat((float) matrixWCounts.getCount(taxaOne, taxaTwo));
                        writer.write(countsBuffer.array());
                    }
                }
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("WriteDistanceMatrix: saveBinMultiBlupCounts: Problem writing: " + filename + ".  " + e.getMessage());
            }

        } else {

            int numTaxa = matrix.numberOfTaxa();

            try (BufferedOutputStream writer = Utils.getBufferedOutputStream(filename)) {
                ByteBuffer countsBuffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
                countsBuffer.putFloat(1.0f);
                byte[] counts = countsBuffer.array();
                for (int taxaOne = 0; taxaOne < numTaxa; taxaOne++) {
                    for (int taxaTwo = 0; taxaTwo <= taxaOne; taxaTwo++) {
                        writer.write(counts);
                    }
                }
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("WriteDistanceMatrix: saveBinMultiBlupCounts: Problem writing: " + filename + ".  " + e.getMessage());
            }

        }

    }

    public static void saveDARwinMatrix(DistanceMatrix matrix, String saveFile) {
        String[] filenames = DistanceMatrixUtils.getDARwinFilenames(saveFile);
        saveDARwinMatrix(matrix, filenames[0], filenames[1]);
    }

    public static void saveDARwinMatrix(DistanceMatrix matrix, String donFile, String disFile) {

        if ((disFile == null) || (donFile == null)) {
            throw new IllegalArgumentException("WriteDistanceMatrix: saveDARwinMatrix: No file specified.");
        }

        saveDARwinIDs(donFile, matrix.getTaxaList());

        int numTaxa = matrix.numberOfTaxa();

        try (BufferedWriter writer = Utils.getBufferedWriter(disFile)) {

            writer.write("@DARwin 5.0 - DIS\n");
            writer.write(String.valueOf(numTaxa));
            writer.write("\n");

            for (int t = 0; t < numTaxa - 1; t++) {
                writer.write("\t");
                writer.write(String.valueOf(t + 1));
            }
            writer.write("\n");

            for (int t = 1; t < numTaxa; t++) {
                writer.write(String.valueOf(t + 1));
                for (int t2 = 0; t2 < t; t2++) {
                    writer.write("\t");
                    writer.write(String.valueOf(matrix.getDistance(t, t2)));
                }
                writer.write("\n");
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("WriteDistanceMatrix: saveDARwinMatrix: problem writing: " + disFile);
        }
        
        myLogger.info("saveDARwinIDs: wrote file: " + disFile);

    }

    public static void saveDARwinIDs(String filename, TaxaList taxa) {

        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {
            writer.write("@DARwin 5.0 - DON\n");
            writer.write(String.valueOf(taxa.numberOfTaxa()));
            writer.write("\t1\n");
            writer.write("Unit\tTaxa\n");
            int index = 1;
            for (Taxon taxon : taxa) {
                writer.write(String.valueOf(index++));
                writer.write("\t");
                writer.write(taxon.getName());
                writer.write("\n");
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("WriteDistanceMatrix: saveDARwinIDs: Problem writing: " + filename + ".  " + e.getMessage());
        }
        
        myLogger.info("saveDARwinIDs: wrote file: " + filename);

    }

}
