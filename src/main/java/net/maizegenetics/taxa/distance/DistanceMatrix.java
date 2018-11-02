// DistanceMatrix.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.taxa.distance;

import net.maizegenetics.util.FormattedOutput;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

/**
 * Storage for pairwise distance matrices. Only stores half the matrix as it is
 * symmetrical.<p>
 *
 * For best performance, iterate over matrix this way.
 * <blockquote><pre>
 * DistanceMatrix matrix;
 * for (int i = 0; i < myNumTaxa; i++) {
 *     for (int j = 0; j <= i; j++) {
 *         matrix.getDistance(i, j);
 *     }
 * }
 * </pre></blockquote>
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 * @author Terry Casstevens
 */
public class DistanceMatrix implements TableReport {

    private final TaxaList myTaxaList;
    private final int myNumTaxa;
    private final GeneralAnnotation myAnnotations;
    private final float[][] myDistances;

    /**
     * Use DistanceMatrixBuilder instead of this.
     *
     * @see DistanceMatrixBuilder
     */
    DistanceMatrix(float[][] distances, TaxaList taxa, GeneralAnnotation annotations) {
        myDistances = distances;
        myTaxaList = taxa;
        myNumTaxa = myTaxaList.numberOfTaxa();
        myAnnotations = annotations;
    }

    /**
     * Constructor taking distances array and taxa list. Use
     * DistanceMatrixBuilder instead of this.
     *
     * @see DistanceMatrixBuilder
     */
    public DistanceMatrix(double[][] distance, TaxaList taxa) {
        this(distance, taxa, null);
    }

    /**
     * Use DistanceMatrixBuilder instead of this.
     *
     * @see DistanceMatrixBuilder
     */
    public DistanceMatrix(double[][] distances, TaxaList taxa, GeneralAnnotation annotations) {
        myNumTaxa = taxa.numberOfTaxa();
        if ((distances == null) || (distances.length != myNumTaxa) || (distances[0].length != myNumTaxa)) {
            throw new IllegalArgumentException("DistanceMatrix: init: dimensions of distances aren't correct.");
        }
        myDistances = new float[myNumTaxa][];
        for (int i = 0; i < myNumTaxa; i++) {
            myDistances[i] = new float[i + 1];
        }
        for (int x = 0; x < myNumTaxa; x++) {
            for (int y = 0; y <= x; y++) {
                if (Math.abs(distances[x][y] - distances[y][x]) > 0.0000001) {
                    throw new IllegalStateException("DistanceMatrix: init: values passed in are not symmetrical: " + distances[x][y] + " and: " + distances[y][x]);
                }
                myDistances[x][y] = (float) distances[x][y];
            }
        }
        myTaxaList = taxa;
        myAnnotations = annotations;
    }

    /**
     * Constructor that clones a distance matrix.
     */
    public DistanceMatrix(DistanceMatrix dm) {
        myNumTaxa = dm.numberOfTaxa();
        myDistances = new float[myNumTaxa][];
        for (int i = 0; i < myNumTaxa; i++) {
            myDistances[i] = new float[i + 1];
        }
        for (int x = 0; x < myNumTaxa; x++) {
            for (int y = 0; y <= x; y++) {
                myDistances[x][y] = dm.myDistances[x][y];
            }
        }
        myTaxaList = dm.myTaxaList;
        myAnnotations = dm.myAnnotations;
    }

    /**
     * Constructor that clones a distance matrix and for only the specified
     * taxa.
     */
    public DistanceMatrix(DistanceMatrix dm, TaxaList subset) {

        myNumTaxa = subset.numberOfTaxa();
        myDistances = new float[myNumTaxa][];
        for (int i = 0; i < myNumTaxa; i++) {
            myDistances[i] = new float[i + 1];
        }

        for (int i = 0; i < myNumTaxa; i++) {
            int index1 = dm.whichIdNumber(subset.taxaName(i));
            myDistances[i][i] = dm.myDistances[index1][index1];
            for (int j = 0; j < i; j++) {
                int index2 = dm.whichIdNumber(subset.taxaName(j));
                myDistances[i][j] = dm.getDistance(index1, index2);
            }
        }
        myTaxaList = subset;
        myAnnotations = dm.myAnnotations;

    }

    /**
     * print alignment (PHYLIP format)
     */
    public void printPHYLIP(PrintWriter out) throws IOException {
        // PHYLIP header line
        out.println("  " + myNumTaxa);
        FormattedOutput format = FormattedOutput.getInstance();

        for (int i = 0; i < myNumTaxa; i++) {
            format.displayLabel(out,
                    myTaxaList.taxaName(i), 10);
            out.print("      ");

            for (int j = 0; j < myNumTaxa; j++) {
                // Chunks of 6 blocks each
                if (j % 6 == 0 && j != 0) {
                    out.println();
                    out.print("                ");
                }

                out.print("  ");
                format.displayDecimal(out, getDistance(i, j), 5);
            }
            out.println();
        }
    }

    /**
     * returns representation of this alignment as a string
     */
    @Override
    public String toString() {

        StringWriter sw = new StringWriter();
        try {
            printPHYLIP(new PrintWriter(sw));
        } catch (Exception e) {
            e.printStackTrace();
        }

        return sw.toString();

    }

    /**
     * compute squared distance to second distance matrix
     */
    public double squaredDistance(DistanceMatrix mat, boolean weighted) {
        double sum = 0;
        for (int i = 0; i < myNumTaxa - 1; i++) {
            for (int j = 0; j < i; j++) {
                double diff = myDistances[i][j] - mat.getDistance(i, j);
                double weight;
                if (weighted) {
                    // Fitch-Margoliash weight
                    // (variances proportional to distances)
                    float distance = myDistances[i][j];
                    weight = 1.0 / distance * distance;
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
     * compute absolute distance to second distance matrix
     */
    public double absoluteDistance(DistanceMatrix mat) {
        double sum = 0;
        for (int i = 0; i < myNumTaxa - 1; i++) {
            for (int j = 0; j < i; j++) {
                double diff = Math.abs(myDistances[i][j] - mat.getDistance(i, j));
                sum += diff;
            }
        }

        return 2.0 * sum; // we counted only half the matrix
    }

    /**
     * Returns the number of taxa which is also the number of rows and columns
     * that the distance matrix has.
     */
    public int getSize() {
        return myNumTaxa;
    }

    /**
     * Returns the distances as a 2-dimensional array of doubles. Matrix is
     * cloned first so it can be altered freely.
     */
    public final double[][] getClonedDistances() {
        double[][] copy = new double[myNumTaxa][myNumTaxa];
        for (int i = 0; i < myNumTaxa; i++) {
            for (int j = 0; j <= i; j++) {
                copy[i][j] = myDistances[i][j];
                copy[j][i] = copy[i][j];
            }
        }
        return copy;
    }

    /**
     * Returns the distances as a 2-dimensional array of doubles (in the actual
     * array used to store the distances)
     */
    public final double[][] getDistances() {
        return getClonedDistances();
    }

    public final float getDistance(final int row, final int col) {
        if (row > col) {
            return myDistances[row][col];
        } else {
            return myDistances[col][row];
        }
    }

    /**
     * Returns the mean pairwise distance of this matrix
     */
    public double meanDistance() {
        double dist = 0.0;
        int count = 0;
        for (int i = 1; i < myNumTaxa; i++) {
            for (int j = 0; j < i; j++) {
                float distance = myDistances[i][j];
                if (!Float.isNaN(distance)) {
                    dist += distance;
                    count++;
                }
            }
        }
        return dist / (double) count;
    }

    public Taxon getTaxon(int i) {
        return myTaxaList.get(i);
    }

    public int numberOfTaxa() {
        return myTaxaList.numberOfTaxa();
    }

    public int whichIdNumber(String name) {
        return myTaxaList.indexOf(name);
    }

    public int whichIdNumber(Taxon id) {
        return myTaxaList.indexOf(id);
    }

    /**
     * Return TaxaList of this alignment.
     */
    public TaxaList getTaxaList() {
        return myTaxaList;
    }

    /**
     * test whether this matrix is a symmetric distance matrix
     *
     */
    public boolean isSymmetric() {
        for (int i = 0; i < myNumTaxa; i++) {
            if (myDistances[i][i] != 0) {
                return false;
            }
        }
        return true;
    }

    private boolean isIn(int value, int[] set) {
        if (set == null) {
            return false;
        }
        for (int i = 0; i < set.length; i++) {
            if (set[i] == value) {
                return true;
            }
        }
        return false;
    }

    /**
     * @param fromIndex the index of the thing (taxa,sequence) from which we
     * want to find the closest (excluding self)
     * @param exclusion indexes of things that should not be considered, may be
     * null
     * @return the index of the member closes to the specified
     */
    public int getClosestIndex(int fromIndex, int[] exclusion) {
        float min = Float.POSITIVE_INFINITY;
        int index = -1;
        for (int i = 0; i < myNumTaxa; i++) {
            if (i != fromIndex && !isIn(i, exclusion)) {
                float d = getDistance(fromIndex, i);
                if (d < min) {
                    min = d;
                    index = i;
                }
            }
        }
        return index;
    }

    public static DistanceMatrix hadamardProduct(DistanceMatrix m0, DistanceMatrix m1) {

        int n = m0.numberOfTaxa();
        if (m1.numberOfTaxa() != n) {
            throw new IllegalArgumentException("Matrices must be of the same dimensions to compute a Hadamard product.");
        }

        DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(m0.getTaxaList());
        for (int r = 0; r < n; r++) {
            for (int c = 0; c <= r; c++) {
                builder.set(r, c, m0.myDistances[r][c] * m1.myDistances[r][c]);
            }
        }

        return builder.build();

    }

    @Override
    public Object[] getTableColumnNames() {
        String[] colNames = new String[getSize() + 1];
        colNames[0] = "Taxa";
        for (int i = 0; i < myNumTaxa; i++) {
            colNames[i + 1] = getTaxon(i).toString();
        }
        return colNames;
    }

    /**
     * Returns specified row.
     *
     * @param rowLong row number
     *
     * @return row
     */
    @Override
    public Object[] getRow(long rowLong) {

        int row = (int) rowLong;
        Object[] result = new Object[myNumTaxa + 1];
        result[0] = getTaxon(row);
        for (int j = 1; j <= myNumTaxa; j++) {
            result[j] = String.valueOf(getDistance(row, j - 1));
        }

        return result;

    }

    @Override
    public String getTableTitle() {
        return "Distance Matrix";
    }

    @Override
    public long getRowCount() {
        return myNumTaxa;
    }

    @Override
    public long getElementCount() {
        return getRowCount() * getColumnCount();
    }

    @Override
    public int getColumnCount() {
        return myNumTaxa + 1;
    }

    @Override
    public Object getValueAt(long rowIndex, int columnIndex) {
        if (columnIndex == 0) {
            return getTaxon((int) rowIndex);
        }
        return getDistance((int) rowIndex, columnIndex - 1);
    }

    public String getColumnName(int col) {
        if (col == 0) {
            return "Taxa";
        }
        return getTaxon(col - 1).toString();
    }

    public GeneralAnnotation annotations() {
        return myAnnotations;
    }

}
