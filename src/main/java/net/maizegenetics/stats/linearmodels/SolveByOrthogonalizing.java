package net.maizegenetics.stats.linearmodels;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.FDistribution;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.matrixalgebra.decomposition.SingularValueDecomposition;

import org.apache.commons.math3.exception.OutOfRangeException;

public class SolveByOrthogonalizing {
    private List<ModelEffect> myBaseModel;
    private List<double[]> myBasisVectors;
    private List<double[]> myData;
    private List<double[]> myOrthogonalizedData;
    private List<double[]> UColumns = null;
    private SingularValueDecomposition baseSvd = null;
    private final static double tol = 1e-10;

    private SolveByOrthogonalizing() {

    }

    public static SolveByOrthogonalizing getInstanceFromModel(List<ModelEffect> baseModel, List<double[]> dataList) {
        SolveByOrthogonalizing sbo = new SolveByOrthogonalizing();
        sbo.myBaseModel = baseModel;
        sbo.myData = dataList;
        DoubleMatrix[][] design = sbo.createDesignMatricesFromModel();
        sbo.computeBaseSvd(design);
        sbo.OrthogonalizeData();
        return sbo;
    }

    public static SolveByOrthogonalizing getInstanceFromVectors(List<double[]> basisVectors, List<double[]> dataList) {
        SolveByOrthogonalizing sbo = new SolveByOrthogonalizing();
        sbo.myBasisVectors = basisVectors;
        sbo.myData = dataList;
        sbo.computeBaseSvd(sbo.createDesignMatricesFromVectors());
        sbo.OrthogonalizeData();
        return sbo;
    }

    public SolveByOrthogonalizing.Marker solveForR(SolveByOrthogonalizing.Marker marker) {
        if (marker.vector2 == null)
            return solveForR(marker.position(), marker.vector1());
        return solveForR(marker.position(), marker.vector1(), marker.vector2);
    }

    public SolveByOrthogonalizing.Marker solveForR(Position pos, double[] values) {

        double[] val1 = center(values);
        double[] val2 = orthogonalizeByBase(val1);
        double[] val3 = centerAndScale(val2);
        double[] orthogonalValues = centerAndScale(orthogonalizeByBase(center(values)));

        
        if (orthogonalValues == null) {
            double[] rValues =
                    IntStream.range(0, myOrthogonalizedData.size()).mapToDouble(i -> Double.NaN).toArray();
            return new SolveByOrthogonalizing.Marker(pos, rValues, 0);
        }

        double[] rValues = new double[myOrthogonalizedData.size()];
        int count = 0;
        for (double[] data : myOrthogonalizedData) {
            double ip = innerProduct(data, orthogonalValues);
            rValues[count++] = ip * ip;
        }
        return new SolveByOrthogonalizing.Marker(pos, rValues, 1);
    }

    public SolveByOrthogonalizing.Marker solveForR(Position pos, double[] add, double[] dom) {
        if (dom == null)
            return solveForR(pos, add);
        double[] orthogonalAdd = orthogonalizeByBase(center(add));
        double[] orthogonalDom = orthogonalizeByBase(center(dom));

        //orthogonalize dom with respect to add
        double mult =
                innerProduct(orthogonalAdd, orthogonalDom) / innerProduct(orthogonalAdd, orthogonalAdd);
        int n = orthogonalDom.length;
        for (int i = 0; i < n; i++) {
            orthogonalDom[i] -= mult * orthogonalAdd[i];
        }

        //center and scale
        double[] v1 = centerAndScale(orthogonalAdd);
        double[] v2 = centerAndScale(orthogonalDom);

        if (v1 == null) {
            if (v2 == null) {
                double[] rValues =
                        IntStream.range(0, myOrthogonalizedData.size()).mapToDouble(i -> Double.NaN).toArray();
                return new SolveByOrthogonalizing.Marker(pos, rValues, 0);
            }
            double[] rValues =
                    myOrthogonalizedData.stream().mapToDouble(d -> innerProduct(d, v2)).map(d -> d
                            * d).toArray();
            return new SolveByOrthogonalizing.Marker(pos, rValues, 1);
        }
        if (v2 == null) {
            double[] rValues =
                    myOrthogonalizedData.stream().mapToDouble(d -> innerProduct(d, v1)).map(d -> d
                            * d).toArray();
            return new SolveByOrthogonalizing.Marker(pos, rValues, 1);
        }

        double[] rValues = myOrthogonalizedData.stream()
                .mapToDouble(d -> {
                    double r1 = innerProduct(v1, d);
                    double r2 = innerProduct(v2, d);
                    return r1 * r1 + r2 * r2;
                }).toArray();

        return new SolveByOrthogonalizing.Marker(pos, rValues, 2);
    }

    private DoubleMatrix[][] createDesignMatricesFromModel() {
        DoubleMatrix[][] designMatrices = new DoubleMatrix[1][];
        designMatrices[0] = myBaseModel.stream()
                .filter(a -> !a.getID().toString().toLowerCase().equals("mean"))
                .map(me -> me.getX())
                .toArray(DoubleMatrix[]::new);
        return designMatrices;
    }

    private DoubleMatrix[][] createDesignMatricesFromVectors() {
        DoubleMatrix[][] designMatrices = new DoubleMatrix[1][];
        designMatrices[0] = myBasisVectors.stream()
                .map(d -> DoubleMatrixFactory.DEFAULT.make(d.length, 1, d))
                .toArray(DoubleMatrix[]::new);
        return designMatrices;
    }

    private void computeBaseSvd(DoubleMatrix[][] designMatrices) {
        if (designMatrices[0].length == 0 || designMatrices[0][0] == null)
            return;

        DoubleMatrix X = DoubleMatrixFactory.DEFAULT.compose(designMatrices);

        //center the columns of X
        int nrows = X.numberOfRows();
        double dblnrows = nrows;
        int ncols = X.numberOfColumns();
        for (int c = 0; c < ncols; c++) {
            double mean = X.columnSum(c) / dblnrows;
            for (int r = 0; r < nrows; r++) {
                X.set(r, c, X.get(r, c) - mean);
            }
        }

        baseSvd = X.getSingularValueDecomposition();
        DoubleMatrix U = baseSvd.getU(false);
        UColumns = new ArrayList<>();
        int ncol = U.numberOfColumns();
        for (int i = 0; i < ncol; i++) {
            UColumns.add(U.column(i).to1DArray());
        }

    }

    private void OrthogonalizeData() {
        if (baseSvd == null) {
            myOrthogonalizedData = myData.stream()
                    .map(d -> centerAndScale(Arrays.copyOf(d, d.length)))
                    .collect(Collectors.toList());

        } else {
            myOrthogonalizedData = myData.stream()
                    .map(d -> center(Arrays.copyOf(d, d.length)))
                    .map(d -> orthogonalizeByBase(d))
                    .map(d -> centerAndScale(d))
                    .collect(Collectors.toList());
        }

    }

    private double[] orthogonalizeByBase(double[] vector) {
        if (baseSvd == null) {
            return Arrays.copyOf(vector, vector.length);
        }

        int nrows = vector.length;
        double[] result = Arrays.copyOf(vector, nrows);

        DoubleMatrix U = baseSvd.getU(false);
        int ncol = U.numberOfColumns();
        for (int i = 0; i < ncol; i++) {
            double[] u = U.column(i).to1DArray();
            double ip = innerProduct(vector, u);
            for (int j = 0; j < nrows; j++)
                result[j] -= ip * u[j];
        }

        return result;
    }

    /**
     * @return the df in the base model, including 1 df for the mean
     */
    public int baseDf() {
        int df = 1;  //the mean
        if (baseSvd != null) {
            double[] singularValues = baseSvd.getSingularValues();

            int n = singularValues.length;
            for (int i = 0; i < n; i++) {
                if (singularValues[i] > tol)
                    df++;
            }
        }
        return df;
    }

    public static double innerProduct(double[] x, double[] y) {
        int n = x.length;
        double sumprod = 0;
        for (int i = 0; i < n; i++) {
            sumprod += x[i] * y[i];
        }
        return sumprod;
    }

    public static double[] center(double[] values) {
        int n = values.length;
        double mean = Arrays.stream(values).sum() / n;
        for (int i = 0; i < n; i++)
            values[i] -= mean;
        return values;
    }

    public static double[] scale(double[] values) {
        int n = values.length;
        double divisor = Math.sqrt(innerProduct(values, values));
        for (int i = 0; i < n; i++)
            values[i] /= divisor;
        return values;
    }

    public static double[] centerAndScale(double[] values) {
        int n = values.length;
        double sum = 0;
        double sumsq = 0;
        for (int i = 0; i < n; i++) {
            double val = values[i];
            sum += val;
        }
        double mean = sum / n;

        for (int i = 0; i < n; i++) {
            values[i] = values[i] - mean;
            sumsq += values[i] * values[i];
        }

        double divisor = Math.sqrt(sumsq);
        if (divisor < tol)
            return null;
        for (int i = 0; i < n; i++)
            values[i] /= divisor;
        return values;
    }

    public static double calculateFfromR2(double r2, double markerDf, double errorDf) {
        return r2 / (1 - r2) * errorDf / markerDf;
    }

    public static double calculateP(double F, double markerDf, double errorDf) {
        if (!Double.isFinite(F))
            return Double.NaN;
        double p;
        try {
            p = LinearModelUtils.Ftest(F, markerDf, errorDf);
        } catch (Exception e) {
            p = Double.NaN;
        }
        return p;
    }

    public static double calculateR2Fromp(double alpha, double modelDf, double errorDf) {
        //returns the value of R^2 corresponding to the value of F, f for which P(F>f) = alpha

        FDistribution fdist = new FDistribution(modelDf, errorDf);
        try {
            double p = 1 - alpha;
            double F = fdist.inverseCumulativeProbability(p);
            double Fme = F * modelDf / errorDf;
            return Fme / (1 + Fme);
        } catch (OutOfRangeException e) {
            e.printStackTrace();
            return Double.NaN;
        }

    }

    public List<double[]> getOrthogonalizedData() {
        return myOrthogonalizedData;
    }
    
    public List<double[]> copyOrthogonalizedData() {
        return myOrthogonalizedData.stream().map(darray -> Arrays.copyOf(darray, darray.length)).collect(Collectors.toList());
    }
    
    public List<double[]> getUColumns() {
        return UColumns;
    }
    
    public List<double[]> copyUColumns() {
        return UColumns.stream().map(darray -> Arrays.copyOf(darray, darray.length)).collect(Collectors.toList());
    }
    
    public static class Marker {
        public final Position myPosition;
        public final double[] vector1;
        public final double[] vector2;
        public int df;

        public Marker(Position pos, double[] values, int df) {
            myPosition = pos;
            vector1 = values;
            vector2 = null;
            this.df = df;
        }

        public Marker(Position pos, double[] additive, double[] dominant, int df) {
            myPosition = pos;
            vector1 = additive;
            vector2 = dominant;
            this.df = df;
        }

        public Position position() {
            return myPosition;
        }

        public double[] vector1() {
            return vector1;
        }

        public double[] vector2() {
            return vector2;
        }

        public boolean hasTwoVectors() {
            return vector2 != null;
        }

        public int degreesOfFreedom() {
            if (vector2 == null)
                return 1;
            return 2;
        }
    }

}
