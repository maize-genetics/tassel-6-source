/*
 *  DistanceMatrixWithCounts
 * 
 *  Created on Jan 14, 2016
 */
package net.maizegenetics.taxa.distance;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.GeneralAnnotation;

/**
 *
 * @author Terry Casstevens
 */
public class DistanceMatrixWithCounts extends DistanceMatrix {

    private final int[][] myCounts;

    DistanceMatrixWithCounts(float[][] distances, TaxaList taxa, GeneralAnnotation annotations, int[][] counts) {
        super(distances, taxa, annotations);
        myCounts = counts;
    }

    DistanceMatrixWithCounts(double[][] distance, TaxaList taxaList, GeneralAnnotation annotations, int[][] counts) {
        super(distance, taxaList, annotations);
        myCounts = counts;
    }

    public int getCount(int x, int y) {
        if (x > y) {
            return myCounts[x][y];
        } else {
            return myCounts[y][x];
        }
    }

}
