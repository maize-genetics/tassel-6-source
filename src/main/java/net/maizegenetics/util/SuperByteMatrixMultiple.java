/*
 *  SuperByteMatrixMultiple
 */
package net.maizegenetics.util;

import java.util.Arrays;
import java.util.Spliterator;
import static java.util.Spliterator.IMMUTABLE;
import static java.util.Spliterator.ORDERED;
import static java.util.Spliterator.SIZED;
import static java.util.Spliterator.SUBSIZED;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrixMultiple implements SuperByteMatrix {

    private final byte[][] myData;
    private final int myNumRows;
    private final int myNumColumns;
    private final int myNumRowsPerSingleDimArray;
    private final long myNumElements;
    private final int myNumElementsPerSingleDimArray;

    SuperByteMatrixMultiple(int rows, int columns) {

        myNumRows = rows;
        myNumColumns = columns;

        myNumElements = (long) myNumRows * (long) myNumColumns;
        myNumRowsPerSingleDimArray = Integer.MAX_VALUE / myNumColumns;
        myNumElementsPerSingleDimArray = myNumRowsPerSingleDimArray * myNumColumns;
        int numSingleDimArrays = (int) (myNumElements / (long) myNumElementsPerSingleDimArray);
        int numRemaining = (int) (myNumElements % (long) myNumElementsPerSingleDimArray);
        if (numRemaining != 0) {
            myData = new byte[numSingleDimArrays + 1][];
            for (int i = 0; i < numSingleDimArrays; i++) {
                myData[i] = new byte[myNumElementsPerSingleDimArray];
            }
            myData[numSingleDimArrays] = new byte[numRemaining];
        } else {
            myData = new byte[numSingleDimArrays][];
            for (int i = 0; i < numSingleDimArrays; i++) {
                myData[i] = new byte[myNumElementsPerSingleDimArray];
            }
        }

    }

    @Override
    public void set(int row, int column, byte value) {
        myData[getFirstIndex(row)][getSecondIndex(row, column)] = value;
    }

    @Override
    public void setAll(byte value) {
        int numSingleDimArrays = myData.length;
        for (int i = 0; i < numSingleDimArrays; i++) {
            Arrays.fill(myData[i], value);
        }
    }

    @Override
    public byte get(int row, int column) {
        return myData[getFirstIndex(row)][getSecondIndex(row, column)];
    }

    @Override
    public byte[] getAllColumns(int row) {

        if ((row < 0) || (row >= myNumRows)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixMultiple: getAllColumns: row: " + row);
        }

        int start = getSecondIndex(row, 0);
        byte[] result = new byte[myNumColumns];
        System.arraycopy(myData[getFirstIndex(row)], start, result, 0, myNumColumns);
        return result;

    }

    @Override
    public byte[] getColumnRange(int row, int start, int end) {

        if ((row < 0) || (row >= myNumRows)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixMultiple: getColumnRange: row: " + row);
        }

        if ((start < 0) || (start >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixMultiple: getColumnRange: start: " + start);
        }

        if ((end < 0) || (end >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixMultiple: getColumnRange: end: " + end);
        }

        if (end < start) {
            throw new IllegalArgumentException("SuperByteMatrixMultiple: getColumnRange: end: " + end + " less than start: " + start);
        }

        int startIndex = getSecondIndex(row, start);
        int numElements = end - start;
        byte[] result = new byte[numElements];
        System.arraycopy(myData[getFirstIndex(row)], startIndex, result, 0, numElements);
        return result;

    }

    @Override
    public byte[] getAllRows(int column) {

        if ((column < 0) || (column >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixMultiple: getAllRows: column: " + column);
        }

        byte[] result = new byte[myNumRows];
        for (int row = 0; row < myNumRows; row++) {
            result[row] = get(row, column);
        }
        return result;

    }

    private int getFirstIndex(int row) {
        return row / myNumRowsPerSingleDimArray;
    }

    private int getSecondIndex(int row, int column) {
        return (row % myNumRowsPerSingleDimArray) * myNumColumns + column;
    }

    @Override
    public int getNumRows() {
        return myNumRows;
    }

    @Override
    public int getNumColumns() {
        return myNumColumns;
    }

    @Override
    public boolean isColumnInnerLoop() {
        return true;
    }

    @Override
    public void reorderRows(int[] newIndices) {

        if (newIndices.length != myNumRows) {
            throw new IllegalArgumentException("SuperByteMatrixMultiple: reorderRows: index array size: " + newIndices.length + " doesn't equal num rows in matrix: " + myNumRows);
        }

        int[] tempIndices = new int[newIndices.length];
        System.arraycopy(newIndices, 0, tempIndices, 0, myNumRows);

        int currentRow = 0;
        byte[] temp = new byte[myNumColumns];

        while (currentRow < myNumRows) {

            while (currentRow < myNumRows) {
                if ((tempIndices[currentRow] == currentRow) || (tempIndices[currentRow] == -1)) {
                    tempIndices[currentRow] = -1;
                } else {
                    break;
                }
                currentRow++;
            }

            if (currentRow < myNumRows) {

                System.arraycopy(myData[getFirstIndex(currentRow)], getSecondIndex(currentRow, 0), temp, 0, myNumColumns);

                int srcRow = tempIndices[currentRow];
                int destRow = currentRow;
                while (srcRow != currentRow) {
                    System.arraycopy(myData[getFirstIndex(srcRow)], getSecondIndex(srcRow, 0), myData[getFirstIndex(destRow)], getSecondIndex(destRow, 0), myNumColumns);
                    tempIndices[destRow] = -1;
                    destRow = srcRow;
                    srcRow = tempIndices[destRow];
                }

                System.arraycopy(temp, 0, myData[getFirstIndex(destRow)], getSecondIndex(destRow, 0), myNumColumns);
                tempIndices[destRow] = -1;

            }

        }

    }

    @Override
    public void reorderColumns(int[] newIndices) {

        if (newIndices.length != myNumColumns) {
            throw new IllegalArgumentException("SuperByteMatrixMultiple: reorderColumns: index array size: " + newIndices.length + " doesn't equal num columns in matrix: " + myNumColumns);
        }

        int[] tempIndices = new int[newIndices.length];
        System.arraycopy(newIndices, 0, tempIndices, 0, myNumColumns);

        int currentRow = 0;
        byte[] temp = new byte[myNumRows];

        while (currentRow < myNumColumns) {

            while (currentRow < myNumColumns) {
                if ((tempIndices[currentRow] == currentRow) || (tempIndices[currentRow] == -1)) {
                    tempIndices[currentRow] = -1;
                } else {
                    break;
                }
                currentRow++;
            }

            if (currentRow < myNumColumns) {

                for (int r = 0; r < myNumRows; r++) {
                    temp[r] = get(r, currentRow);
                }

                int srcColumn = tempIndices[currentRow];
                int destColumn = currentRow;
                while (srcColumn != currentRow) {
                    for (int r = 0; r < myNumRows; r++) {
                        set(r, destColumn, get(r, srcColumn));
                    }
                    tempIndices[destColumn] = -1;
                    destColumn = srcColumn;
                    srcColumn = tempIndices[destColumn];
                }

                for (int r = 0; r < myNumRows; r++) {
                    set(r, destColumn, temp[r]);
                }
                tempIndices[destColumn] = -1;

            }

        }

    }

    @Override
    public void setHetsTo(byte value) {
        for (byte[] current : myData) {
            for (int i = 0; i < current.length; i++) {
                if (((current[i] >>> 4) & 0xf) != (current[i] & 0xf)) {
                    current[i] = value;
                }
            }
        }
    }

    @Override
    public void arraycopy(int row, byte[] src, int startColumn) {
        System.arraycopy(src, 0, myData[getFirstIndex(row)], getSecondIndex(row, startColumn), src.length);
    }

    @Override
    public Stream<Byte> stream() {
        return StreamSupport.stream(spliterator(), true);
    }

    @Override
    public Stream<Byte> stream(int row) {
        long start = (long) row * (long) myNumColumns;
        return StreamSupport.stream(new SuperByteMatrixMultipleSpliterator<>(start, start + (long) myNumColumns), true);
    }

    public Spliterator<Byte> spliterator() {
        return new SuperByteMatrixMultipleSpliterator<>(0, myNumElements);
    }

    class SuperByteMatrixMultipleSpliterator<T extends Byte> implements Spliterator<Byte> {

        private long myCurrentIndex;
        private final long myFence;

        SuperByteMatrixMultipleSpliterator(long currentIndex, long fence) {
            myCurrentIndex = currentIndex;
            myFence = fence;
        }

        private int firstIndex(long index) {
            return (int) (index / myNumElementsPerSingleDimArray);
        }

        private int secondIndex(long index) {
            return (int) (index % myNumElementsPerSingleDimArray);
        }

        @Override
        public void forEachRemaining(Consumer<? super Byte> action) {

            int firstIndex = firstIndex(myCurrentIndex);
            int secondIndex = secondIndex(myCurrentIndex);
            int fenceFirst = firstIndex(myFence);
            int fenceSecond = secondIndex(myFence);

            if (firstIndex == fenceFirst) {
                for (int i = secondIndex; i < fenceSecond; i++) {
                    action.accept(Byte.valueOf(myData[firstIndex][i]));
                }
            } else {
                for (int i = secondIndex; i < myNumElementsPerSingleDimArray; i++) {
                    action.accept(Byte.valueOf(myData[firstIndex][i]));
                }
                for (int i = firstIndex + 1; i < fenceFirst; i++) {
                    for (int j = 0; j < myNumElementsPerSingleDimArray; j++) {
                        action.accept(Byte.valueOf(myData[i][j]));
                    }
                }
                for (int i = 0; i < fenceSecond; i++) {
                    action.accept(Byte.valueOf(myData[fenceFirst][i]));
                }
            }

            myCurrentIndex = myFence;

        }

        @Override
        public boolean tryAdvance(Consumer<? super Byte> action) {
            if (myCurrentIndex < myFence) {
                action.accept(Byte.valueOf(myData[firstIndex(myCurrentIndex)][secondIndex(myCurrentIndex)]));
                myCurrentIndex++;
                return true;
            } else {
                return false;
            }
        }

        @Override
        public Spliterator<Byte> trySplit() {
            long lo = myCurrentIndex;
            long mid = (lo + myFence) >>> 1;
            if (lo < mid) {
                myCurrentIndex = mid;
                return new SuperByteMatrixMultipleSpliterator<>(lo, mid);
            } else {
                return null;
            }
        }

        @Override
        public long estimateSize() {
            return myFence - myCurrentIndex;
        }

        @Override
        public int characteristics() {
            return ORDERED | SIZED | IMMUTABLE | SUBSIZED;
        }
    }

}
