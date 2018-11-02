/*
 *  SuperByteMatrixSingleTranspose
 */
package net.maizegenetics.util;

import java.util.stream.Stream;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrixTranspose implements SuperByteMatrix {

    private final SuperByteMatrix myMatrix;

    SuperByteMatrixTranspose(int rows, int columns) {
        myMatrix = SuperByteMatrixBuilder.getInstance(columns, rows);
    }

    @Override
    public int getNumRows() {
        return myMatrix.getNumColumns();
    }

    @Override
    public int getNumColumns() {
        return myMatrix.getNumRows();
    }

    @Override
    public void set(int row, int column, byte value) {
        myMatrix.set(column, row, value);
    }

    @Override
    public void setAll(byte value) {
        myMatrix.setAll(value);
    }

    @Override
    public byte get(int row, int column) {
        return myMatrix.get(column, row);
    }

    @Override
    public byte[] getAllColumns(int row) {
        return myMatrix.getAllRows(row);
    }

    @Override
    public byte[] getColumnRange(int row, int start, int end) {
        int length = end - start;
        byte[] result = new byte[length];
        for (int i = 0; i < length; i++) {
            result[i] = get(row, i);
        }
        return result;
    }

    @Override
    public byte[] getAllRows(int column) {
        return myMatrix.getAllColumns(column);
    }

    @Override
    public boolean isColumnInnerLoop() {
        return false;
    }

    @Override
    public void reorderRows(int[] newIndices) {
        myMatrix.reorderColumns(newIndices);
    }

    @Override
    public void reorderColumns(int[] newIndices) {
        myMatrix.reorderRows(newIndices);
    }

    @Override
    public void setHetsTo(byte value) {
        myMatrix.setHetsTo(value);
    }

    @Override
    public void arraycopy(int row, byte[] src, int startColumn) {
        for (int i = 0; i < src.length; i++) {
            set(row, startColumn + i, src[i]);
        }
    }

    @Override
    public Stream<Byte> stream() {
        return myMatrix.stream();
    }

    @Override
    public Stream<Byte> stream(int row) {
        return myMatrix.stream(row);
    }
}
