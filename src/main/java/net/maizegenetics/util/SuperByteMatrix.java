/*
 *  SuperByteMatrix
 */
package net.maizegenetics.util;

import java.util.stream.Stream;

/**
 *
 * @author Terry Casstevens
 */
public interface SuperByteMatrix {

    /**
     * Return number of rows.
     *
     * @return number of rows
     */
    public int getNumRows();

    /**
     * Return number of columns.
     *
     * @return number of columns
     */
    public int getNumColumns();

    /**
     * Sets value at given row and column.
     *
     * @param row row
     * @param column column
     * @param value value
     */
    public void set(int row, int column, byte value);

    /**
     * Sets values at given row and starting column.
     *
     * @param row row
     * @param src values
     * @param startColumn start column
     */
    public void arraycopy(int row, byte[] src, int startColumn);

    /**
     * Sets value for all elements.
     *
     * @param value value
     */
    public void setAll(byte value);

    /**
     * Gets value at given row and column.
     *
     * @param row row
     * @param column column
     *
     * @return value
     */
    public byte get(int row, int column);

    /**
     * Get all values for given row.
     *
     * @param row row
     *
     * @return values
     */
    public byte[] getAllColumns(int row);

    /**
     * Get values for given row from start column (inclusive) to end column
     * (exclusive).
     *
     * @param row row
     * @param start start
     * @param end end
     *
     * @return values
     */
    public byte[] getColumnRange(int row, int start, int end);

    /**
     * Get all values for give column.
     *
     * @param column column
     *
     * @return values
     */
    public byte[] getAllRows(int column);

    /**
     * Returns true if the matrix stored for better performance when column loop
     * inside row loop. False if matrix stored for better performance when row
     * loop inside column loop.
     *
     * @return true if the matrix stored for better performance when column loop
     * inside row loop. False if matrix stored for better performance when row
     * loop inside column loop.
     */
    public boolean isColumnInnerLoop();

    /**
     * Reorders rows of this matrix based on the given indices.
     *
     * @param newIndices new indices.
     */
    public void reorderRows(int[] newIndices);

    /**
     * Reorders columns of this matrix based on the given indices.
     *
     * @param newIndices new indices.
     */
    public void reorderColumns(int[] newIndices);

    /**
     * Changes all heterozygous values to give value.
     *
     * @param value value
     */
    public void setHetsTo(byte value);

    /**
     * Returns a Stream over the bytes of this matrix.
     *
     * @return Stream over the bytes of this matrix
     */
    public Stream<Byte> stream();
    
    public Stream<Byte> stream(int row);
}
