/*
 * AbstractTableReport
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public abstract class AbstractTableReport implements TableReport {

    private long currentRowNumber = -1;
    private Object[] currentRow = null;

    @Override
    public Object getValueAt(long row, int col) {
        if (row != currentRowNumber) {
            currentRowNumber = row;
            currentRow = getRow(row);
        }
        return currentRow[col];
    }

}
