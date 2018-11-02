package net.maizegenetics.util;

import java.io.Serializable;

/**
 * Created by IntelliJ IDEA. User: ed Date: Sep 28, 2006 Time: 9:37:46 PM
 */
public class SimpleTableReport extends AbstractTableReport implements Serializable, TableReport {

    private final Object[][] myData;
    private final Object[] myColumnNames;
    private final String myName;

    public SimpleTableReport(String theName, Object[] columnNames, Object[][] theData) {
        myData = theData;
        myColumnNames = columnNames;
        myName = theName;
    }

    public SimpleTableReport(TableReport tr) {
        int numRows = (int) tr.getRowCount();
        if ((long) numRows != tr.getRowCount()) {
            throw new IllegalArgumentException("SimpleTableReport: init: This implementation can't support more rows than: " + Integer.MAX_VALUE);
        }
        myData = new Object[numRows][tr.getColumnCount()];
        for (int i = 0; i < numRows; i++) {
            System.arraycopy(tr.getRow(i), 0, myData[i], 0, numRows);
        }
        myColumnNames = tr.getTableColumnNames();
        myName = tr.getTableTitle();
    }

    /**
     * Return column names for the table
     */
    @Override
    public Object[] getTableColumnNames() {
        return myColumnNames;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    @Override
    public Object[] getRow(long row) {
        return myData[(int) row];
    }

    /**
     * Return the name for the title of the ANOVA
     */
    @Override
    public String getTableTitle() {
        return myName;
    }

    @Override
    public long getRowCount() {
        return myData.length;
    }

    @Override
    public long getElementCount() {
        return getRowCount() * getColumnCount();
    }

    @Override
    public int getColumnCount() {
        return myColumnNames.length;
    }

    @Override
    public Object getValueAt(long row, int col) {
        return myData[(int) row][col];
    }
}
