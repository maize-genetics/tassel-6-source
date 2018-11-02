/*
 *  TableReportBuilder
 * 
 *  Created on Aug 9, 2014
 */
package net.maizegenetics.util;

import java.io.BufferedWriter;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class TableReportBuilder {

    private static final Logger myLogger = Logger.getLogger(TableReportBuilder.class);

    private static final String DELIMITER = "\t";

    private final String myTableName;
    private final Object[] myColumnNames;
    private final List<Object[]> myData;
    private final String myFilename;
    private final BufferedWriter myWriter;
    private final boolean myInMemory;
    private final int myNumColumns;

    private TableReportBuilder(String tableName, Object[] columnNames) {
        myTableName = tableName;
        myColumnNames = columnNames;
        myNumColumns = columnNames.length;
        myData = new ArrayList<>();
        myFilename = null;
        myWriter = null;
        myInMemory = true;
    }
    
    private TableReportBuilder(String tableName, int numColumns) {
        myTableName = tableName;
        myColumnNames = null;
        myNumColumns = numColumns;
        myData = new ArrayList<>();
        myFilename = null;
        myWriter = null;
        myInMemory = true;
    }

    private TableReportBuilder(String tableName, Object[] columnNames, String filename) {
        myTableName = tableName;
        myColumnNames = columnNames;
        myNumColumns = columnNames.length;
        myData = null;
        myFilename = filename;
        myWriter = Utils.getBufferedWriter(filename);
        myInMemory = false;
        try {
            writeRow(columnNames);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TableReportBuilder: init: Problem adding headers to file: " + myFilename + ": " + e.getMessage());
        }
    }

    public static TableReport readTableReport(String saveFile) {
        return TableReportUtils.readDelimitedTableReport(saveFile, DELIMITER);
    }

    public static TableReportBuilder getInstance(String tableName, Object[] columnNames) {
        return new TableReportBuilder(tableName, columnNames);
    }
    
    public static TableReportBuilder getInstance(String tableName, int numColumns) {
        return new TableReportBuilder(tableName, numColumns);
    }

    public static TableReportBuilder getInstance(String tableName, Object[] columnNames, String filename) {
        return new TableReportBuilder(tableName, columnNames, filename);
    }

    public void add(Object[] row) {

        if (myNumColumns != row.length) {
            throw new IllegalArgumentException("TableReportBuilder: add: number of row elements: " + row.length + " doesn't equal number of headers: " + myNumColumns);
        }

        if (myInMemory) {
            myData.add(row);
        } else {
            writeRow(row);
        }

    }

    public void addElements(Object... rowElements) {
        //this uses reflection to append the arrays together.
        Object[] list = new Object[myNumColumns];
        int index = 0;
        for (Object rowElement : rowElements) {
            if (rowElement.getClass().isArray()) {
                for (int i = 0; i < Array.getLength(rowElement); i++) {
                    list[index++] = Array.get(rowElement, i);
                }
            } else {
                list[index++] = rowElement;
            }
        }
        add(list);
    }

    private void writeRow(Object[] row) {
        try {
            for (int i = 0; i < row.length; i++) {
                if (i != 0) {
                    myWriter.write(DELIMITER);
                }
                myWriter.write(row[i].toString());
            }
            myWriter.write("\n");
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TableReportBuilder: writeRow: Problem adding row to file: " + myFilename + ": " + e.getMessage());
        }
    }

    public TableReport build() {
        if (myInMemory) {
            return new SimpleTableReport(myTableName, myColumnNames, myData.toArray(new Object[myData.size()][]));
        } else {
            try {
                myWriter.close();
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("TableReportBuilder: build: Problem closing file: " + myFilename + ": " + e.getMessage());
            }
            return null;
        }
    }

}
