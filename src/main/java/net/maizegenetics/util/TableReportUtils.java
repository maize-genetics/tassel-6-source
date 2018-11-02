package net.maizegenetics.util;

import com.google.common.collect.ImmutableTable;
import com.google.common.collect.Table;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

/**
 * @author Terry Casstevens
 */
public class TableReportUtils {

    private static final Logger myLogger = Logger.getLogger(TableReportUtils.class);

    /**
     * Saves Table Report to file delimited by tabs.
     *
     * @param theTableReport table report
     * @param saveFile file
     */
    public static void saveDelimitedTableReport(TableReport theTableReport, File saveFile) {
        saveDelimitedTableReport(theTableReport, "\t", saveFile);
    }

    /**
     * Saves Table Report to file delimited by specified delimiter.
     *
     * @param theTableSource table report
     * @param delimit delimiter
     * @param saveFile file
     */
    public static void saveDelimitedTableReport(TableReport theTableSource, String delimit, File saveFile) {
        try (BufferedWriter bw = Utils.getBufferedWriter(saveFile)) {
            saveDelimitedTableReport(theTableSource, delimit, bw, true);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TableReportUtils: saveDelimitedTabeReport: problem saving file: " + saveFile.getName());
        }
    }

    /**
     * Saves Table Report to file delimited by specified delimiter.
     *
     * @param theTableSource table report
     * @param delimit delimiter
     * @param bw writer
     */
    public static void saveDelimitedTableReport(TableReport theTableSource, String delimit, BufferedWriter bw, boolean includeHeader) throws IOException {

        if (bw == null) {
            throw new IllegalArgumentException("TableReportUtils: saveDelimitedTableReport: no buffered writer specified.");
        }

        if (includeHeader) {
            Object[] colNames = theTableSource.getTableColumnNames();
            for (int j = 0; j < colNames.length; j++) {
                if (j != 0) {
                    bw.write(delimit);
                }
                bw.write(colNames[j].toString());
            }
            bw.write("\n");
        }

        for (long r = 0, n = theTableSource.getRowCount(); r < n; r++) {
            Object[] theRow = theTableSource.getRow(r);
            for (int i = 0; i < theRow.length; i++) {
                if (i != 0) {
                    bw.write(delimit);
                }
                if (theRow[i] == null) {
                    // do nothing
                } else if (theRow[i] instanceof Double) {
                    bw.write(DoubleFormat.format((Double) theRow[i]));
                } else {
                    bw.write(theRow[i].toString());
                }
            }
            bw.write("\n");
        }

    }

    public static TableReport readDelimitedTableReport(String saveFile, String delimit) {

        myLogger.info("readDelimitedTableReport: Reading: " + saveFile);
        int numLines = Utils.getNumberLines(saveFile) - 1;
        myLogger.info("readDelimitedTableReport: Num Lines (Not including header): " + numLines);

        Pattern delimitPattern = Pattern.compile(delimit);
        try (BufferedReader br = Utils.getBufferedReader(saveFile)) {
            String[] columnHeaders = delimitPattern.split(br.readLine().trim());

            ForkJoinPool pool = new ForkJoinPool();
            Object[][] data = new Object[numLines][columnHeaders.length];
            int maxNumLinesPerThread = 100000;
            for (int i = 0; i < numLines; i += maxNumLinesPerThread) {
                int numLinesForThread = Math.min(maxNumLinesPerThread, numLines - i);
                String[] lines = new String[numLinesForThread];
                for (int j = 0; j < numLinesForThread; j++) {
                    lines[j] = br.readLine().trim();
                }
                pool.execute(new SplitTableReportString(data, i, lines, delimitPattern, columnHeaders.length));
            }
            pool.shutdown();
            if (!pool.awaitTermination(6000, TimeUnit.SECONDS)) {
                throw new IllegalStateException("TableReportUtils: readDelimitedTableReport: processing threads timed out.");
            }
            return new SimpleTableReport(saveFile, columnHeaders, data);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalArgumentException("Problem creating TableReport: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
        }

    }

    private static class SplitTableReportString implements Runnable {

        private final Object[][] myData;
        private int myLineNum;
        private final String[] myLines;
        private final Pattern myPattern;
        private final int myNumColumns;

        public SplitTableReportString(Object[][] data, int lineNum, String[] lines, Pattern pattern, int numColumns) {
            myData = data;
            myLineNum = lineNum;
            myLines = lines;
            myPattern = pattern;
            myNumColumns = numColumns;
        }

        @Override
        public void run() {
            for (int i = 0, n = myLines.length; i < n; i++) {
                String[] tokens = myPattern.split(myLines[i]);
                if (tokens.length != myNumColumns) {
                    throw new IllegalStateException("TableReportUtils: SplitTableReportString: Number of values don't equal number of columns in line: " + myLineNum);
                }
                for (int c = 0; c < myNumColumns; c++) {
                    try {
                        try {
                            myData[myLineNum][c] = Integer.valueOf(tokens[c]);
                        } catch (Exception ex) {
                            myData[myLineNum][c] = Double.valueOf(tokens[c]);
                        }
                    } catch (Exception e) {
                        myData[myLineNum][c] = tokens[c];
                    }
                }
                myLineNum++;
            }
        }
    }

    public static Table<Integer, String, Object> convertTableReportToGuavaTable(TableReport tr) {
        ImmutableTable.Builder<Integer, String, Object> result = new ImmutableTable.Builder<>();
        String[] colNames = new String[tr.getColumnCount()];
        for (int i = 0; i < colNames.length; i++) {
            colNames[i] = tr.getTableColumnNames()[i].toString();
        }
        for (int i = 0; i < tr.getRowCount(); i++) {
            for (int j = 0; j < tr.getColumnCount(); j++) {
                result.put(i, colNames[j], tr.getValueAt(i, j));
            }
        }
        return result.build();
    }
}
