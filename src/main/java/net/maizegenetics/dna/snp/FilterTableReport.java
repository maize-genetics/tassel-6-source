/*
 *  FilterTableReport
 * 
 *  Created on Jul 1, 2014
 */
package net.maizegenetics.dna.snp;

import java.util.Map;
import net.maizegenetics.dna.snp.FilterSite.FILTER_SITES_ATTRIBUTES;
import net.maizegenetics.dna.snp.FilterTaxa.FILTER_TAXA_ATTRIBUTES;
import net.maizegenetics.util.TableReport;

/**
 *
 * @author Terry Casstevens
 */
public class FilterTableReport implements TableReport {

    private final static String[] COLUMN_HEADERS = new String[]{"Key", "Value"};

    private final int myNumRows;
    private final String[] myRowHeaders;
    private final Object[] myFilterAttributes;

    public FilterTableReport(FilterList filters) {

        int numRows = 0;
        for (Filter filter : filters) {
            numRows += filter.numAttributes();
            numRows++;
        }
        myNumRows = numRows;

        myFilterAttributes = new Object[myNumRows];
        myRowHeaders = new String[myNumRows];

        int count = 0;
        for (Filter filter : filters) {
            if (filter instanceof FilterSite) {
                for (Map.Entry<FILTER_SITES_ATTRIBUTES, Object> current : ((FilterSite) filter).attributes().entrySet()) {
                    myRowHeaders[count] = current.getKey().name();
                    myFilterAttributes[count] = current.getValue();
                    count++;
                }
            } else if (filter instanceof FilterTaxa) {
                for (Map.Entry<FILTER_TAXA_ATTRIBUTES, Object> current : ((FilterTaxa) filter).attributes().entrySet()) {
                    myRowHeaders[count] = current.getKey().name();
                    myFilterAttributes[count] = current.getValue();
                    count++;
                }
            }

            count++;
        }

    }

    @Override
    public Object[] getTableColumnNames() {
        return COLUMN_HEADERS;
    }

    @Override
    public String getTableTitle() {
        return "Filter";
    }

    @Override
    public int getColumnCount() {
        return COLUMN_HEADERS.length;
    }

    @Override
    public long getRowCount() {
        return myNumRows;
    }

    @Override
    public long getElementCount() {
        return getColumnCount() * getRowCount();
    }

    @Override
    public Object[] getRow(long row) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object getValueAt(long row, int col) {
        switch (col) {
            case 0:
                return myRowHeaders[(int) row];
            case 1:
                return myFilterAttributes[(int) row];
            default:
                throw new IllegalArgumentException("FilterTableReport: getValueAt: unknown column: " + col);
        }
    }

}
