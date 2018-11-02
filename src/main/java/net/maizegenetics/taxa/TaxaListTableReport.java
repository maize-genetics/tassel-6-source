/*
 * TaxaListTableReport
 */
package net.maizegenetics.taxa;

import com.google.common.base.Joiner;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.util.TableReport;

/**
 *
 * @author Terry Casstevens
 */
public class TaxaListTableReport implements TableReport {

    private static final String[] DEFAULT_COLUMN_HEADINGS = new String[]{"Taxa"};

    private final TaxaList myTaxaList;
    private final String[] myColumnHeadings;

    public TaxaListTableReport(TaxaList taxaList) {
        myTaxaList = taxaList;
        List<String> annotationColumns = new ArrayList<>();
        for (Taxon current : myTaxaList) {
            for (String key : current.getAnnotation().getAnnotationKeys()) {
                if (!annotationColumns.contains(key)) {
                    annotationColumns.add(key);
                }
            }
        }
        int totalHeadings = DEFAULT_COLUMN_HEADINGS.length + annotationColumns.size();
        myColumnHeadings = new String[totalHeadings];
        for (int i = 0; i < DEFAULT_COLUMN_HEADINGS.length; i++) {
            myColumnHeadings[i] = DEFAULT_COLUMN_HEADINGS[i];
        }
        for (int j = DEFAULT_COLUMN_HEADINGS.length; j < totalHeadings; j++) {
            myColumnHeadings[j] = annotationColumns.get(j - DEFAULT_COLUMN_HEADINGS.length);
        }
    }

    @Override
    public Object[] getTableColumnNames() {
        return myColumnHeadings;
    }

    @Override
    public String getTableTitle() {
        return "Taxa List";
    }

    @Override
    public int getColumnCount() {
        return myColumnHeadings.length;
    }

    @Override
    public long getRowCount() {
        return myTaxaList.numberOfTaxa();
    }

    @Override
    public long getElementCount() {
        return getColumnCount() * getRowCount();
    }

    @Override
    public Object[] getRow(long row) {
        int numColumns = getColumnCount();
        Object[] result = new Object[numColumns];
        for (int c = 0; c < numColumns; c++) {
            result[c] = getValueAt(row, c);
        }
        return result;
    }

    @Override
    public Object getValueAt(long row, int col) {
        switch (col) {
            case 0:
                return myTaxaList.get((int) row).getName();
            default:
                String[] annotations = myTaxaList.get((int) row).getAnnotation().getTextAnnotation(myColumnHeadings[col]);
                if ((annotations != null) && (annotations.length != 0)) {
                    return Joiner.on(",").join(annotations);
                } else {
                    return null;
                }
        }
    }

    public TaxaList getTaxaList() {
        return myTaxaList;
    }

}
