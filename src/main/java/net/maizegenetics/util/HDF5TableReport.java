/*
 *  HDF5TableReport
 *
 *  Created on Jun 27, 2014
 */
package net.maizegenetics.util;

import ch.systemsx.cisd.hdf5.HDF5DataSetInformation;
import ch.systemsx.cisd.hdf5.HDF5DataTypeInformation;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Terry Casstevens
 */
public class HDF5TableReport implements TableReport {

    private final IHDF5Reader myReader;
    private final Map<Integer, List<String>> myData = new HashMap<>();
    private final String[] myColumnHeadings;

    public HDF5TableReport(String filename) {
        myReader = HDF5Factory.openForReading(filename);
        traverseNode("/", 0);
        myColumnHeadings = new String[myData.size()];
        for (int i = 0; i < myData.size(); i++) {
            myColumnHeadings[i] = String.valueOf(i);
        }
    }

    @Override
    public Object[] getTableColumnNames() {
        return myColumnHeadings;
    }

    @Override
    public String getTableTitle() {
        return "HDF5 Schema";
    }

    @Override
    public int getColumnCount() {
        return myData.size();
    }

    @Override
    public long getRowCount() {
        return myData.get(0).size();
    }

    @Override
    public long getElementCount() {
        return getColumnCount() * getRowCount();
    }

    @Override
    public Object[] getRow(long row) {
        Object[] result = new Object[getColumnCount()];
        for (int i = 0; i < getColumnCount(); i++) {
            result[i] = getValueAt(row, i);
        }
        return result;
    }

    @Override
    public Object getValueAt(long row, int col) {
        return myData.get(col).get((int) row);
    }

    private void traverseNode(String node, int column) {
        addAttributes(node, column);
        List<String> members = myReader.object().getAllGroupMembers(node);
        for (String str : members) {
            StringBuilder builder = new StringBuilder();
            builder.append(str);
            String path = node + "/" + str;

            if (myReader.isGroup(path)) {
                setValue(column, builder.toString());
                traverseNode(path, column + 1);
            } else if (myReader.object().isDataSet(path)) {
                HDF5DataSetInformation info = myReader.getDataSetInformation(path);
                builder.append(" : ").append(myReader.getDataSetInformation(path).getTypeInformation().tryGetJavaType().getName());
                long[] dimensions = info.getDimensions();
                if (dimensions.length != 0) {
                    builder.append(" [");
                    boolean first = true;
                    for (long dimension : dimensions) {
                        if (!first) {
                            builder.append(", ");
                        } else {
                            first = false;
                        }
                        builder.append(dimension);
                    }
                    builder.append("]");
                }
                setValue(column, builder.toString());
                addAttributes(path, column + 1);
            } else {
                setValue(column, builder.toString());
            }
        }
    }

    private void addAttributes(String path, int column) {
        List<String> attributes = myReader.object().getAttributeNames(path);
        for (String attribute : attributes) {
            StringBuilder builder = new StringBuilder();
            HDF5DataTypeInformation info = myReader.object().getAttributeInformation(path, attribute);
            builder.append(attribute);

            Class<?> clazz = info.tryGetJavaType();
            if (clazz.isAssignableFrom(int.class)) {
                int value = myReader.int32().getAttr(path, attribute);
                builder.append(" : ").append(value);
            } else if (clazz.isAssignableFrom(boolean.class)) {
                boolean value = myReader.bool().getAttr(path, attribute);
                builder.append(" : ").append(value);
            } else if (clazz.isAssignableFrom(String.class)) {
                String value = myReader.string().getAttr(path, attribute);
                builder.append(" : ").append(value);
            }

            builder.append(" (").append(clazz.getName()).append(" attribute)");
            setValue(column, builder.toString());
        }
    }

    private void setValue(int column, String value) {

        List<String> currentColumn = myData.get(column);
        if (currentColumn == null) {
            currentColumn = new ArrayList<>();
            if (!myData.isEmpty()) {
                currentColumn.add("");
                for (int i = 1; i < myData.get(0).size(); i++) {
                    currentColumn.add(null);
                }
            }
        }

        myData.put(column, currentColumn);

        for (List<String> currentList : myData.values()) {
            if (currentList == currentColumn) {
                currentColumn.add(value);
            } else {
                currentList.add(null);
            }
        }

    }

}
