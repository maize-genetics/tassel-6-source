/*
 * CombineDataSetsPlugin.java
 *
 * Created on January 5, 2007, 2:25 AM
 *
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginEvent;

import java.util.*;

/**
 * @author Terry Casstevens
 */
public class CombineDataSetsPlugin extends AbstractPlugin {

    private final Map<Plugin, DataSet> myDataSets = new LinkedHashMap<>();
    private final Map<Plugin, DataSet> myOnceDataSets = new LinkedHashMap<>();

    /**
     * Creates a new instance of CombineDataSetsPlugin
     */
    public CombineDataSetsPlugin() {
        super(false);
    }

    /**
     * Returns combined data set if all inputs have been received.
     *
     * @param dataSet Not used. All input received from listening to other
     * plugins.
     */
    @Override
    public DataSet performFunction(DataSet dataSet) {

        try {

            List<DataSet> dataSets;
            synchronized (myDataSets) {
                if ((myDataSets.containsValue(null)) || myOnceDataSets.containsValue(null)) {
                    return null;
                }

                dataSets = new ArrayList<>();

                dataSets.addAll(myDataSets.values());
                dataSets.addAll(myOnceDataSets.values());

                reset();
            }

            DataSet result = DataSet.getDataSet(dataSets, this);
            fireDataSetReturned(result);

            return result;

        } finally {
            fireProgress(100);
        }

    }

    public void reset() {

        // Clear only values.
        // Method dataSetReturned knows what inputs
        // to expect based on keys stored here.
        Set<Plugin> keys = myDataSets.keySet();
        for (Iterator<Plugin> itr = keys.iterator(); itr.hasNext(); ) {
            myDataSets.put(itr.next(), null);
        }

    }

    @Override
    public String getToolTipText() {
        return "Combine Datasets";
    }

    @Override
    public String getButtonName() {
        return "Combine";
    }

    @Override
    public void dataSetReturned(PluginEvent event) {

        DataSet input = (DataSet) event.getSource();
        Plugin creator = input.getCreator();

        if (myOnceDataSets.containsKey(creator)) {
            Object value = myOnceDataSets.get(creator);
            if (value != null) {
                throw new IllegalStateException("CombineDataSetsPlugin: dataSetReturned: this plugin should only return data once: " + creator);
            } else {
                myOnceDataSets.put(creator, input);
            }
        } else if (myDataSets.containsKey(creator)) {
            Object value = myDataSets.get(creator);
            if (value != null) {
                throw new IllegalStateException("CombineDataSetsPlugin: dataSetReturned: this plugin should only return data once per iteration: " + creator);
            } else {
                myDataSets.put(creator, input);
            }
        } else {
            throw new IllegalStateException("CombineDataSetsPlugin: dataSetReturned: can not receive data from unknown plugin: " + creator);
        }

        performFunction(null);

    }

    /**
     * Add given plugin as source to receive data sets only once and use that
     * data set in every resulting output.
     *
     * @param plugin plugin
     */
    public void receiveDataSetOnceFrom(Plugin plugin) {
        super.receiveInput(plugin);
        myOnceDataSets.put(plugin, null);
    }

    /**
     * Add given plugin as source to receive data sets iteratively.
     *
     * @param plugin plugin
     */
    public void receiveDataSetFrom(Plugin plugin) {
        super.receiveInput(plugin);
        myDataSets.put(plugin, null);
    }

    /**
     * Sets up this plugin to receive input from another plugin.
     *
     * @param input input
     */
    @Override
    public void receiveInput(Plugin input) {
        receiveDataSetFrom(input);
    }

    @Override
    public String toString() {

        StringBuilder str = new StringBuilder();

        Iterator<DataSet> itr = myDataSets.values().iterator();
        while (itr.hasNext()) {
            DataSet current = itr.next();
            if (current != null) {
                str.append(current.toString());
            }
        }

        return str.toString();

    }
}
