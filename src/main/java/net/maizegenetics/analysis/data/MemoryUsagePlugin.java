/*
 *  MemoryUsagePlugin
 *
 *  Created on April 18, 2016
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.Sizeof;

/**
 * @author Terry Casstevens
 */
public class MemoryUsagePlugin extends AbstractPlugin {

    public MemoryUsagePlugin(boolean isInteractive) {
        super(isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        Sizeof.printMemoryUse();
        return input;
    }

    @Override
    public String getButtonName() {
        return "Memory Usage";
    }

    @Override
    public String getToolTipText() {
        return "Memory Usage";
    }

}
