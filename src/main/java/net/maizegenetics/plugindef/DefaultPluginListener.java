/*
 *  DefaultPluginListener
 * 
 *  Created on Sep 23, 2015
 */
package net.maizegenetics.plugindef;

import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class DefaultPluginListener implements PluginListener {

    private static final Logger myLogger = Logger.getLogger(DefaultPluginListener.class);

    private static DefaultPluginListener SINGLETON = null;
    private final Map<Plugin, Integer> myProgressValues = new HashMap<>();

    private DefaultPluginListener() {
    }

    public static DefaultPluginListener getInstance() {
        if (SINGLETON == null) {
            SINGLETON = new DefaultPluginListener();
        }
        return SINGLETON;
    }

    @Override
    public void dataSetReturned(PluginEvent event) {
        // do nothing
    }

    @Override
    public void progress(PluginEvent event) {
        DataSet ds = (DataSet) event.getSource();
        if (ds != null) {
            List<Datum> percentage = ds.getDataOfType(Integer.class);
            Plugin plugin = ds.getCreator();
            Integer lastValue = myProgressValues.get(plugin);
            if (lastValue == null) {
                lastValue = 0;
            }

            if (percentage.size() > 0) {
                Datum datum = percentage.get(0);
                Integer percent = (Integer) datum.getData();
                if (percent >= lastValue) {
                    LocalDateTime time = LocalDateTime.now();
                    String timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
                    myLogger.info(ds.getCreator().getClass().getName() + ": time: " + timeStr + ": progress: " + percent + "%");
                    lastValue += 10;
                    myProgressValues.put(plugin, lastValue);
                }
                if (percent >= 100) {
                    myProgressValues.remove(plugin);
                }

            }
        }
    }

}
