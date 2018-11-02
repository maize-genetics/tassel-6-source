/*
 * PluginEvent.java
 *
 */
package net.maizegenetics.plugindef;

import java.util.EventObject;

/**
 *
 * @author Terry Casstevens
 */
public class PluginEvent extends EventObject {

    private final Object myMetaData;

    /**
     * Creates a new instance of PluginEvent
     */
    public PluginEvent(DataSet source) {
        this(source, null);
    }

    /**
     * Creates a new instance of PluginEvent
     */
    public PluginEvent(DataSet source, Object metaData) {
        super(source == null ? new DataSet((Datum) null, null) : source);
        myMetaData = metaData;
    }

    public Object getMetaData() {
        return myMetaData;
    }

}
