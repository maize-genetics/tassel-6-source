/*
 * Plugin.java
 *
 */
package net.maizegenetics.plugindef;

import net.maizegenetics.util.ProgressListener;
import org.apache.log4j.Logger;

import java.lang.reflect.Constructor;

/**
 * @author Terry Casstevens
 */
public interface Plugin extends PluginListener, ProgressListener, Runnable {

    /**
     * Performs function of this plugin.
     *
     * @param input input
     *
     * @return resulting data set or null.
     */
    public DataSet performFunction(DataSet input);

    /**
     * For the new Generic Plugin Parameter design, performFunction() will automatically call this. Therefore, coders of
     * Plugins should override this instead of performFunction().
     *
     * @param input input
     *
     * @return resulting data set or null.
     */
    public DataSet processData(DataSet input);

    /**
     * Returns parameter value for given parameter key.
     *
     * @param key key
     *
     * @return value
     */
    public Object getParameter(Enum key);

    /**
     * Returns parameter value for given parameter key.
     *
     * @param key key
     *
     * @return value
     */
    public Object getParameter(String key);

    /**
     * Sets parameter value for a given plugin
     *
     * @param param parameter
     * @param value value
     *
     * @return this plugin
     */
    public Plugin setParameter(PluginParameter<?> param, Object value);

    /**
     * Sets parameter value for given parameter key.
     *
     * @param key key
     * @param value value
     *
     * @return this plugin
     */
    public Plugin setParameter(String key, Object value);

    /**
     * Sets parameter value for given parameter key.
     *
     * @param key key
     * @param value value
     *
     * @return this plugin
     */
    public Plugin setParameter(String key, String value);

    /**
     * Sets all parameter values to default.
     */
    public void setParametersToDefault();

    /**
     * Sets up this plugin to receive input from another plugin.
     *
     * @param input input
     */
    public void receiveInput(Plugin input);

    /**
     * If interactive = true, the plugin will create dialogs and panels to interacts with the user
     *
     * @return boolean
     */
    public boolean isInteractive();

    /**
     * Icon filename for this plugin to be used in buttons, etc.
     *
     * @return Icon filename
     */
    public String icon();

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName();

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText();

    /**
     * Adds listener to this plugin.
     *
     * @param listener listener to add
     */
    public void addListener(PluginListener listener);

    /**
     * Set whether this plugin is threaded.
     *
     * @param threaded whether to be threaded.
     */
    public void setThreaded(boolean threaded);

    /**
     * Attempt to cancel processing.
     *
     * @return true if plugin will cancel itself.
     */
    public boolean cancel();

    /**
     * Allows self-describing Plugins to use args to set parameters specific to itself.
     *
     * @param args arguments
     */
    public void setParameters(String[] args);

    /**
     * Returns Citation for this plugin.
     *
     * @return Citation
     */
    public String getCitation();

    /**
     * Returns description of the plugin.
     *
     * @return description
     */
    public String pluginDescription();

    /**
     * Returns URL to User Manual.
     *
     * @return URL
     */
    public String pluginUserManualURL();

    /**
     * Gets the Usage Statement for this Plugin.
     *
     * @return Usage Statement
     */
    public String getUsage();

    public boolean wasCancelled();

    Logger myLogger = Logger.getLogger(Plugin.class);

    public static Plugin getPluginInstance(String className) {
        return getPluginInstance(className, false);
    }

    /**
     * Gets instance of Plugin
     *
     * @param className class name
     * @param isInteractive is interactive
     *
     * @return Plugin
     */
    public static Plugin getPluginInstance(String className, boolean isInteractive) {
        try {
            Class currentMatch = Class.forName(className);
            Constructor constructor = currentMatch.getConstructor(boolean.class);
            return (Plugin) constructor.newInstance(isInteractive);
        } catch (Exception ex) {
            try {
                Class currentMatch = Class.forName(className);
                Constructor constructor = currentMatch.getConstructor();
                return (Plugin) constructor.newInstance();
            } catch (NoSuchMethodException nsme) {
                myLogger.warn("Self-describing Plugins should implement this constructor: " + className);
                myLogger.warn("public Plugin(Frame parentFrame, boolean isInteractive) {");
                myLogger.warn("   super(parentFrame, isInteractive);");
                myLogger.warn("}");
                return null;
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                return null;
            }
        }
    }

    /**
     * Returns whether given class name is Plugin.
     *
     * @param className class name
     *
     * @return true if class is Plugin
     */
    public static boolean isPlugin(String className) {
        try {
            Class currentMatch = Class.forName(className);
            return !currentMatch.isInterface() && Plugin.class.isAssignableFrom(currentMatch);
        } catch (Exception | Error ex) {
            return false;
        }
    }
}
