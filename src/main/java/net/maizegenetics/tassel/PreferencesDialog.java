package net.maizegenetics.tassel;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.ParameterCache;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.prefs.TasselPrefs;

/**
 * @author Terry Casstevens
 */
public class PreferencesDialog extends AbstractPlugin {

    private PluginParameter<Boolean> myRetainRareAlleles = new PluginParameter.Builder<>("retainRareAlleles", TasselPrefs.ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT, Boolean.class)
            .description("True if rare alleles should be retained.  This has no effect on Nucleotide Data as all alleles will be retained regardless.")
            .build();

    private PluginParameter<Boolean> mySendLogToConsole = new PluginParameter.Builder<>("sendLogToConsole", TasselPrefs.TASSEL_LOG_SEND_TO_CONSOLE_DEFAULT, Boolean.class)
            .description("Flag whether to send logging to the console.")
            .build();

    private PluginParameter<String> myConfigFile = new PluginParameter.Builder<>("configFile", TasselPrefs.TASSEL_CONFIG_FILE_DEFAULT, String.class)
            .description("Global configuration file")
            .required(false)
            .inFile()
            .build();

    public PreferencesDialog(boolean isInteractive) {
        super(isInteractive);
        // Load global parameter / value is config file specified
        ParameterCache.load(TasselPrefs.getConfigFile());
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        setParameter(myRetainRareAlleles, TasselPrefs.getAlignmentRetainRareAlleles());
        setParameter(mySendLogToConsole, TasselPrefs.getLogSendToConsole());
        setParameter(myConfigFile, TasselPrefs.getConfigFile());
    }

    @Override
    public DataSet processData(DataSet input) {

        TasselPrefs.putAlignmentRetainRareAlleles(retainRareAlleles());

        TasselPrefs.putLogSendToConsole(sendLogToConsole());
        TasselLogging.updateLoggingLocation();

        TasselPrefs.putConfigFile(configFile());
        ParameterCache.load(TasselPrefs.getConfigFile());
        TASSELGUI.instance.updatePluginsWithGlobalConfigParameters();

        return null;

    }

    /**
     * Retain Rare Alleles
     *
     * @return Retain Rare Alleles
     */
    public Boolean retainRareAlleles() {
        return myRetainRareAlleles.value();
    }

    /**
     * Set Retain Rare Alleles. Retain Rare Alleles
     *
     * @param value Retain Rare Alleles
     *
     * @return this plugin
     */
    public PreferencesDialog retainRareAlleles(Boolean value) {
        myRetainRareAlleles = new PluginParameter<>(myRetainRareAlleles, value);
        return this;
    }

    /**
     * Flag whether to send logging to the console.
     *
     * @return Send Log To Console
     */
    public Boolean sendLogToConsole() {
        return mySendLogToConsole.value();
    }

    /**
     * Set Send Log To Console. Flag whether to send logging to the console.
     *
     * @param value Send Log To Console
     *
     * @return this plugin
     */
    public PreferencesDialog sendLogToConsole(Boolean value) {
        mySendLogToConsole = new PluginParameter<>(mySendLogToConsole, value);
        return this;
    }

    /**
     * Global configuration file
     *
     * @return Config File
     */
    public String configFile() {
        return myConfigFile.value();
    }

    /**
     * Set Config File. Global configuration file
     *
     * @param value Config File
     *
     * @return this plugin
     */
    public PreferencesDialog configFile(String value) {
        myConfigFile = new PluginParameter<>(myConfigFile, value);
        return this;
    }

    @Override
    public String icon() {
        return "/images/preferences.gif";
    }

    @Override
    public String getButtonName() {
        return "Preferences";
    }

    @Override
    public String getToolTipText() {
        return "Preferences";
    }
}
