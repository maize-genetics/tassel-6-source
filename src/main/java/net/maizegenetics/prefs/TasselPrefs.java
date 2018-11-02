/*
 * TasselPrefs.java
 *
 * Created on August 5, 2007, 6:58 PM
 *
 */
package net.maizegenetics.prefs;

import java.util.HashMap;
import java.util.Map;
import java.util.prefs.Preferences;

/**
 * @author Terry Casstevens
 */
public class TasselPrefs {

    private static boolean PERSIST_PREFERENCES = false;
    private static final Map<String, Object> TEMP_CACHED_VALUES = new HashMap<>();
    //
    // Top level preferences
    //
    public static final String TASSEL_TOP = "/tassel";
    public static final String TASSEL_SAVE_DIR = "saveDir";
    public static final String TASSEL_SAVE_DIR_DEFAULT = "";
    public static final String TASSEL_OPEN_DIR = "openDir";
    public static final String TASSEL_OPEN_DIR_DEFAULT = "";
    public static final String TASSEL_X_DIM = "xDimension";
    public static final int TASSEL_X_DIM_DEFAULT = -1;
    public static final String TASSEL_Y_DIM = "yDimension";
    public static final int TASSEL_Y_DIM_DEFAULT = -1;
    public static final String TASSEL_LOG_SEND_TO_CONSOLE = "logToConsole";
    public static final boolean TASSEL_LOG_SEND_TO_CONSOLE_DEFAULT = false;
    public static final String TASSEL_LOG_DEBUG = "logDebug";
    public static final boolean TASSEL_LOG_DEBUG_DEFAULT = false;
    public static final String TASSEL_LOG_X_DIM = "logxDimension";
    public static final int TASSEL_LOG_X_DIM_DEFAULT = -1;
    public static final String TASSEL_LOG_Y_DIM = "logyDimension";
    public static final int TASSEL_LOG_Y_DIM_DEFAULT = -1;
    public static final String TASSEL_MAX_THREADS = "maxThreads";
    public static final int TASSEL_MAX_THREADS_DEFAULT = Math.max(Runtime.getRuntime().availableProcessors() - 1, 1);
    public static final String TASSEL_CONFIG_FILE = "configFile";
    public static final String TASSEL_CONFIG_FILE_DEFAULT = "";
    //
    // ExportPlugin preferences
    //
    public static final String EXPORT_PLUGIN_TOP = "/tassel/plugins/export";
    // Export as Diploids
    public static final String EXPORT_PLUGIN_EXPORT_DIPLOIDS = "exportDiploids";
    public static final boolean EXPORT_PLUGIN_EXPORT_DIPLOIDS_DEFAULT = false;
    // Include Taxa Annotations
    public static final String EXPORT_PLUGIN_INCLUDE_TAXA_ANNOTATIONS = "includeTaxaAnnotations";
    public static final boolean EXPORT_PLUGIN_INCLUDE_TAXA_ANNOTATIONS_DEFAULT = true;
    //
    // FilterAlignmentPlugin preferences
    //
    public static final String FILTER_ALIGN_PLUGIN_TOP = "/tassel/plugins/filterAlign";
    // Min. frequency for filtering sites.
    public static final String FILTER_ALIGN_PLUGIN_MIN_FREQ = "minFreq";
    public static final double FILTER_ALIGN_PLUGIN_MIN_FREQ_DEFAULT = 0.0;
    // Max. frequency for filtering sites.
    public static final String FILTER_ALIGN_PLUGIN_MAX_FREQ = "maxFreq";
    public static final double FILTER_ALIGN_PLUGIN_MAX_FREQ_DEFAULT = 1.0;
    // Min. frequency for filtering sites.
    public static final String FILTER_ALIGN_PLUGIN_MIN_COUNT = "minCount";
    public static final int FILTER_ALIGN_PLUGIN_MIN_COUNT_DEFAULT = 1;
    //
    // FilterTaxaPropertiesPlugin preferences
    //
    public static final String FILTER_TAXA_PROPS_PLUGIN_TOP = "/tassel/plugins/filterTaxaAlign";
    // Min. Not Missing Gametes Proportion
    public static final String FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING = "minNotMissingFreq";
    public static final double FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING_DEFAULT = 0.0;
    //Min. Heterozygotes Proportion
    public static final String FILTER_TAXA_PROPS_PLUGIN_MIN_HET = "minHetFreq";
    public static final double FILTER_TAXA_PROPS_PLUGIN_MIN_HET_DEFAULT = 0.0;
    //Max. Heterozygotes Proportion
    public static final String FILTER_TAXA_PROPS_PLUGIN_MAX_HET = "maxHetFreq";
    public static final double FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT = 1.0;
    //
    // Alignment preferences
    //
    public static final String ALIGNMENT_TOP = "/tassel/alignment";
    // Retain Rare Alleles
    public static final String ALIGNMENT_RETAIN_RARE_ALLELES = "retainRareAlleles";
    public static final boolean ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT = false;
    //
    // GOBII preferences
    //
    public static final String GOBII_TOP = "/tassel/gobii";
    // Postgres
    public static final String GOBII_USER = "user";
    public static final String GOBII_USER_DEFAULT = "";
    public static final String GOBII_DB = "db";
    public static final String GOBII_DB_DEFAULT = "";
    // BMS
    public static final String BMS_USER = "bmsuser";
    public static final String BMS_USER_DEFAULT = "";
    public static final String BMS_HOST = "bmshost";
    public static final String BMS_HOST_DEFAULT = "localhost";
    public static final String BMS_DB = "bmsdb";
    public static final String BMS_DB_DEFAULT = "";

    /**
     * Creates a new instance of TasselPrefs
     */
    private TasselPrefs() {
    }

    public static boolean getPersistPreferences() {
        return PERSIST_PREFERENCES;
    }

    /**
     * Whether to Persist Preferences. Preference changes should be persisted when executing GUI and set only
     * temporarily from Command Line Flags. Also getting preferences should use stored values when executing GUI. And
     * should use default values (if not temporarily set) when executing from Command Line.
     *
     * @param persist whether to persist preferences
     */
    public static void setPersistPreferences(boolean persist) {
        PERSIST_PREFERENCES = persist;
    }

    public static String getPref(String path, String key, String def) {
        String pref = path + "/" + key;
        String result = (String) TEMP_CACHED_VALUES.get(pref);
        if (result != null) {
            return result;
        }
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            result = node.get(key, def);
        } else {
            result = def;
        }
        TEMP_CACHED_VALUES.put(pref, result);
        return result;
    }

    public static void putPref(String path, String key, String value) {
        String pref = path + "/" + key;
        TEMP_CACHED_VALUES.put(pref, value);
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            node.put(key, value);
        }
    }

    public static double getDoublePref(String path, String key, double def) {
        String pref = path + "/" + key;
        Double result = (Double) TEMP_CACHED_VALUES.get(pref);
        if (result != null) {
            return result;
        }
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            result = node.getDouble(key, def);
        } else {
            result = def;
        }
        TEMP_CACHED_VALUES.put(pref, result);
        return result;
    }

    public static void putDoublePref(String path, String key, double value) {
        String pref = path + "/" + key;
        TEMP_CACHED_VALUES.put(pref, value);
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            node.putDouble(key, value);
        }
    }

    public static int getIntPref(String path, String key, int def) {
        String pref = path + "/" + key;
        Integer result = (Integer) TEMP_CACHED_VALUES.get(pref);
        if (result != null) {
            return result;
        }
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            result = node.getInt(key, def);
        } else {
            result = def;
        }
        TEMP_CACHED_VALUES.put(pref, result);
        return result;
    }

    public static void putIntPref(String path, String key, int value) {
        String pref = path + "/" + key;
        TEMP_CACHED_VALUES.put(pref, value);
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            node.putInt(key, value);
        }
    }

    public static boolean getBooleanPref(String path, String key, boolean def) {
        String pref = path + "/" + key;
        Boolean result = (Boolean) TEMP_CACHED_VALUES.get(pref);
        if (result != null) {
            return result;
        }
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            result = node.getBoolean(key, def);
        } else {
            result = def;
        }
        TEMP_CACHED_VALUES.put(pref, result);
        return result;
    }

    public static void putBooleanPref(String path, String key, boolean value) {
        String pref = path + "/" + key;
        TEMP_CACHED_VALUES.put(pref, value);
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            node.putBoolean(key, value);
        }
    }

    //
    // Top level preferences
    //
    public static String getSaveDir() {
        return getPref(TASSEL_TOP, TASSEL_SAVE_DIR, TASSEL_SAVE_DIR_DEFAULT);
    }

    public static void putSaveDir(String value) {
        putPref(TASSEL_TOP, TASSEL_SAVE_DIR, value);
    }

    public static String getOpenDir() {
        return getPref(TASSEL_TOP, TASSEL_OPEN_DIR, TASSEL_OPEN_DIR_DEFAULT);
    }

    public static void putOpenDir(String value) {
        putPref(TASSEL_TOP, TASSEL_OPEN_DIR, value);
    }

    public static int getXDim() {
        return getIntPref(TASSEL_TOP, TASSEL_X_DIM, TASSEL_X_DIM_DEFAULT);
    }

    public static void putXDim(int value) {
        putIntPref(TASSEL_TOP, TASSEL_X_DIM, value);
    }

    public static int getYDim() {
        return getIntPref(TASSEL_TOP, TASSEL_Y_DIM, TASSEL_Y_DIM_DEFAULT);
    }

    public static void putYDim(int value) {
        putIntPref(TASSEL_TOP, TASSEL_Y_DIM, value);
    }

    public static boolean getLogSendToConsole() {
        return getBooleanPref(TASSEL_TOP, TASSEL_LOG_SEND_TO_CONSOLE, TASSEL_LOG_SEND_TO_CONSOLE_DEFAULT);
    }

    public static void putLogSendToConsole(boolean value) {
        putBooleanPref(TASSEL_TOP, TASSEL_LOG_SEND_TO_CONSOLE, value);
    }

    public static boolean getLogDebug() {
        return getBooleanPref(TASSEL_TOP, TASSEL_LOG_DEBUG, TASSEL_LOG_DEBUG_DEFAULT);
    }

    public static void putLogDebug(boolean value) {
        putBooleanPref(TASSEL_TOP, TASSEL_LOG_DEBUG, value);
    }

    public static int getLogXDim() {
        return getIntPref(TASSEL_TOP, TASSEL_LOG_X_DIM, TASSEL_LOG_X_DIM_DEFAULT);
    }

    public static void putLogXDim(int value) {
        putIntPref(TASSEL_TOP, TASSEL_LOG_X_DIM, value);
    }

    public static int getLogYDim() {
        return getIntPref(TASSEL_TOP, TASSEL_LOG_Y_DIM, TASSEL_LOG_Y_DIM_DEFAULT);
    }

    public static void putLogYDim(int value) {
        putIntPref(TASSEL_TOP, TASSEL_LOG_Y_DIM, value);
    }

    public static int getMaxThreads() {
        return getIntPref(TASSEL_TOP, TASSEL_MAX_THREADS, TASSEL_MAX_THREADS_DEFAULT);
    }

    public static void putMaxThreads(int value) {
        if (value <= 0) {
            return;
        }
        putIntPref(TASSEL_TOP, TASSEL_MAX_THREADS, value);
    }

    public static String getConfigFile() {
        return getPref(TASSEL_TOP, TASSEL_CONFIG_FILE, TASSEL_CONFIG_FILE_DEFAULT);
    }

    public static void putConfigFile(String value) {
        putPref(TASSEL_TOP, TASSEL_CONFIG_FILE, value);
    }

    //
    // FilterAlignmentPlugin preferences
    //
    public static boolean getExportPluginExportDiploids() {
        return getBooleanPref(EXPORT_PLUGIN_TOP, EXPORT_PLUGIN_EXPORT_DIPLOIDS, EXPORT_PLUGIN_EXPORT_DIPLOIDS_DEFAULT);
    }

    public static void putExportPluginExportDiploids(boolean value) {
        putBooleanPref(EXPORT_PLUGIN_TOP, EXPORT_PLUGIN_EXPORT_DIPLOIDS, value);
    }

    public static boolean getExportPluginIncludeTaxaAnnotations() {
        return getBooleanPref(EXPORT_PLUGIN_TOP, EXPORT_PLUGIN_INCLUDE_TAXA_ANNOTATIONS, EXPORT_PLUGIN_INCLUDE_TAXA_ANNOTATIONS_DEFAULT);
    }

    public static void putExportPluginIncludeTaxaAnnotations(boolean value) {
        putBooleanPref(EXPORT_PLUGIN_TOP, EXPORT_PLUGIN_INCLUDE_TAXA_ANNOTATIONS, value);
    }

    //
    // FilterAlignmentPlugin preferences
    //
    public static double getFilterAlignPluginMinFreq() {
        return getDoublePref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MIN_FREQ, FILTER_ALIGN_PLUGIN_MIN_FREQ_DEFAULT);
    }

    public static void putFilterAlignPluginMinFreq(double value) {
        putDoublePref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MIN_FREQ, value);
    }

    public static double getFilterAlignPluginMaxFreq() {
        return getDoublePref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MAX_FREQ, FILTER_ALIGN_PLUGIN_MAX_FREQ_DEFAULT);
    }

    public static void putFilterAlignPluginMaxFreq(double value) {
        putDoublePref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MAX_FREQ, value);
    }

    public static int getFilterAlignPluginMinCount() {
        return getIntPref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MIN_COUNT, FILTER_ALIGN_PLUGIN_MIN_COUNT_DEFAULT);
    }

    public static void putFilterAlignPluginMinCount(int value) {
        putIntPref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MIN_COUNT, value);
    }

    //
    // FilterTaxaPropertiesPlugin preferences
    //
    public static double getFilterTaxaPropsMinNotMissingFreq() {
        return getDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING, FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING_DEFAULT);
    }

    public static void putFilterTaxaPropsMinNotMissingFreq(double value) {
        putDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING, value);
    }

    public static double getFilterTaxaPropsMinHetFreq() {
        return getDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MIN_HET, FILTER_TAXA_PROPS_PLUGIN_MIN_HET_DEFAULT);
    }

    public static void putFilterTaxaPropsMinHetFreq(double value) {
        putDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MIN_HET, value);
    }

    public static double getFilterTaxaPropsMaxHetFreq() {
        return getDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MAX_HET, FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT);
    }

    public static void putFilterTaxaPropsMaxHetFreq(double value) {
        putDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MAX_HET, value);
    }

    //
    // Alignment preferences
    //
    public static boolean getAlignmentRetainRareAlleles() {
        return getBooleanPref(ALIGNMENT_TOP, ALIGNMENT_RETAIN_RARE_ALLELES, ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT);
    }

    public static void putAlignmentRetainRareAlleles(boolean value) {
        putBooleanPref(ALIGNMENT_TOP, ALIGNMENT_RETAIN_RARE_ALLELES, value);
    }

    //
    // GOBII preferences
    //
    public static String getGOBIIDB() {
        return getPref(GOBII_TOP, GOBII_DB, GOBII_DB_DEFAULT);
    }

    public static void putGOBIIDB(String value) {
        putPref(GOBII_TOP, GOBII_DB, value);
    }

    public static String getGOBIIUser() {
        return getPref(GOBII_TOP, GOBII_USER, GOBII_USER_DEFAULT);
    }

    public static void putGOBIIUser(String value) {
        putPref(GOBII_TOP, GOBII_USER, value);
    }

    public static String getBMSHost() {
        return getPref(GOBII_TOP, BMS_HOST, BMS_HOST_DEFAULT);
    }

    public static void putBMSHost(String value) {
        putPref(GOBII_TOP, BMS_HOST, value);
    }

    public static String getBMSDB() {
        return getPref(GOBII_TOP, BMS_DB, BMS_DB_DEFAULT);
    }

    public static void putBMSDB(String value) {
        putPref(GOBII_TOP, BMS_DB, value);
    }

    public static String getBMSUser() {
        return getPref(GOBII_TOP, BMS_USER, BMS_USER_DEFAULT);
    }

    public static void putBMSUser(String value) {
        putPref(GOBII_TOP, BMS_USER, value);
    }
}
