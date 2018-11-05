/*
 * AbstractPlugin.java
 *
 * Created on December 22, 2006, 5:03 PM
 *
 */
package net.maizegenetics.plugindef;

import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundFill;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.web.WebEngine;
import javafx.scene.web.WebView;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.Modality;
import javafx.stage.Stage;
import javafx.stage.StageStyle;
import javafx.stage.Window;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.io.JSONUtils;
import net.maizegenetics.gui.AlertUtils;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.ReadDistanceMatrix;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.lang.reflect.Field;
import java.math.BigDecimal;
import java.text.NumberFormat;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.regex.Pattern;

/**
 * @author Terry Casstevens
 */
abstract public class AbstractPlugin implements Plugin {

    private static final Logger myLogger = Logger.getLogger(AbstractPlugin.class);

    public static final String DEFAULT_CITATION = "Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. (2007) TASSEL: Software for association mapping of complex traits in diverse samples. Bioinformatics 23:2633-2635.";

    public static final String POSITION_LIST_NONE = "None";
    public static final String TAXA_LIST_NONE = "None";

    private final List<PluginListener> myListeners = new ArrayList<>();
    private final List<Plugin> myInputs = new ArrayList<>();
    private DataSet myCurrentInputData = null;
    private final boolean myIsInteractive;
    private boolean myTrace = false;
    private boolean myThreaded = false;
    protected boolean myWasCancelled = false;

    /**
     * Creates a new instance of AbstractPlugin
     */
    public AbstractPlugin() {
        this(true);
    }

    /**
     * Creates a new instance of AbstractPlugin
     */
    public AbstractPlugin(boolean isInteractive) {
        myIsInteractive = isInteractive;
    }

    @Override
    public DataSet performFunction(DataSet input) {

        LocalDateTime time = LocalDateTime.now();
        String timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
        myLogger.info("Starting " + getClass().getName() + ": time: " + timeStr);

        myCurrentInputData = input;

        try {

            preProcessParameters(input);

            if (isInteractive()) {
                if (!setParametersViaGUI()) {
                    return null;
                }
            }

            checkRequiredParameters();
            postProcessParameters();
            printParameterValues();
            checkParameters();

            DataSet output = processData(input);
            time = LocalDateTime.now();
            timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
            myLogger.info("Finished " + getClass().getName() + ": time: " + timeStr);
            fireDataSetReturned(new PluginEvent(output, getClass()));
            return output;

        } catch (Exception e) {

            if (isInteractive()) {
                myLogger.debug(e.getMessage(), e);
                AlertUtils.showError(e.getMessage() + "\n");
            } else {
                myLogger.debug(e.getMessage(), e);
                printUsage();
                myLogger.error(e.getMessage());
                System.exit(1);
            }
            return null;

        } finally {
            fireProgress(100);
        }

    }

    protected void preProcessParameters(DataSet input) {
        // do nothing
    }

    protected void postProcessParameters() {
        // do nothing
    }

    @Override
    public DataSet processData(DataSet input) {
        throw new UnsupportedOperationException();
    }

    protected List<Field> getParameterFields() {

        List<Field> result = new ArrayList<>();
        Field[] fields = getClass().getDeclaredFields();
        for (Field current : fields) {
            if (current.getType().isAssignableFrom(PluginParameter.class)) {
                current.setAccessible(true);
                result.add(current);
            }
        }

        return result;
    }

    private List<PluginParameter<?>> getParameterInstances() {

        List<PluginParameter<?>> result = new ArrayList<>();
        Field[] fields = getClass().getDeclaredFields();
        for (Field current : fields) {
            if (current.getType().isAssignableFrom(PluginParameter.class)) {
                current.setAccessible(true);
                try {
                    PluginParameter<?> parameter = (PluginParameter) current.get(this);
                    result.add(parameter);
                } catch (Exception e) {
                    e.printStackTrace();
                    throw new IllegalArgumentException("AbstractPlugin: getParameterInstances: problem getting parameter instances");
                }

            }
        }

        return result;
    }

    private Field getParameterField(String key) {

        Field[] fields = getClass().getDeclaredFields();
        for (Field current : fields) {
            if (current.getType().isAssignableFrom(PluginParameter.class)) {
                try {
                    current.setAccessible(true);
                    PluginParameter<?> parameter = (PluginParameter) current.get(this);
                    if (parameter.cmdLineName().equals(key)) {
                        return current;
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    throw new IllegalArgumentException("AbstractPlugin: getParameterField: problem with key: " + key);
                }

            }
        }

        throw new IllegalArgumentException("AbstractPlugin: getParameterField: unknown key: " + key);
    }

    private PluginParameter<?> getParameterInstance(String key) {
        try {
            Field field = getParameterField(key);
            return (PluginParameter) field.get(this);
        } catch (Exception e) {
            return null;
        }
    }

    public static <T> T convert(String input, Class<T> outputClass) {
        try {
            if ((input == null) || (input.length() == 0)) {
                return null;
            } else if (outputClass.isEnum()) {
                return (T) Enum.valueOf((Class<Enum>) outputClass, input);
            } else if (outputClass.isAssignableFrom(String.class)) {
                return (T) input;
            } else if (outputClass.isAssignableFrom(Integer.class)) {
                input = input.replace(",", "");
                return (T) new Integer(new BigDecimal(input).intValueExact());
            } else if (outputClass.isAssignableFrom(Double.class)) {
                input = input.replace(",", "");
                return (T) new Double(new BigDecimal(input).doubleValue());
            } else if (outputClass.isAssignableFrom(List.class)) {
                return (T) getListFromString(input);
            } else if (outputClass.isAssignableFrom(PositionList.class)) {
                return (T) JSONUtils.importPositionListFromJSON(input);
            } else if (outputClass.isAssignableFrom(TaxaList.class)) {
                String test = input.trim().substring(Math.max(0, input.length() - 8)).toLowerCase();
                if (test.endsWith(".json") || test.endsWith(".json.gz")) {
                    return (T) JSONUtils.importTaxaListFromJSON(input);
                } else if (test.endsWith(".txt")) {
                    TaxaListBuilder builder = new TaxaListBuilder();
                    try (BufferedReader br = Utils.getBufferedReader(input)) {
                        String line = br.readLine();
                        Pattern sep = Pattern.compile("\\s+");

                        while (line != null) {
                            line = line.trim();
                            String[] parsedline = sep.split(line);
                            for (int i = 0; i < parsedline.length; i++) {
                                if ((parsedline[i] != null) || (parsedline[i].length() != 0)) {
                                    builder.add(parsedline[i]);
                                }
                            }
                            line = br.readLine();
                        }
                    }

                    return (T) builder.build();
                } else {
                    String[] taxa = input.trim().split(",");
                    return (T) new TaxaListBuilder().addAll(taxa).build();
                }
            } else if (outputClass.isAssignableFrom(DistanceMatrix.class)) {
                return (T) ReadDistanceMatrix.readDistanceMatrix(input);
            } else if (outputClass.isAssignableFrom(Character.class)) {
                if (input.length() != 1) {
                    throw new IllegalArgumentException("Should be one character");
                }
                return (T) new Character(input.charAt(0));
            } else {
                return outputClass.getConstructor(String.class).newInstance(input);
            }
        } catch (Exception nfe) {
            myLogger.debug(nfe.getMessage(), nfe);
            String message = nfe.getMessage();
            if (message == null) {
                throw new IllegalArgumentException("Problem converting: " + input + " to " + Utils.getBasename(outputClass.getName()));
            } else {
                throw new IllegalArgumentException(message + " Problem converting: " + input + " to " + Utils.getBasename(outputClass.getName()));
            }
        }
    }

    private static List<String> getListFromString(String str) {

        if ((str == null) || (str.length() == 0)) {
            return null;
        }
        String[] tokens = str.split(",");
        List<String> result = new ArrayList<>();
        for (String current : tokens) {
            current = current.trim();
            if (current.length() != 0) {
                result.add(current);
            }
        }
        return result;

    }

    @Override
    public void setParameters(String[] args) {

        if (args != null) {

            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (args[i].startsWith("-")) {
                    arg = arg.substring(1);
                    PluginParameter<?> parameter = getParameterInstance(arg);
                    if (parameter == null) {
                        myLogger.error("Unrecognized argument: " + args[i]);
                        printUsage();
                        System.exit(1);
                    }
                    if ((i == args.length - 1) || (args[i + 1]).startsWith("-")) {
                        if (parameter.valueType().isAssignableFrom(Boolean.class)) {
                            setParameter(arg, Boolean.TRUE);
                        } else if (Number.class.isAssignableFrom(parameter.valueType())) {
                            setParameter(arg, args[i + 1]);
                            i++;
                        } else {
                            myLogger.error("Parameter requires a value: " + args[i]);
                            printUsage();
                            System.exit(1);
                        }
                    } else {
                        setParameter(arg, args[i + 1]);
                        i++;
                    }
                } else {
                    myLogger.error("Argument expected to start with dash(-): " + args[i]);
                    printUsage();
                    System.exit(1);
                }
            }

        }

    }

    private void setFieldsToConfigParameters(Map<String, Node> parameterFields) {

        final List<PluginParameter<?>> parameterInstances = getParameterInstances();
        if (parameterInstances.isEmpty()) {
            return;
        }

        for (final PluginParameter<?> current : parameterInstances) {
            Node component = parameterFields.get(current.cmdLineName());
            setFieldToConfigParameters(component, current);
        }

    }

    private void setFieldToConfigParameters(Node component, PluginParameter<?> parameter) {

        Optional<String> configValue = ParameterCache.value(this, parameter.cmdLineName());
        if (!configValue.isPresent()) {
            return;
        }
        try {
            if (component instanceof TextField) {
                ((TextField) component).setText(configValue.get());
            } else if (component instanceof CheckBox) {
                Boolean value = convert(configValue.get(), Boolean.class);
                ((CheckBox) component).setSelected(value);
            } else if (component instanceof ComboBox) {
                Object value = convert(configValue.get(), parameter.valueType());
                ((ComboBox) component).getSelectionModel().select(value);
            }
        } catch (Exception e) {
            myLogger.warn("setFieldToConfigParameters: problem with configuration key: " + this.getClass().getName() + "." + parameter.cmdLineName() + "  value: " + configValue.get() + "\n" + e.getMessage());
        }

    }

    public void setConfigParameters() {

        if (ParameterCache.hasValues()) {

            for (PluginParameter<?> parameter : getParameterInstances()) {
                Optional<String> value = ParameterCache.value(this, parameter.cmdLineName());
                if (value.isPresent()) {
                    try {
                        setParameter(parameter, value.get());
                    } catch (Exception e) {
                        myLogger.warn("setConfigParameters: problem with configuration key: " + this.getClass().getName() + "." + parameter.cmdLineName() + "  value: " + value.get() + "\n" + e.getMessage());
                    }
                }
            }

        }

    }

    private void checkRequiredParameters() {

        List<String> cmdLineNames = new ArrayList<>();
        for (PluginParameter<?> current : getParameterInstances()) {
            if (cmdLineNames.contains(current.cmdLineName())) {
                if (isInteractive()) {
                    throw new IllegalStateException(current.cmdLineName() + " exist multiple times for this plugin.");
                } else {
                    myLogger.error("-" + current.cmdLineName() + " exist multiple times for this plugin.\n");
                    printUsage();
                    System.exit(1);
                }
            } else {
                cmdLineNames.add(current.cmdLineName());
            }

            if (current.required()) {
                if (current.isEmpty()) {
                    if (isInteractive()) {
                        throw new IllegalStateException(current.guiName() + " must be defined.");
                    } else {
                        myLogger.error("-" + current.cmdLineName() + " is required.\n");
                        printUsage();
                        System.exit(1);
                    }
                }
            }
        }

    }

    /**
     * Verification checks of parameters.
     */
    private void checkParameters() {

        for (PluginParameter<?> current : getParameterInstances()) {

            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.IN_FILE) {
                if (!current.isEmpty()) {
                    String filename = current.value().toString();
                    File theFile = new File(filename);
                    if (!theFile.exists()) {
                        if (isInteractive()) {
                            throw new IllegalStateException(current.guiName() + ": " + filename + " doesn't exist.");
                        } else {
                            myLogger.error("-" + current.cmdLineName() + ": " + filename + " doesn't exist\n");
                            printUsage();
                            System.exit(1);
                        }
                    }
                    if (!theFile.isFile()) {
                        if (isInteractive()) {
                            throw new IllegalStateException(current.guiName() + ": " + filename + " isn't a file.");
                        } else {
                            myLogger.error("-" + current.cmdLineName() + ": " + filename + " isn't a file\n");
                            printUsage();
                            System.exit(1);
                        }
                    }
                }
            }

            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.OUT_FILE) {
                if (!current.isEmpty()) {
                    String filename = current.value().toString();
                    String outFolder = Utils.getDirectory(filename);
                    File outDir = new File(outFolder);
                    if (!outDir.isDirectory()) {
                        if (isInteractive()) {
                            throw new IllegalStateException(current.guiName() + ": Output Directory: " + outFolder + " doesn't exist.");
                        } else {
                            myLogger.error("-" + current.cmdLineName() + ": Output Directory: " + outFolder + " doesn't exist\n");
                            printUsage();
                            System.exit(1);
                        }
                    }
                }
            }

            if ((current.parameterType() == PluginParameter.PARAMETER_TYPE.IN_DIR)
                    || (current.parameterType() == PluginParameter.PARAMETER_TYPE.OUT_DIR)) {
                if (!current.isEmpty()) {
                    String dirname = current.value().toString();
                    File directory = new File(dirname);
                    if (!directory.isDirectory()) {
                        if (isInteractive()) {
                            throw new IllegalStateException(current.guiName() + ": Directory: " + dirname + " doesn't exist.");
                        } else {
                            myLogger.error("-" + current.cmdLineName() + ": Directory: " + dirname + " doesn't exist\n");
                            printUsage();
                            System.exit(1);
                        }
                    }
                }
            }

        }

    }

    protected void printParameterValues() {
        List<PluginParameter<?>> parameters = getParameterInstances();
        if ((parameters == null) || (parameters.isEmpty())) {
            return;
        }
        StringBuilder builder = new StringBuilder();
        builder.append("\n");
        builder.append(Utils.getBasename(getClass().getName()));
        builder.append(" Parameters\n");
        for (PluginParameter<?> current : parameters) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                continue;
            }
            builder.append(current.cmdLineName());
            builder.append(": ");
            Object value = current.value();
            if (value instanceof PositionList) {
                builder.append(((PositionList) value).numberOfSites());
                builder.append(" positions");
            } else if (value instanceof List) {
                builder.append((Arrays.toString(((List) value).toArray())));
            } else if (value instanceof DistanceMatrix) {
                DistanceMatrix temp = (DistanceMatrix) value;
                builder.append(temp.getColumnCount());
                builder.append(" columns x ");
                builder.append(temp.getRowCount());
                builder.append(" rows");
            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.PASSWORD) {
                builder.append("?????");
            } else {
                builder.append(value);
            }
            builder.append("\n");
        }
        myLogger.info(builder.toString());
    }

    private void printUsage() {

        StringBuilder builder = new StringBuilder();
        String description = pluginDescription();
        if (description != null) {
            builder.append("\n");
            builder.append(Utils.getBasename(getClass().getName())).append(" Description...\n");
            builder.append(description);
            builder.append("\n");
        }
        builder.append("\nUsage:\n");
        builder.append(Utils.getBasename(getClass().getName())).append(" <options>\n");
        for (PluginParameter<?> current : getParameterInstances()) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                continue;
            }
            builder.append("-");
            builder.append(current.cmdLineName());
            builder.append(" ");
            if (current.valueType().isAssignableFrom(Boolean.class)) {
                builder.append("<true | false>");
            } else {
                builder.append("<");
                builder.append(current.guiName());
                builder.append(">");
            }
            builder.append(" : ");
            builder.append(current.description());
            if (current.hasRange()) {
                builder.append(" ");
                builder.append(current.rangeToString());
            }
            if (current.defaultValue() != null) {
                builder.append(" (Default: ");
                builder.append(current.defaultValue());
                builder.append(")");
            }
            if (current.required()) {
                builder.append(" (required)");
            }
            builder.append("\n");
        }

        myLogger.info(builder.toString());
    }

    @Override
    public String getUsage() {

        StringBuilder builder = new StringBuilder();
        builder.append(Utils.getBasename(getClass().getName()));
        builder.append("\n");
        String description = pluginDescription();
        if (description != null) {
            builder.append("\nDescription: ");
            builder.append(description);
            builder.append("\n\n");
        }
        for (PluginParameter<?> current : getParameterInstances()) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                continue;
            }
            builder.append("\n");
            builder.append(current.guiName());
            builder.append(" : ");
            builder.append(current.description());
            if (current.hasRange()) {
                builder.append(" ");
                builder.append(current.rangeToString());
            }
            if (current.defaultValue() != null) {
                builder.append(" (Default: ");
                builder.append(current.defaultValue());
                builder.append(")");
            }
            if (current.required()) {
                builder.append(" (required)");
            }
            builder.append("\n");
        }

        return builder.toString();
    }

    @Override
    public Object getParameter(Enum key) {
        return getParameterInstance(key.toString()).value();
    }

    @Override
    public Object getParameter(String key) {
        return getParameterInstance(key).value();
    }

    @Override
    public Plugin setParameter(PluginParameter<?> param, Object value) {
        if (value == null) {
            setParameter(param.cmdLineName(), value);
        } else if (value instanceof String) {
            setParameter(param.cmdLineName(), (String) value);
        } else {
            setParameter(param.cmdLineName(), value);
        }
        return this;
    }

    @Override
    public Plugin setParameter(String key, Object value) {

        try {

            Field field = getParameterField(key);
            PluginParameter parameter = (PluginParameter) field.get(this);
            if (parameter == null) {
                throw new IllegalArgumentException("setParameter: Unknown Parameter: " + key);
            }
            if ((parameter.hasRange()) && (!parameter.acceptsValue(value))) {
                throw new IllegalArgumentException("setParameter: " + parameter.cmdLineName() + " value: " + value.toString() + " outside range: " + parameter.rangeToString());
            }
            PluginParameter newParameter = new PluginParameter<>(parameter, value);
            field.set(this, newParameter);

        } catch (Exception e) {
            if (isInteractive()) {
                try {
                    throw e;
                } catch (IllegalAccessException ex) {
                    myLogger.error(ex.getMessage(), ex);
                }
            } else {
                myLogger.error(key + ": " + e.getMessage());
                printUsage();
                myLogger.debug(e.getMessage(), e);
                System.exit(1);
            }
        }

        return this;
    }

    @Override
    public Plugin setParameter(String key, String value) {

        try {
            PluginParameter parameter = getParameterInstance(key);
            return setParameter(key, convert(value, parameter.valueType()));
        } catch (Exception e) {
            if (isInteractive()) {
                throw new IllegalArgumentException(getParameterInstance(key).guiName() + ": " + e.getMessage());
            } else {
                myLogger.error(key + ": " + e.getMessage());
                printUsage();
                System.exit(1);
            }
        }
        return this;
    }

    private static final int TEXT_FIELD_WIDTH = 25;

    private boolean parametersAreSet = true;

    /**
     * Generates dialog based on this plugins define parameters.
     *
     * @return true if OK clicked, false if canceled
     */
    private boolean setParametersViaGUI() {

        final List<PluginParameter<?>> parameterInstances = getParameterInstances();
        if (parameterInstances.isEmpty()) {
            return true;
        }

        final Stage dialog = new Stage();
        dialog.initModality(Modality.APPLICATION_MODAL);
        dialog.setOnCloseRequest(event -> {
            parametersAreSet = false;
            dialog.close();
        });

        final Map<String, Node> parameterFields = new HashMap<>();

        parametersAreSet = true;

        Button okButton = new Button("Ok");
        okButton.setOnAction(event -> {
            try {
                for (final PluginParameter<?> current : parameterInstances) {
                    Node component = parameterFields.get(current.cmdLineName());
                    if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                        // do nothing
                    } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.GENOTYPE_TABLE) {
                        GenotypeWrapper input = (GenotypeWrapper) ((ComboBox) component).getSelectionModel().getSelectedItem();
                        if (input != null) {
                            setParameter(current.cmdLineName(), input.myObj);
                        }
                    } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.POSITION_LIST) {
                        if (component instanceof ComboBox) {
                            Object temp = ((ComboBox) component).getSelectionModel().getSelectedItem();
                            if (temp == POSITION_LIST_NONE) {
                                setParameter(current.cmdLineName(), null);
                            } else {
                                setParameter(current.cmdLineName(), ((Datum) temp).getData());
                            }
                        } else {
                            String input = ((TextField) component).getText().trim();
                            setParameter(current.cmdLineName(), input);
                        }
                    } else if (TaxaList.class.isAssignableFrom(current.valueType())) {
                        if (component instanceof ComboBox) {
                            Object temp = ((ComboBox) component).getSelectionModel().getSelectedItem();
                            if (temp == TAXA_LIST_NONE) {
                                setParameter(current.cmdLineName(), null);
                            } else {
                                setParameter(current.cmdLineName(), ((Datum) temp).getData());
                            }
                        } else {
                            String input = ((TextField) component).getText().trim();
                            setParameter(current.cmdLineName(), input);
                        }
                    } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.DISTANCE_MATRIX) {
                        if (component instanceof ComboBox) {
                            Object temp = ((ComboBox) component).getSelectionModel().getSelectedItem();
                            if (temp == null) {
                                throw new IllegalArgumentException("setParametersViaGUI: must specify a distance matrix.");
                            } else {
                                setParameter(current.cmdLineName(), ((Datum) temp).getData());
                            }
                        } else {
                            String input = ((TextField) component).getText().trim();
                            setParameter(current.cmdLineName(), input);
                        }
                    } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.OBJECT_LIST_SINGLE_SELECT) {
                        Object selectedObjects = ((ComboBox<?>) component).getSelectionModel().getSelectedItem();
                        setParameter(current.cmdLineName(), selectedObjects);
                    } else if (component instanceof TextField) {
                        String input = ((TextField) component).getText().trim();
                        setParameter(current.cmdLineName(), input);
                    } else if (component instanceof CheckBox) {
                        if (((CheckBox) component).isSelected()) {
                            setParameter(current.cmdLineName(), Boolean.TRUE);
                        } else {
                            setParameter(current.cmdLineName(), Boolean.FALSE);
                        }
                    } else if (component instanceof ComboBox) {
                        Enum temp = (Enum) ((ComboBox) component).getSelectionModel().getSelectedItem();
                        setParameter(current.cmdLineName(), temp);
                    }
                }
            } catch (Exception ex) {
                myLogger.debug(ex.getMessage(), ex);
                StringBuilder builder = new StringBuilder();
                builder.append("Problem Setting Parameters: ");
                builder.append("\n");
                builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(ex), 50));
                String str = builder.toString();
                DialogUtils.showError(str, null);
                return;
            }
            dialog.close();
        });

        Button cancelButton = new Button("Cancel");
        cancelButton.setOnAction(event -> {
            parametersAreSet = false;
            dialog.close();
        });

        Button defaultsButton = new Button("Defaults");
        defaultsButton.setOnAction(event -> {
            setFieldsToDefault(parameterFields);
            setFieldsToConfigParameters(parameterFields);
        });

        Button userManualButton = new Button("User Manual");
        userManualButton.setOnAction(event -> {
            try {
                initUserManualWindow();
                webEngine.load(pluginUserManualURL());
                userManual.show();
            } catch (Exception ex) {
                myLogger.debug(ex.getMessage(), ex);
            }
        });

        VBox panel = new VBox(15.0);
        panel.setPadding(new Insets(20.0));
        panel.setAlignment(Pos.CENTER);

        boolean show_citation = !DEFAULT_CITATION.equals(getCitation());
        TextArea citationText = null;
        if (show_citation) {
            citationText = new TextArea();
            citationText.setPadding(new Insets(5.0));
            //citationText.setContentType("text/html");
            citationText.setEditable(false);
            ScrollPane scroll = new ScrollPane(citationText);
            scroll.setHbarPolicy(ScrollPane.ScrollBarPolicy.NEVER);
            scroll.setVbarPolicy(ScrollPane.ScrollBarPolicy.AS_NEEDED);
            scroll.setPrefSize(scroll.getWidth(), 45.0);
            panel.getChildren().add(scroll);
        }

        for (final PluginParameter<?> current : getParameterInstances()) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.GENOTYPE_TABLE) {
                Datum datum = getGenotypeTable();
                ComboBox<GenotypeWrapper> menu = new ComboBox<>();
                if (datum != null) {
                    String name = datum.getName();
                    GenotypeTable table;
                    if (datum.getData() instanceof GenotypeTable) {
                        table = (GenotypeTable) datum.getData();
                    } else if (datum.getData() instanceof GenotypePhenotype) {
                        table = ((GenotypePhenotype) datum.getData()).genotypeTable();
                    } else {
                        throw new IllegalStateException("AbstractPlugin: setParametersViaGUI: unknown GenotypeTable type: " + datum.getData().getClass().getName());
                    }

                    if (current.acceptsValue(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype) && table.hasGenotype()) {
                        menu.getItems().add(new GenotypeWrapper(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, "Genotype (" + name + ")"));
                    }
                    if (current.acceptsValue(GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability) && table.hasReferenceProbablity()) {
                        menu.getItems().add(new GenotypeWrapper(GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability, "Reference Probability (" + name + ")"));
                    }
                    if (current.acceptsValue(GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability) && table.hasAlleleProbabilities()) {
                        menu.getItems().add(new GenotypeWrapper(GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability, "Allele Probability (" + name + ")"));
                    }
                    if (current.acceptsValue(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Depth) && table.hasDepth()) {
                        menu.getItems().add(new GenotypeWrapper(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Depth, "Depth (" + name + ")"));
                    }
                    if (current.acceptsValue(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Dosage) && table.hasDosage()) {
                        menu.getItems().add(new GenotypeWrapper(GenotypeTable.GENOTYPE_TABLE_COMPONENT.Dosage, "Dosage (" + name + ")"));
                    }
                    menu.getSelectionModel().clearAndSelect(0);
                }
                createEnableDisableAction(current, parameterFields, menu);
                FlowPane temp = new FlowPane(Orientation.HORIZONTAL);
                temp.getChildren().add(new javafx.scene.control.Label(current.guiName()));
                temp.getChildren().add(menu);
                //temp.setToolTipText(getToolTip(current));
                panel.getChildren().add(temp);
                parameterFields.put(current.cmdLineName(), menu);
            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.POSITION_LIST) {
                Datum datum = getPositionList();
                if (datum != null) {
                    ComboBox menu = new ComboBox();
                    menu.getItems().add(POSITION_LIST_NONE);
                    menu.getItems().add(datum);
                    menu.getSelectionModel().select(0);
                    createEnableDisableAction(current, parameterFields, menu);
                    FlowPane temp = new FlowPane(Orientation.HORIZONTAL);
                    temp.getChildren().add(new Label(current.guiName()));
                    temp.getChildren().add(menu);
                    //temp.setToolTipText(getToolTip(current));
                    panel.getChildren().add(temp);
                    parameterFields.put(current.cmdLineName(), menu);
                } else {
                    TextField field = new TextField();
                    Button browse = getOpenFile(dialog, field);
                    HBox line = getLine(current.guiName(), field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new Node[]{field, browse}, field);
                    panel.getChildren().add(line);
                    parameterFields.put(current.cmdLineName(), field);
                }
            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.DISTANCE_MATRIX) {
                List<Datum> matrices = getDistanceMatrices();
                if (!matrices.isEmpty()) {
                    ComboBox menu = new ComboBox();
                    for (Datum matrix : matrices) {
                        menu.getItems().add(matrix);
                    }
                    menu.getSelectionModel().select(0);
                    createEnableDisableAction(current, parameterFields, menu);
                    FlowPane temp = new FlowPane(Orientation.HORIZONTAL);
                    temp.getChildren().add(new Label(current.guiName()));
                    temp.getChildren().add(menu);
                    //temp.setToolTipText(getToolTip(current));
                    panel.getChildren().add(temp);
                    parameterFields.put(current.cmdLineName(), menu);
                } else {
                    TextField field = new TextField();
                    Button browse = getOpenFile(dialog, field);
                    HBox line = getLine(current.guiName(), field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new Node[]{field, browse}, field);
                    panel.getChildren().add(line);
                    parameterFields.put(current.cmdLineName(), field);
                }
            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.OBJECT_LIST_SINGLE_SELECT) {
                GridPane temp = new GridPane();
                temp.setAlignment(Pos.CENTER);
                temp.setHgap(20.0);
                temp.add(new Label(current.guiName()), 0, 0);
                ComboBox list = new ComboBox();
                list.getItems().addAll(current.possibleValues().toArray());
                list.getSelectionModel().select(0);
                createEnableDisableAction(current, parameterFields, list);
                temp.add(list, 1, 0);
                //listPanel.setToolTipText(getToolTip(current));
                panel.getChildren().add(temp);
                parameterFields.put(current.cmdLineName(), list);
            } else if (current.valueType().isEnum()) {
                ComboBox menu = new ComboBox();
                Object[] values = current.valueType().getEnumConstants();
                for (Object item : values) {
                    menu.getItems().add(item);
                }
                menu.getSelectionModel().select(current.value());
                createEnableDisableAction(current, parameterFields, menu);
                GridPane temp = new GridPane();
                temp.setAlignment(Pos.CENTER);
                temp.setHgap(20.0);
                //HBox.setHgrow(temp, Priority.ALWAYS);
                temp.add(new Label(current.guiName()), 0, 0);
                temp.add(menu, 1, 0);
                //temp.getChildren().add(new Label(current.guiName()));
                //temp.getChildren().add(menu);
                //temp.setToolTipText(getToolTip(current));
                panel.getChildren().add(temp);
                parameterFields.put(current.cmdLineName(), menu);
            } else if (Boolean.class.isAssignableFrom(current.valueType())) {
                CheckBox check = new CheckBox(current.guiName());
                check.setTooltip(new Tooltip(getToolTip(current)));
                check.setSelected((Boolean) current.value());
                createEnableDisableAction(current, parameterFields, check);
                panel.getChildren().add(check);
                parameterFields.put(current.cmdLineName(), check);
            } else if (TaxaList.class.isAssignableFrom(current.valueType())) {

                List<Datum> datums = getTaxaListDatum();
                if (datums != null) {
                    ComboBox menu = new ComboBox();
                    menu.getItems().add(TAXA_LIST_NONE);
                    for (Datum datum : datums) {
                        menu.getItems().add(datum);
                    }
                    menu.getSelectionModel().select(0);
                    createEnableDisableAction(current, parameterFields, menu);
                    FlowPane temp = new FlowPane(Orientation.HORIZONTAL);
                    temp.getChildren().add(new Label(current.guiName()));
                    temp.getChildren().add(menu);
                    //temp.setToolTipText(getToolTip(current));
                    panel.getChildren().add(temp);
                    parameterFields.put(current.cmdLineName(), menu);
                } else {
                    TaxaList taxa = getTaxaList();
                    TextField field;
                    if (taxa == null) {
                        field = new TextField();
                    } else {
                        field = new TextField();
                    }
                    if (current.value() != null) {
                        field.setText(current.value().toString());
                    }
                    HBox taxaPanel = getTaxaListPanel(current.guiName(), field, current.description(), dialog, taxa);
                    panel.getChildren().add(taxaPanel);
                    parameterFields.put(current.cmdLineName(), field);
                }

            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.SITE_NAME_LIST) {
                PositionList positions = getSiteNameList();
                TextField field;
                if (positions == null) {
                    field = new TextField();
                    field.setMinWidth(TEXT_FIELD_WIDTH);
                } else {
                    field = new TextField();
                    field.setMinWidth(TEXT_FIELD_WIDTH - 7);
                }
                if (current.value() != null) {
                    field.setText(current.value().toString());
                }
                HBox positionsPanel = new HBox();
                //JPanel positionsPanel = getPositionListPanel(current.guiName(), field, current.description(), dialog, positions);
                List<Node> componentList = new ArrayList<>();
                //for (Component component : positionsPanel.getComponents()) {
                //    if (component instanceof JComponent) {
                //        componentList.add((JComponent) component);
                //    }
                //}
                createEnableDisableAction(current, parameterFields, componentList.toArray(new Node[0]), field);
                panel.getChildren().add(positionsPanel);
                parameterFields.put(current.cmdLineName(), field);
            } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                FlowPane temp = new FlowPane(Orientation.HORIZONTAL);
                Label label = new Label(current.guiName());
                label.setFont(new javafx.scene.text.Font("Dialog", 14.0));
                //label.setFont(new Font("Dialog", Font.BOLD, 14));
                temp.getChildren().add(label);
                panel.getChildren().add(temp);
            } else {
                final TextField field;
                if (current.parameterType() == PluginParameter.PARAMETER_TYPE.PASSWORD) {
                    field = new PasswordField();
                    field.setMinWidth(TEXT_FIELD_WIDTH);
                } else if (current.parameterType() != PluginParameter.PARAMETER_TYPE.NA) {
                    field = new TextField();
                    field.setMinWidth(TEXT_FIELD_WIDTH - 8);
                } else {
                    field = new TextField();
                    field.setMinWidth(TEXT_FIELD_WIDTH);
                }

                if (current.value() != null) {
                    if (Integer.class.isAssignableFrom(current.valueType())) {
                        field.setText(NumberFormat.getInstance().format(current.value()));
                    } else {
                        field.setText(current.value().toString());
                    }
                }

                field.focusedProperty().addListener((observable, oldValue, newValue) -> {
                    String input = field.getText().trim();
                    try {
                        if (!current.acceptsValue(input)) {
                            Alert alert = new Alert(Alert.AlertType.ERROR, current.guiName() + " range: " + current.rangeToString(), ButtonType.OK);
                            alert.showAndWait();
                            field.setText(getParameterInstance(current.cmdLineName()).value().toString());
                        }
                        if (Integer.class.isAssignableFrom(current.valueType())) {
                            Integer temp = convert(field.getText(), Integer.class);
                            if (temp == null) {
                                field.setText(null);
                            } else {
                                field.setText(NumberFormat.getInstance().format(temp.intValue()));
                            }
                        }
                    } catch (Exception ex) {
                        myLogger.debug(ex.getMessage(), ex);
                        Alert alert = new Alert(Alert.AlertType.ERROR, current.guiName() + ": " + ex.getMessage(), ButtonType.OK);
                        alert.showAndWait();
                        field.setText(getParameterInstance(current.cmdLineName()).value().toString());
                    }
                });

                String label = null;
                if (current.required()) {
                    label = current.guiName() + "*";
                } else {
                    label = current.guiName();
                }

                HBox line;
                if (current.parameterType() == PluginParameter.PARAMETER_TYPE.IN_FILE) {
                    Button browse = getOpenFile(dialog, field);
                    line = getLine(label, field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new Node[]{field, browse}, field);
                } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.OUT_FILE) {
                    Button browse = getSaveFile(dialog, field);
                    line = getLine(label, field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new Node[]{field, browse}, field);
                } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.IN_DIR) {
                    Button browse = getOpenDir(dialog, field);
                    line = getLine(label, field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new Node[]{field, browse}, field);
                } else if (current.parameterType() == PluginParameter.PARAMETER_TYPE.OUT_DIR) {
                    Button browse = getSaveDir(dialog, field);
                    line = getLine(label, field, browse, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, new Node[]{field, browse}, field);
                } else {
                    line = getLine(label, field, null, getToolTip(current));
                    createEnableDisableAction(current, parameterFields, field);
                }
                panel.getChildren().add(line);

                parameterFields.put(current.cmdLineName(), field);
            }
        }

        TabPane tabbedPane = new TabPane();
        ScrollPane scroll = new ScrollPane(panel);
        scroll.setMinSize(500.0, 200.0);
        scroll.setFitToWidth(true);
        scroll.setFitToHeight(true);
        scroll.setVbarPolicy(ScrollPane.ScrollBarPolicy.AS_NEEDED);
        scroll.setHbarPolicy(ScrollPane.ScrollBarPolicy.NEVER);
        tabbedPane.getTabs().add(new Tab(getButtonName(), scroll));

        FlowPane pnlButtons = new FlowPane(Orientation.HORIZONTAL);
        pnlButtons.setAlignment(Pos.CENTER);
        pnlButtons.setHgap(15.0);
        pnlButtons.setVgap(15.0);
        pnlButtons.getChildren().add(okButton);
        pnlButtons.getChildren().add(cancelButton);
        pnlButtons.getChildren().add(defaultsButton);
        pnlButtons.getChildren().add(userManualButton);

        BorderPane main = new BorderPane();
        main.setPadding(new Insets(15.0));
        main.setCenter(tabbedPane);
        main.setBottom(pnlButtons);
        BorderPane.setMargin(pnlButtons, new Insets(20.0, 0.0, 10.0, 0.0));

        WebView browser = new WebView();
        browser.prefHeightProperty().bind(dialog.heightProperty());
        browser.prefWidthProperty().bind(dialog.widthProperty());
        WebEngine webEngine = browser.getEngine();
        webEngine.loadContent(getUsageHTML());
        tabbedPane.getTabs().add(new Tab("Help", browser));
        if (show_citation) {
            citationText.setText(getCitationHTML((int) (dialog.getWidth() / 9.0)));
            //dialog.setMinimumSize(null);
        }
        //Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        //if (screenSize.getHeight() - 125 < dialog.getHeight()) {
        //dialog.setSize(Math.max(dialog.getWidth(), 550), (int) screenSize.getHeight() - 125);
        //} else {
        //dialog.setSize(Math.max(dialog.getWidth(), 550), Math.max(dialog.getHeight(), 250));
        //}

        dialog.setResizable(false);

        Scene scene = new Scene(main);
        scene.getStylesheets().add("/javafx/AppStyle.css");
        dialog.setScene(scene);

        dialog.sizeToScene();

        dialog.showAndWait();
        return parametersAreSet;

    }

    private static Stage userManual;
    private static WebEngine webEngine;

    private static void initUserManualWindow() {

        if (userManual != null) {
            return;
        }

        userManual = new Stage(StageStyle.DECORATED);
        userManual.initModality(Modality.NONE);
        userManual.setTitle("User Manual");
        userManual.setResizable(false);

        WebView browser = new WebView();
        browser.setPrefWidth(Double.MAX_VALUE);
        VBox.setVgrow(browser, Priority.ALWAYS);
        HBox.setHgrow(browser, Priority.ALWAYS);
        webEngine = browser.getEngine();

        Button close = new Button("Close");
        close.setOnAction(event1 -> userManual.close());

        VBox main = new VBox();
        main.setAlignment(Pos.CENTER);
        main.setPadding(new Insets(10.0));
        main.setSpacing(10.0);
        main.getChildren().add(browser);
        main.getChildren().add(close);

        userManual.setScene(new Scene(main, 750, 1000));

    }

    private void setFieldsToDefault(Map<String, Node> parameterFields) {

        final List<PluginParameter<?>> parameterInstances = getParameterInstances();
        if (parameterInstances.isEmpty()) {
            return;
        }

        for (final PluginParameter<?> current : parameterInstances) {
            Node component = parameterFields.get(current.cmdLineName());
            setFieldToDefault(component, current);
        }

    }

    private void setFieldToDefault(Node component, PluginParameter<?> current) {
        if (component instanceof TextField) {
            Object defaultValue = current.defaultValue();
            if (defaultValue == null) {
                ((TextField) component).setText(null);
            } else {
                ((TextField) component).setText(defaultValue.toString());
            }
            setParameter(current.cmdLineName(), defaultValue);
        } else if (component instanceof CheckBox) {
            Boolean value = (Boolean) current.defaultValue();
            ((CheckBox) component).setSelected(value);
            setParameter(current.cmdLineName(), value);
        } else if (component instanceof ComboBox) {
            ((ComboBox) component).getSelectionModel().select(current.defaultValue());
            setParameter(current.cmdLineName(), current.defaultValue());
        }
    }

    @Override
    public void setParametersToDefault() {

        final List<PluginParameter<?>> parameterInstances = getParameterInstances();
        if (parameterInstances.isEmpty()) {
            return;
        }

        for (final PluginParameter<?> current : parameterInstances) {
            setParameter(current.cmdLineName(), current.defaultValue());
        }

    }

    private String getCitationHTML(int lineWidth) {
        String citation = getCitation();
        int count = 10;
        StringBuilder builder = new StringBuilder();
        builder.append("<html><center>Citation: ");
        for (int i = 0, n = citation.length(); i < n; i++) {
            count++;
            if (citation.charAt(i) == '\n') {
                builder.append("<br>");
                count = 0;
            } else if ((count > lineWidth) && (citation.charAt(i) == ' ')) {
                builder.append("<br>");
                count = 0;
            } else {
                builder.append(citation.charAt(i));
            }
        }
        builder.append("</center></html>");
        return builder.toString();
    }

    public String getUsageHTML() {

        StringBuilder builder = new StringBuilder();
        builder.append("<html><center><strong>");
        builder.append(Utils.getBasename(getClass().getName()));
        builder.append("</strong>");
        String description = pluginDescription();
        if (description != null) {
            builder.append("<br><br><strong>Description:</strong> ");
            builder.append(description);
        }
        builder.append("<br><br>");
        builder.append("<table border='1'>");
        builder.append("<tr><th>Parameter</th><th>Description</th><th>Values</th><th>Default</th></tr>");
        for (PluginParameter<?> current : getParameterInstances()) {
            if (current.parameterType() == PluginParameter.PARAMETER_TYPE.LABEL) {
                continue;
            }
            builder.append("<tr>");

            builder.append("<th>");
            if (current.required()) {
                builder.append("<font color='red'>");
            }
            builder.append(current.guiName());
            if (current.required()) {
                builder.append("</font>");
            }
            builder.append("</th>");

            builder.append("<td>");
            builder.append(current.description());
            builder.append("</td>");

            builder.append("<td>");
            if (current.hasPossibleValues()) {
                String range = current.possibleValuesString(true);
                if ((range.charAt(0) == '[') && (range.charAt(range.length() - 1) == ']')) {
                    range = range.substring(1, range.length() - 1);
                }
                StringBuilder buildRange = new StringBuilder();
                for (char rangeChr : range.toCharArray()) {
                    if (rangeChr == '_') {
                        buildRange.append("\n");
                    }
                    buildRange.append(rangeChr);
                }
                builder.append(buildRange.toString());
            } else if (current.hasRange()) {
                String range = current.rangeToString(true);
                if ((range.charAt(0) == '[') && (range.charAt(range.length() - 1) == ']')) {
                    range = range.substring(1, range.length() - 1);
                }
                StringBuilder buildRange = new StringBuilder();
                for (char rangeChr : range.toCharArray()) {
                    if (rangeChr == '_') {
                        buildRange.append("\n");
                    }
                    buildRange.append(rangeChr);
                }
                builder.append(buildRange.toString());
            } else if (current.valueType().isAssignableFrom(Boolean.class)) {
                builder.append("true, false");
            }
            builder.append("</td>");

            builder.append("<td>");
            if (current.defaultValue() != null) {
                String defaultValue = current.defaultValue().toString();
                StringBuilder buildDefault = new StringBuilder();
                for (char defaultChr : defaultValue.toCharArray()) {
                    if (defaultChr == '_') {
                        buildDefault.append("\n");
                    }
                    buildDefault.append(defaultChr);
                }
                builder.append(buildDefault.toString());
            }
            builder.append("</td>");

            builder.append("</tr>");
        }
        builder.append("</table>");

        builder.append("</center>");

        builder.append("<br><font color='red'>* parameters in red are required</font>");

        builder.append("</html>");

        return builder.toString();
    }

    private void createEnableDisableAction(PluginParameter<?> current, Map<String, Node> parameterFields, final Node component) {
        createEnableDisableAction(current, parameterFields, new Node[]{component}, component);
    }

    private void createEnableDisableAction(final PluginParameter<?> current, Map<String, Node> parameterFields, final Node[] components, final Node input) {

        if (current.dependentOnParameter() != null) {
            Node depends = parameterFields.get(current.dependentOnParameter().cmdLineName());
            if (depends instanceof CheckBox) {
                final CheckBox checkBox = (CheckBox) depends;

                for (Node component : components) {
                    if (checkBox.isSelected() == (Boolean) current.dependentOnParameterValue()[0]) {
                        component.setDisable(false);
                    } else {
                        component.setDisable(true);
                    }
                }

                checkBox.setOnAction(event -> {
                    for (Node component : components) {
                        if (checkBox.isSelected() == (Boolean) current.dependentOnParameterValue()[0]) {
                            component.setDisable(false);
                        } else {
                            component.setDisable(true);
                        }
                    }
                });

            } else if (depends instanceof ComboBox) {
                final ComboBox comboBox = (ComboBox) depends;

                for (Node component : components) {
                    Object[] values = current.dependentOnParameterValue();
                    component.setDisable(true);
                    for (Object value : values) {
                        if (comboBox.getSelectionModel().getSelectedItem() == value) {
                            component.setDisable(false);
                            break;
                        }
                    }
                }

                comboBox.setOnAction(event -> {
                    for (Node component : components) {
                        Object[] values = current.dependentOnParameterValue();
                        component.setDisable(true);
                        for (Object value : values) {
                            if (comboBox.getSelectionModel().getSelectedItem() == value) {
                                component.setDisable(false);
                                break;
                            }
                        }
                    }
                });

            }
        }

    }

    private static final int DEFAULT_TOOL_TIP_LINE_LENGTH = 50;

    private String getToolTip(PluginParameter<?> parameter) {
        String description = parameter.description();
        int count = 0;
        StringBuilder builder = new StringBuilder();
        builder.append("<html>");
        for (int i = 0, n = description.length(); i < n; i++) {
            count++;
            if (description.charAt(i) == '\n') {
                builder.append("<br>");
                count = 0;
            } else if ((count > DEFAULT_TOOL_TIP_LINE_LENGTH) && (description.charAt(i) == ' ')) {
                builder.append("<br>");
                count = 0;
            } else {
                builder.append(description.charAt(i));
            }
        }
        builder.append("</html>");
        return builder.toString();
    }

    private HBox getLine(String label, TextField ref, Button button, String description) {

        HBox result = new HBox();
        result.setAlignment(Pos.CENTER);
        result.setPadding(new Insets(10.0));
        result.setSpacing(15.0);

        Label labelObj = new Label(label);
        labelObj.setTooltip(new Tooltip(description));
        result.getChildren().add(labelObj);

        ref.setEditable(true);
        //ref.setMaximumSize(ref.getPreferredSize());
        result.getChildren().add(ref);

        if (button != null) {
            result.getChildren().add(button);
        }

        return result;

    }

    private HBox getTaxaListPanel(String label, final TextField ref, String description, final Window parent, final TaxaList taxa) {

        HBox result = new HBox();
        //result.setToolTipText(description);

        result.getChildren().add(new Label(label));
        ref.setEditable(true);
        ref.setAlignment(Pos.CENTER_LEFT);
        result.getChildren().add(ref);

        if (taxa != null) {
            //final SelectFromAvailableDialog dialog = new SelectFromAvailableDialog(null, "Taxa Filter", new TaxaAvailableListModel(taxa));
//            Button taxaButton = new Button(new AbstractAction() {
//
//                @Override
//                public void actionPerformed(ActionEvent e) {
//                    dialog.setLocationRelativeTo(parent);
//                    dialog.setVisible(true);
//                    if (!dialog.isCanceled()) {
//                        int[] indicesToKeep = dialog.getDesiredIndices();
//                        StringBuilder builder = new StringBuilder();
//                        for (int i = 0; i < indicesToKeep.length; i++) {
//                            if (i != 0) {
//                                builder.append(",");
//                            }
//                            builder.append(taxa.taxaName(indicesToKeep[i]));
//                        }
//                        ref.setText(builder.toString());
//                    }
//                    dialog.setVisible(false);
//                }
//            });
            //taxaButton.setText("Select");
            //result.getChildren().add(taxaButton);
        }

        Button browse = getOpenFile(parent, ref);
        result.getChildren().add(browse);

        return result;

    }

//    private JPanel getPositionListPanel(String label, final JTextField ref, String description, final JDialog parent, final PositionList positions) {
//
//        JPanel result = new JPanel(new FlowLayout(FlowLayout.RIGHT));
//        result.setToolTipText(description);
//
//        result.add(new JLabel(label));
//        ref.setEditable(true);
//        ref.setHorizontalAlignment(JTextField.LEFT);
//        ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
//        ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
//        ref.setMaximumSize(ref.getPreferredSize());
//        result.add(ref);
//
//        if (positions != null) {
//            final SelectFromAvailableDialog dialog = new SelectFromAvailableDialog(null, "Site Name Filter", new SiteNamesAvailableListModel(positions));
//            JButton siteNamesButton = new JButton(new AbstractAction() {
//
//                @Override
//                public void actionPerformed(ActionEvent e) {
//                    dialog.setLocationRelativeTo(parent);
//                    dialog.setVisible(true);
//                    if (!dialog.isCanceled()) {
//                        int[] indicesToKeep = dialog.getDesiredIndices();
//                        StringBuilder builder = new StringBuilder();
//                        for (int i = 0; i < indicesToKeep.length; i++) {
//                            if (i != 0) {
//                                builder.append(",");
//                            }
//                            builder.append(positions.siteName(indicesToKeep[i]));
//                        }
//                        ref.setText(builder.toString());
//                    }
//                    dialog.setVisible(false);
//                }
//            });
//            siteNamesButton.setText("Select");
//            result.add(siteNamesButton);
//        }
//
//        return result;
//
//    }

    private Button getOpenFile(final Window parent, final TextField textField) {

        final FileChooser fileChooser = new FileChooser();
        fileChooser.setInitialDirectory(new File(TasselPrefs.getOpenDir()));

        Button result = new Button("Browse");

        result.setOnAction(event -> {
            File file = fileChooser.showOpenDialog(parent);
            if (file != null) {
                textField.setText(file.getPath());
                TasselPrefs.putOpenDir(Utils.getDirectory(file.getAbsolutePath()));
            }
        });

        return result;

    }

    private Button getSaveFile(final Window parent, final TextField textField) {

        final FileChooser fileChooser = new FileChooser();
        fileChooser.setInitialDirectory(new File(TasselPrefs.getSaveDir()));

        Button result = new Button("Browse");

        result.setOnAction(event -> {
            File file = fileChooser.showSaveDialog(parent);
            if (file != null) {
                textField.setText(file.getPath());
                TasselPrefs.putSaveDir(Utils.getDirectory(file.getAbsolutePath()));
            }
        });

        return result;
    }

    private Button getOpenDir(final Window parent, final TextField textField) {

        final DirectoryChooser fileChooser = new DirectoryChooser();
        fileChooser.setInitialDirectory(new File(Utils.getDirectory(TasselPrefs.getOpenDir())));

        Button result = new Button("Browse");

        result.setOnAction(event -> {
            File file = fileChooser.showDialog(parent);
            if (file != null) {
                textField.setText(file.getPath());
                TasselPrefs.putOpenDir(file.getPath());
            }
        });

        return result;
    }

    private Button getSaveDir(final Window parent, final TextField textField) {

        final DirectoryChooser fileChooser = new DirectoryChooser();
        fileChooser.setInitialDirectory(new File(Utils.getDirectory(TasselPrefs.getSaveDir())));

        Button result = new Button("Browse");

        result.setOnAction(event -> {
            File file = fileChooser.showDialog(parent);
            if (file != null) {
                textField.setText(file.getPath());
                TasselPrefs.putSaveDir(file.getPath());
            }
        });

        return result;
    }

    private class GenotypeWrapper {

        private final Object myObj;
        private final String myName;

        public GenotypeWrapper(Object obj, String name) {
            myObj = obj;
            myName = name;
        }

        @Override
        public String toString() {
            return myName;
        }

        public Object getObject() {
            return myObj;
        }

    }

    private Datum getGenotypeTable() {

        if (myCurrentInputData == null) {
            return null;
        }

        List<Datum> genotypeTables = myCurrentInputData.getDataOfType(GenotypeTable.class);
        if (!genotypeTables.isEmpty()) {
            return genotypeTables.get(0);
        }

        genotypeTables = myCurrentInputData.getDataOfType(GenotypePhenotype.class);
        if (!genotypeTables.isEmpty()) {
            return genotypeTables.get(0);
        }

        return null;
    }

    private TaxaList getTaxaList() {

        if (myCurrentInputData == null) {
            return null;
        }

        List<Datum> taxaList = myCurrentInputData.getDataOfType(GenotypeTable.class);
        if (!taxaList.isEmpty()) {
            return ((GenotypeTable) taxaList.get(0).getData()).taxa();
        }

        taxaList = myCurrentInputData.getDataOfType(TaxaList.class);
        if (!taxaList.isEmpty()) {
            return (TaxaList) taxaList.get(0).getData();
        }

        return null;
    }

    private List<Datum> getTaxaListDatum() {

        if (myCurrentInputData == null) {
            return null;
        }

        List<Datum> taxaList = myCurrentInputData.getDataOfType(TaxaList.class);
        if (!taxaList.isEmpty()) {
            return taxaList;
        }

        return null;

    }

    private PositionList getSiteNameList() {

        if (myCurrentInputData == null) {
            return null;
        }

        List<Datum> positionList = myCurrentInputData.getDataOfType(GenotypeTable.class);
        if (!positionList.isEmpty()) {
            return ((GenotypeTable) positionList.get(0).getData()).positions();
        }

        positionList = myCurrentInputData.getDataOfType(PositionList.class);
        if (!positionList.isEmpty()) {
            return (PositionList) positionList.get(0).getData();
        }

        return null;
    }

    private Datum getPositionList() {

        if (myCurrentInputData == null) {
            return null;
        }

        List<Datum> positionList = myCurrentInputData.getDataOfType(PositionList.class);
        if (!positionList.isEmpty()) {
            return positionList.get(0);
        }

        return null;
    }

    private List<Datum> getDistanceMatrices() {
        if (myCurrentInputData != null) {
            return myCurrentInputData.getDataOfType(DistanceMatrix.class);
        } else {
            return null;
        }
    }

    /**
     * Sets up this plugin to receive input from another plugin.
     *
     * @param input input
     */
    @Override
    public void receiveInput(Plugin input) {

        if (input == null) {
            throw new IllegalArgumentException("AbstractPlugin: receiveInput: input can not be null.");
        }

        if (!myInputs.contains(input)) {
            myInputs.add(input);
        }

        input.addListener(this);

    }

    /**
     * If interactive = true, the plugin will create dialogs and panels to interacts with the user
     *
     * @return boolean
     */
    @Override
    public boolean isInteractive() {
        return myIsInteractive;
    }

    /**
     * Adds listener to this plugin.
     *
     * @param listener listener to add
     */
    @Override
    public void addListener(PluginListener listener) {

        synchronized (myListeners) {
            if ((listener != null) && (!myListeners.contains(listener))) {
                myListeners.add(listener);
            }
        }

    }

    public List<PluginListener> getListeners() {
        return Collections.unmodifiableList(myListeners);
    }

    public boolean hasListeners() {
        return !myListeners.isEmpty();
    }

    public List<Plugin> getInputs() {
        return myInputs;
    }

    /**
     * Returns data set after complete.
     *
     * @param event event
     */
    protected void fireDataSetReturned(PluginEvent event) {

        synchronized (myListeners) {
            Iterator<PluginListener> itr = myListeners.iterator();
            while (itr.hasNext()) {
                try {
                    if (myThreaded) {
                        PluginListener current = itr.next();
                        ThreadedPluginListener thread = new ThreadedPluginListener(current, event);
                        thread.start();
                    } else {
                        PluginListener current = itr.next();
                        current.dataSetReturned(event);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

    }

    /**
     * Returns data set after complete.
     *
     * @param data data set
     */
    protected void fireDataSetReturned(DataSet data) {
        fireDataSetReturned(new PluginEvent(data));
    }

    private static final List<String> myPrintedCitations = new ArrayList<>();

    /**
     * Returns progress of execution.
     *
     * @param event event
     */
    protected void fireProgress(PluginEvent event) {

        synchronized (myListeners) {
            Iterator<PluginListener> itr = myListeners.iterator();
            while (itr.hasNext()) {
                PluginListener current = itr.next();
                current.progress(event);
            }
        }

        DataSet ds = (DataSet) event.getSource();
        if (ds != null) {
            List<Datum> percentage = ds.getDataOfType(Integer.class);

            if (percentage.size() > 0) {
                Datum datum = percentage.get(0);
                Integer percent = (Integer) datum.getData();
                if (percent == 100) {
                    if (!myPrintedCitations.contains(getCitation())) {
                        myLogger.info(getClass().getName() + "  Citation: " + getCitation());
                        myPrintedCitations.add(getCitation());
                    }
                }
            }
        }

    }

    /**
     * Returns progress of execution.
     *
     * @param percent percentage between 0 and 100 inclusive.
     */
    protected void fireProgress(Integer percent) {

        if ((percent < 0) || (percent > 100)) {
            throw new IllegalArgumentException("AbstractPlugin: fireProgress: percent must be between 0 and 100 inclusive.  arg: " + percent);
        }

        Datum percentage = new Datum("Percent", percent, null);
        fireProgress(new PluginEvent(new DataSet(percentage, this)));

    }

    @Override
    public String getCitation() {
        return DEFAULT_CITATION;
    }

    @Override
    public String pluginDescription() {
        return null;
    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual";
    }

    @Override
    public String icon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return this.getClass().getName().replace("Plugin", "");
    }

    @Override
    public String getToolTipText() {
        return getButtonName();
    }

    //
    // Methods for PluginListener.
    //

    /**
     * Returns data set after complete.
     *
     * @param event event
     */
    @Override
    public void dataSetReturned(PluginEvent event) {

        DataSet input = (DataSet) event.getSource();

        performFunction(input);

    }

    /**
     * No operation for this abstract class.
     */
    @Override
    public void progress(PluginEvent event) {
        // The default action of a plugin is to do
        // nothing when another plugin reports its
        // progress.   This is intended to be implemented
        // by GUI applications to show the user the
        // progress of an interactive action.
    }

    public void reverseTrace(int indent) {

        if (myTrace) {
            return;
        }

        indent(indent);
        System.out.println(getClass().getName());

        Iterator<Plugin> itr = myInputs.iterator();
        while (itr.hasNext()) {
            try {
                AbstractPlugin current = (AbstractPlugin) itr.next();
                current.reverseTrace(indent + 3);
            } catch (Exception e) {
                // do nothing
            }
        }

        myTrace = true;

    }

    public void trace(int indent) {

        if (myTrace) {
            return;
        }

        indent(indent);
        System.out.println(getClass().getName());

        Iterator<PluginListener> itr = myListeners.iterator();
        while (itr.hasNext()) {
            try {
                AbstractPlugin current = (AbstractPlugin) itr.next();
                current.trace(indent + 3);
            } catch (Exception e) {
                // do nothing
            }
        }

        myTrace = true;

    }

    private void indent(int indent) {

        for (int i = 0; i < indent; i++) {
            System.out.print(" ");
        }

    }

    @Override
    public void setThreaded(boolean threaded) {
        myThreaded = threaded;
    }

    @Override
    public boolean cancel() {
        return false;
    }

    @Override
    public void run() {
        performFunction(null);
    }

    @Override
    public void progress(int percent, Object meta) {
        fireProgress(percent);
    }

    @Override
    public boolean wasCancelled() {
        return myWasCancelled;
    }

}
