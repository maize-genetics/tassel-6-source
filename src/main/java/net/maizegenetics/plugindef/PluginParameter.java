package net.maizegenetics.plugindef;

import com.google.common.collect.Range;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

import static net.maizegenetics.plugindef.AbstractPlugin.convert;

/**
 * Defines the attributes of parameters to be used in the plugins
 *
 * @param <T>
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public final class PluginParameter<T> {

    private static final Logger myLogger = LogManager.getLogger(PluginParameter.class);

    private final String myGuiName;
    private final String myUnits;
    private final String myCmdLineName;
    private final String myDescription;
    private final List<Range<Comparable<T>>> myRanges;
    private final T myDefaultValue;
    private final T myValue;
    private final boolean myRequired;
    private final Class<T> myClass;
    private final PluginParameter<?> myDependentOnParameter;
    private final Object[] myDependentOnParameterValue;
    private final List<T> myPossibleValues;
    private final boolean myIsNullable;

    public enum PARAMETER_TYPE {

        NA, IN_FILE, OUT_FILE, IN_DIR, OUT_DIR, GENOTYPE_TABLE, SITE_NAME_LIST,
        OBJECT_LIST_SINGLE_SELECT, POSITION_LIST, DISTANCE_MATRIX,
        LABEL, PASSWORD
    }

    private final PARAMETER_TYPE myParameterType;

    private PluginParameter(String guiName, String guiUnits, String cmdLineName,
                            String description, List<Range<Comparable<T>>> ranges, T defaultValue,
                            T value, boolean required, PARAMETER_TYPE fileType,
                            PluginParameter<?> dependentOnParameter, Object[] dependentOnParameterValue,
                            List<T> possibleValues, boolean isNullable, Class<T> type) {
        myGuiName = guiName;
        myUnits = guiUnits;
        myCmdLineName = cmdLineName;
        myDescription = description;
        myRanges = ranges;
        myDefaultValue = defaultValue;
        if (value == null) {
            myValue = defaultValue;
        } else {
            myValue = value;
        }

        if ((hasRange()) && (myValue != null)) {
            if (!acceptsValue(myValue)) {
                StringBuilder builder = new StringBuilder();
                builder.append("PluginParameter: init: " + myCmdLineName + " value: " + myValue.toString() + " outside range: ");
                builder.append(rangeToString());
                throw new IllegalArgumentException(builder.toString());
            }
        }

        myRequired = required;
        if ((myDefaultValue != null) && (myRequired)) {
            throw new IllegalArgumentException("PluginParameter: init: " + myCmdLineName + " shouldn't have default value and be required.");
        }
        myClass = type;
        myParameterType = fileType;
        myDependentOnParameter = dependentOnParameter;
        myDependentOnParameterValue = dependentOnParameterValue;
        if ((possibleValues != null) && (possibleValues.isEmpty())) {
            myPossibleValues = null;
        } else {
            myPossibleValues = possibleValues;
        }

        myIsNullable = isNullable;
    }

    /**
     * Use these to change the value of an existing parameter, e.g. after a user
     * changes the value. Otherwise use the Builder to create the parameter
     *
     * @param oldParameter
     * @param newValue
     */
    public PluginParameter(PluginParameter<T> oldParameter, T newValue) {
        this(oldParameter.myGuiName, oldParameter.myUnits, oldParameter.myCmdLineName,
                oldParameter.myDescription, oldParameter.myRanges, oldParameter.myDefaultValue, newValue,
                oldParameter.myRequired, oldParameter.myParameterType, oldParameter.dependentOnParameter(),
                oldParameter.dependentOnParameterValue(), oldParameter.possibleValues(), oldParameter.myIsNullable, oldParameter.myClass);
    }

    /**
     * Use this to change the possible values of a PluginParameter built as
     * objectListSingleSelect().
     *
     * @param oldParameter old plugin parameter
     * @param possibleValues new values
     */
    public PluginParameter(PluginParameter<T> oldParameter, List<T> possibleValues) {
        this(oldParameter.myGuiName, oldParameter.myUnits, oldParameter.myCmdLineName,
                oldParameter.myDescription, oldParameter.myRanges, oldParameter.myDefaultValue, oldParameter.value(),
                oldParameter.myRequired, oldParameter.myParameterType, oldParameter.dependentOnParameter(),
                oldParameter.dependentOnParameterValue(), possibleValues, oldParameter.myIsNullable, oldParameter.myClass);
    }

    public static PluginParameter<String> getLabelInstance(String label) {
        return new PluginParameter<>(label, null, label, "label", null, null, null, false, PARAMETER_TYPE.LABEL, null, null, null, false, String.class);
    }

    public String guiName() {
        return myGuiName;
    }

    public String units() {
        return myUnits;
    }

    public String cmdLineName() {
        return myCmdLineName;
    }

    public String description() {
        return myDescription;
    }

    public boolean hasRange() {
        if ((myRanges == null) || myRanges.isEmpty()) {
            return false;
        } else {
            return true;
        }
    }

    public String rangeToString() {
        return rangeToString(false);
    }

    public String rangeToString(boolean friendly) {
        if (hasRange()) {
            StringBuilder builder = new StringBuilder();

            if (myRanges.size() == 1) {
                Range<Comparable<T>> current = myRanges.get(0);
                if ((current.hasLowerBound() && current.hasUpperBound())
                        && (current.lowerEndpoint().equals(current.upperEndpoint()))) {
                    builder.append("[");
                    if (!friendly && (current.lowerEndpoint() instanceof Enum)) {
                        builder.append(((Enum) current.lowerEndpoint()).name());
                    } else {
                        builder.append(current.lowerEndpoint().toString());
                    }
                    builder.append("]");
                } else {
                    builder.append(current.toString());
                }
            } else {
                boolean first = true;
                builder.append("[");
                for (Range<Comparable<T>> current : myRanges) {
                    if (!first) {
                        builder.append(", ");
                    } else {
                        first = false;
                    }

                    if ((current.hasLowerBound() && current.hasUpperBound())
                            && (current.lowerEndpoint().equals(current.upperEndpoint()))) {
                        if (!friendly && (current.lowerEndpoint() instanceof Enum)) {
                            builder.append(((Enum) current.lowerEndpoint()).name());
                        } else {
                            builder.append(current.lowerEndpoint().toString());
                        }
                    } else {
                        builder.append(current.toString());
                    }
                }
                builder.append("]");
            }

            return builder.toString();
        } else {
            return "";
        }
    }

    public String possibleValuesString(boolean friendly) {
        if (hasPossibleValues()) {
            StringBuilder builder = new StringBuilder();

            if (myPossibleValues.size() == 1) {
                T current = myPossibleValues.get(0);
                builder.append("[");
                if (!friendly && (current instanceof Enum)) {
                    builder.append(((Enum) current).name());
                } else {
                    builder.append(current.toString());
                }
                builder.append("]");
            } else {
                boolean first = true;
                builder.append("[");
                for (T current : myPossibleValues) {
                    if (!first) {
                        builder.append(", ");
                    } else {
                        first = false;
                    }

                    if (!friendly && (current instanceof Enum)) {
                        builder.append(((Enum) current).name());
                    } else {
                        builder.append(current.toString());
                    }
                }
                builder.append("]");
            }

            return builder.toString();
        } else {
            return "";
        }
    }

    public boolean acceptsValue(Object value) {
        try {
            if (value == null && myIsNullable) {
                return true;
            }
            if (hasRange()) {
                for (Range<Comparable<T>> current : myRanges) {
                    if (current.contains((Comparable<T>) value)) {
                        return true;
                    }
                }
                return false;
            } else {
                return true;
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            return false;
        }
    }

    public boolean acceptsValue(String input) {
        Comparable<T> value = (Comparable<T>) convert(input, valueType());
        return acceptsValue(value);
    }

    public T value() {
        return myValue;
    }

    public T defaultValue() {
        return myDefaultValue;
    }

    public boolean required() {
        return myRequired;
    }

    public Class<T> valueType() {
        return myClass;
    }

    public PARAMETER_TYPE parameterType() {
        return myParameterType;
    }

    public PluginParameter<?> dependentOnParameter() {
        return myDependentOnParameter;
    }

    public Object[] dependentOnParameterValue() {
        return myDependentOnParameterValue;
    }

    public boolean hasPossibleValues() {
        return myPossibleValues != null;
    }

    public List<T> possibleValues() {
        return myPossibleValues;
    }

    public boolean isEmpty() {
        return (myValue == null || ((myValue instanceof String) && ((String) myValue).trim().isEmpty()));
    }

    public static class Builder<T> {

        private String myGuiName;
        private String myUnits = "";
        private final String myCmdLineName;
        private String myDescription = "";
        private List<Range<Comparable<T>>> myRanges = new ArrayList<>();
        private final T myDefaultValue;
        private boolean myIsRequired = false;
        private final Class<T> myClass;
        private PARAMETER_TYPE myParameterType = PARAMETER_TYPE.NA;
        private PluginParameter<?> myDependentOnParameter = null;
        private Object[] myDependentOnParameterValue = null;
        private List<T> myPossibleValues = null;
        private boolean myIsNullable = false;

        public Builder(String cmdLineName, T defaultValue, Class<T> type) {
            myCmdLineName = cmdLineName;
            myDefaultValue = defaultValue;
            myClass = type;
        }

        public Builder<T> units(String units) {
            myUnits = units;
            return this;
        }

        public Builder<T> description(String description) {
            myDescription = description;
            return this;
        }

        public Builder<T> range(Range<Comparable<T>> range) {
            if (range != null) {
                myRanges.add(range);
            }
            return this;
        }

        public Builder<T> range(Comparable<T>[] values) {
            if (values != null) {
                for (Comparable<T> current : values) {
                    myRanges.add(Range.singleton(current));
                }
            }
            return this;
        }

        public Builder<T> required(boolean required) {
            myIsRequired = required;
            return this;
        }

        public Builder<T> guiName(String guiName) {
            myGuiName = guiName;
            return this;
        }

        public Builder<T> inFile() {
            myParameterType = PARAMETER_TYPE.IN_FILE;
            return this;
        }

        public Builder<T> outFile() {
            myParameterType = PARAMETER_TYPE.OUT_FILE;
            return this;
        }

        public Builder<T> inDir() {
            myParameterType = PARAMETER_TYPE.IN_DIR;
            return this;
        }

        public Builder<T> outDir() {
            myParameterType = PARAMETER_TYPE.OUT_DIR;
            return this;
        }

        public Builder<T> genotypeTable() {
            myParameterType = PARAMETER_TYPE.GENOTYPE_TABLE;
            return this;
        }

        public Builder<T> distanceMatrix() {
            myParameterType = PARAMETER_TYPE.DISTANCE_MATRIX;
            return this;
        }

        public Builder<T> siteNameList() {
            myParameterType = PARAMETER_TYPE.SITE_NAME_LIST;
            return this;
        }

        public Builder<T> positionList() {
            myParameterType = PARAMETER_TYPE.POSITION_LIST;
            return this;
        }

        public Builder<T> objectListSingleSelect() {
            myParameterType = PARAMETER_TYPE.OBJECT_LIST_SINGLE_SELECT;
            return this;
        }

        public Builder<T> password() {
            myParameterType = PARAMETER_TYPE.PASSWORD;
            return this;
        }

        public Builder<T> dependentOnParameter(PluginParameter<?> parameter) {
            if (Boolean.class.isAssignableFrom(parameter.valueType())) {
                return dependentOnParameter(parameter, true);
            } else {
                throw new IllegalArgumentException("PluginParameter: dependentOnParameter: no default value for: " + parameter.valueType().getName());
            }
        }

        public Builder<T> dependentOnParameter(PluginParameter<?> parameter, Object value) {
            myDependentOnParameter = parameter;
            myDependentOnParameterValue = new Object[]{value};
            return this;
        }

        public Builder<T> dependentOnParameter(PluginParameter<?> parameter, Object[] values) {
            myDependentOnParameter = parameter;
            Object[] result = new Object[values.length];
            System.arraycopy(values, 0, result, 0, values.length);
            myDependentOnParameterValue = result;
            return this;
        }

        public Builder<T> possibleValues(List<T> possibleValues) {
            myPossibleValues = possibleValues;
            return this;
        }

        public Builder<T> nullable() {
            myIsNullable = true;
            return this;
        }

        public PluginParameter<T> build() {
            if ((myGuiName == null) || (myGuiName.isEmpty())) {
                StringBuilder builder = new StringBuilder();
                builder.append(Character.toUpperCase(myCmdLineName.charAt(0)));
                for (int i = 1; i < myCmdLineName.length(); i++) {
                    char current = myCmdLineName.charAt(i);
                    if (Character.isUpperCase(current)) {
                        builder.append(" ");
                    }
                    builder.append(current);
                }
                myGuiName = builder.toString();
            }
            if (myDescription.isEmpty()) {
                myDescription = myGuiName;
            }
            return new PluginParameter<>(myGuiName, myUnits, myCmdLineName,
                    myDescription, myRanges, myDefaultValue, null, myIsRequired,
                    myParameterType, myDependentOnParameter,
                    myDependentOnParameterValue, myPossibleValues, myIsNullable, myClass);
        }
    }
}
