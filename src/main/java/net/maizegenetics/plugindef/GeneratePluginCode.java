/*
 *  GeneratePluginCode
 */
package net.maizegenetics.plugindef;

import com.google.common.base.CaseFormat;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.lang.reflect.Field;

/**
 * @author Terry Casstevens
 */
public class GeneratePluginCode {

    private static final Logger myLogger = LogManager.getLogger(GeneratePluginCode.class);

    private GeneratePluginCode() {
    }

    public static void generate(Class currentMatch) {
        try {
            Plugin plugin = Plugin.getPluginInstance(currentMatch.getName());
            generate((AbstractPlugin) plugin);
        } catch (Exception ex) {
            myLogger.warn("Self-describing Plugins should implement this constructor: " + currentMatch.getClass().getName());
            myLogger.warn("public Plugin(boolean isInteractive) {");
            myLogger.warn("   super(isInteractive);");
            myLogger.warn("}");
            ex.printStackTrace();
        }
    }

    private static void generate(AbstractPlugin plugin) {
        String clazz = Utils.getBasename(plugin.getClass().getName());

        System.out.println("    // The following getters and setters were auto-generated.");
        System.out.println("    // Please use this method to re-generate.");
        System.out.println("    //");
        System.out.println("    // public static void main(String[] args) {");
        System.out.println("    //     GeneratePluginCode.generate(" + clazz + ".class);");
        System.out.println("    // }");
        System.out.println("");

        System.out.println("    /**");
        System.out.println("     * Convenience method to run plugin with one return object.");
        System.out.println("     */");
        System.out.println("    // TODO: Replace <Type> with specific type.");
        System.out.println("    public <Type> runPlugin(DataSet input) {");
        System.out.println("        return (<Type>) performFunction(input).getData(0).getData();");
        System.out.println("    }\n");

        for (Field field : plugin.getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(plugin);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            //String guiNameAsCamelCase = stringToCamelCase(current.guiName());
            String methodName = removeMyFromString(field.getName());

            //Getter
            System.out.println("    /**");
            System.out.println(createDescription(current.description()));
            System.out.println("     *");
            System.out.println("     * @return " + current.guiName());
            System.out.println("     */");
            System.out.println("    public " + current.valueType().getSimpleName() + " " + methodName + "() {");
            System.out.println("        return " + field.getName() + ".value();");
            System.out.println("    }\n");

            // Setter
            System.out.println("    /**");
            System.out.println(createDescription("Set " + current.guiName() + ". " + current.description()));
            System.out.println("     *");
            System.out.println("     * @param value " + current.guiName());
            System.out.println("     *");
            System.out.println("     * @return this plugin");
            System.out.println("     */");
            System.out.println("    public " + clazz + " " + methodName + "(" + current.valueType().getSimpleName() + " value) {");
            System.out.println("        " + field.getName() + " = new PluginParameter<>(" + field.getName() + ", value);");
            System.out.println("        return this;");
            System.out.println("    }\n");
        }
    }

    private static final int DEFAULT_DESCRIPTION_LINE_LENGTH = 50;

    private static String createDescription(String description) {
        int count = 0;
        StringBuilder builder = new StringBuilder();
        builder.append("     * ");
        for (int i = 0, n = description.length(); i < n; i++) {
            count++;
            if (description.charAt(i) == '\n') {
                builder.append("\n");
                builder.append("     * ");
                count = 0;
            } else if ((count > DEFAULT_DESCRIPTION_LINE_LENGTH) && (description.charAt(i) == ' ')) {
                builder.append("\n");
                builder.append("     * ");
                count = 0;
            } else {
                builder.append(description.charAt(i));
            }
        }
        return builder.toString();
    }

    private static String stringToCamelCase(String str) {
        StringBuilder builder = new StringBuilder();
        builder.append(Character.toLowerCase(str.charAt(0)));
        boolean makeUpper = false;
        for (int i = 1; i < str.length(); i++) {
            char current = str.charAt(i);
            if (current == ' ') {
                makeUpper = true;
            } else if (makeUpper) {
                builder.append(Character.toUpperCase(current));
                makeUpper = false;
            } else {
                builder.append(current);
            }
        }
        return builder.toString();
    }

    private static String removeMyFromString(String str) {
        String lower = str.toLowerCase();
        if (lower.startsWith("my")) {
            str = str.substring(2);
        }
        return CaseFormat.UPPER_CAMEL.to(CaseFormat.LOWER_CAMEL, str);
    }

}
