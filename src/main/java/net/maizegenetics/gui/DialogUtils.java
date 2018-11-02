/*
 * DialogUtils.java
 *
 * Created on April 17, 2006
 */
package net.maizegenetics.gui;

import java.awt.Component;
import javax.swing.JOptionPane;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class DialogUtils {

    private static final Logger myLogger = Logger.getLogger(DialogUtils.class);

    private static final int DEFAULT_MESSAGE_LINE_LENGTH = 50;

    private DialogUtils() {
    }

    public static void showWarning(String str, Component parent) {
        if (parent == null) {
            myLogger.warn(str);
        } else {
            JOptionPane.showMessageDialog(parent, getErrorMessage(str), "Warning", JOptionPane.WARNING_MESSAGE);
        }
    }

    public static void showError(String str, Component parent) {
        if (parent == null) {
            myLogger.error(str);
        } else {
            JOptionPane.showMessageDialog(parent, getErrorMessage(str), "Error", JOptionPane.ERROR_MESSAGE);
        }
    }

    public static void showError(Throwable e, Component parent) {
        String str = Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50, 7);
        if (parent == null) {
            myLogger.error(str);
        } else {
            JOptionPane.showMessageDialog(parent, getErrorMessage(str), "Error", JOptionPane.ERROR_MESSAGE);
        }
    }

    public static void showErrorCause(Throwable e, Component parent) {

        Throwable temp = e.getCause();

        if (temp != null) {
            showError(temp, parent);
        } else {
            showError(e, parent);
        }

    }

    private static String getErrorMessage(String message) {
        if (message.length() <= DEFAULT_MESSAGE_LINE_LENGTH) {
            return message;
        }
        int count = 0;
        StringBuilder builder = new StringBuilder();
        builder.append("<html>");
        for (int i = 0, n = message.length(); i < n; i++) {
            count++;
            if (message.charAt(i) == '\n') {
                builder.append("<br>");
                count = 0;
            } else if ((count > DEFAULT_MESSAGE_LINE_LENGTH) && (message.charAt(i) == ' ')) {
                builder.append("<br>");
                count = 0;
            } else {
                builder.append(message.charAt(i));
            }
        }
        builder.append("</html>");
        return builder.toString();
    }

}
