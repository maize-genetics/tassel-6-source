/*
 *  LoggingUtils
 */
package net.maizegenetics.util;

import net.maizegenetics.prefs.TasselPrefs;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.appender.ConsoleAppender;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.core.config.builder.api.*;
import org.apache.logging.log4j.core.config.builder.impl.BuiltConfiguration;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * @author Terry Casstevens
 */
public class LoggingUtils {

    private static final Logger myLogger = LogManager.getLogger(LoggingUtils.class);
    private static final PrintStream myOriginalOutputStream = System.out;
    private static final PrintStream myOriginalErrStream = System.err;

    private static PrintStream myPrintStream;

    private LoggingUtils() {
        // Utility Class
    }

    /**
     * Sends all logging messages (including log4j) to standard out. The logging level will be INFO, unless executed by
     * the GUI and the user preference has been set level to DEBUG.
     */
    public static void setupLogging() {
        if (TasselPrefs.getLogDebug()) {
            setupDebugLogging();
        } else {
            System.setOut(myOriginalOutputStream);
            System.setErr(myOriginalErrStream);
            sendLog4jToStdout();
        }
    }

    /**
     * Sends all logging messages (including log4j) to standard out. The logging level will be DEBUG.
     */
    public static void setupDebugLogging() {
        System.setOut(myOriginalOutputStream);
        System.setErr(myOriginalErrStream);
        sendDebugLog4jToStdout();
    }

    /**
     * Sends all logging messages (including log4j) to specified PrintStream. The logging level will be INFO, unless
     * executed by the GUI and the user preference has been set level to DEBUG. This is used to direct logging messages
     * to the GUI.
     *
     * @param stream stream
     */
    public static void setupLogging(PrintStream stream) {
        if (TasselPrefs.getLogDebug()) {
            setupDebugLogging(stream);
        } else {
            System.setOut(stream);
            System.setErr(stream);
            sendLog4jToStdout();
        }
    }

    /**
     * Sends all logging messages (including log4j) to specified PrintStream. The logging level will be DEBUG. This is
     * used to direct logging messages to the GUI.
     *
     * @param stream
     */
    public static void setupDebugLogging(PrintStream stream) {
        System.setOut(stream);
        System.setErr(stream);
        sendDebugLog4jToStdout();
    }

    /**
     * Sends all logging messages (including log4j) to specified file. The logging level will be INFO, unless executed
     * by the GUI and the user preference has been set level to DEBUG.
     *
     * @param logFileName file name
     *
     * @throws FileNotFoundException when file doesn't exist
     */
    public static void setupLogfile(String logFileName) throws FileNotFoundException {
        if (TasselPrefs.getLogDebug()) {
            setupDebugLogfile(logFileName);
        } else {
            File logFile = new File(logFileName);
            myLogger.info("Log File: " + logFile.getAbsolutePath());
            myPrintStream = new PrintStream(logFile);
            System.setOut(myPrintStream);
            System.setErr(myPrintStream);
            sendLog4jToStdout();
        }
    }

    /**
     * Sends all logging messages (including log4j) to specified file. The logging level will be DEBUG.
     *
     * @param logFileName file name
     *
     * @throws FileNotFoundException
     */
    public static void setupDebugLogfile(String logFileName) throws FileNotFoundException {
        File logFile = new File(logFileName);
        myLogger.info("Log File: " + logFile.getAbsolutePath());
        myPrintStream = new PrintStream(logFile);
        System.setOut(myPrintStream);
        System.setErr(myPrintStream);
        sendDebugLog4jToStdout();
    }

    private static void sendLog4jToStdout() {

        ConfigurationBuilder<BuiltConfiguration> builder = ConfigurationBuilderFactory.newConfigurationBuilder();
        builder.setStatusLevel(Level.ERROR);

        AppenderComponentBuilder console = builder.newAppender("Stdout", "Console")
                .addAttribute("target", ConsoleAppender.Target.SYSTEM_OUT);
        builder.add(console);

        console.add(builder.newLayout("PatternLayout")
                .addAttribute("pattern", "%r [%t] %p %c %notEmpty{%ndc }- %m%n"));

        RootLoggerComponentBuilder rootLogger = builder.newRootLogger(Level.INFO);
        rootLogger.add(builder.newAppenderRef("Stdout"));
        builder.add(rootLogger);

        LoggerComponentBuilder logger = builder.newLogger("net.maizegenetics", Level.INFO);
        logger.add(builder.newAppenderRef("Stdout"));
        logger.addAttribute("additivity", false);
        builder.add(logger);

        Configurator.reconfigure(builder.build());

    }

    private static void sendDebugLog4jToStdout() {

        ConfigurationBuilder<BuiltConfiguration> builder = ConfigurationBuilderFactory.newConfigurationBuilder();
        builder.setStatusLevel(Level.ERROR);

        AppenderComponentBuilder console = builder.newAppender("Stdout", "Console")
                .addAttribute("target", ConsoleAppender.Target.SYSTEM_OUT);
        builder.add(console);

        console.add(builder.newLayout("PatternLayout")
                .addAttribute("pattern", "%r [%t] %p %c %notEmpty{%ndc }- %m%n"));

        RootLoggerComponentBuilder rootLogger = builder.newRootLogger(Level.DEBUG);
        rootLogger.add(builder.newAppenderRef("Stdout"));
        builder.add(rootLogger);

        LoggerComponentBuilder logger = builder.newLogger("net.maizegenetics", Level.DEBUG);
        logger.add(builder.newAppenderRef("Stdout"));
        logger.addAttribute("additivity", false);
        builder.add(logger);

        Configurator.reconfigure(builder.build());

    }

}
