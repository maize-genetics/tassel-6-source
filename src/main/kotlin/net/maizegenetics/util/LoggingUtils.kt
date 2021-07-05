@file:JvmName("LoggingUtils")

package net.maizegenetics.util

import net.maizegenetics.prefs.TasselPrefs
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import org.apache.logging.log4j.core.appender.ConsoleAppender
import org.apache.logging.log4j.core.config.Configurator
import org.apache.logging.log4j.core.config.builder.api.ConfigurationBuilderFactory
import org.apache.logging.log4j.core.layout.PatternLayout
import java.io.File
import java.io.FileNotFoundException
import java.io.PrintStream

private val myLogger = LogManager.getLogger("net.maizegenetics.util.LoggingUtils")

private val myOriginalOutputStream = System.out
private val myOriginalErrStream = System.err
private var myPrintStream: PrintStream? = null

/**
 * Sends all logging messages (including log4j) to standard out. The logging level will be INFO, unless executed by
 * the GUI and the user preference has been set level to DEBUG.
 */
fun setupLogging() {
    if (TasselPrefs.getLogDebug()) {
        setupDebugLogging()
    } else {
        System.setOut(myOriginalOutputStream)
        System.setErr(myOriginalErrStream)
        sendLog4jToStdout()
    }
}

/**
 * Sends all logging messages (including log4j) to standard out. The logging level will be DEBUG.
 */
fun setupDebugLogging() {
    System.setOut(myOriginalOutputStream)
    System.setErr(myOriginalErrStream)
    sendDebugLog4jToStdout()
}

/**
 * Sends all logging messages (including log4j) to specified PrintStream. The logging level will be INFO, unless
 * executed by the GUI and the user preference has been set level to DEBUG. This is used to direct logging messages
 * to the GUI.
 *
 * @param stream stream
 */
fun setupLogging(stream: PrintStream?) {
    if (TasselPrefs.getLogDebug()) {
        setupDebugLogging(stream)
    } else {
        System.setOut(stream)
        System.setErr(stream)
        sendLog4jToStdout()
    }
}

/**
 * Sends all logging messages (including log4j) to specified PrintStream. The logging level will be DEBUG. This is
 * used to direct logging messages to the GUI.
 *
 * @param stream
 */
fun setupDebugLogging(stream: PrintStream?) {
    System.setOut(stream)
    System.setErr(stream)
    sendDebugLog4jToStdout()
}

/**
 * Sends all logging messages (including log4j) to specified file. The logging level will be INFO, unless executed
 * by the GUI and the user preference has been set level to DEBUG.
 *
 * @param logFileName file name
 *
 * @throws FileNotFoundException when file doesn't exist
 */
@Throws(FileNotFoundException::class)
fun setupLogfile(logFileName: String?) {
    if (TasselPrefs.getLogDebug()) {
        setupDebugLogfile(logFileName)
    } else {
        val logFile = File(logFileName)
        myLogger.info("Log File: " + logFile.absolutePath)
        myPrintStream = PrintStream(logFile)
        System.setOut(myPrintStream)
        System.setErr(myPrintStream)
        sendLog4jToStdout()
    }
}

/**
 * Sends all logging messages (including log4j) to specified file. The logging level will be DEBUG.
 *
 * @param logFileName file name
 *
 * @throws FileNotFoundException
 */
@Throws(FileNotFoundException::class)
fun setupDebugLogfile(logFileName: String?) {
    val logFile = File(logFileName)
    myLogger.info("Log File: " + logFile.absolutePath)
    myPrintStream = PrintStream(logFile)
    System.setOut(myPrintStream)
    System.setErr(myPrintStream)
    sendDebugLog4jToStdout()
}

private fun sendLog4jToStdout() {
    val builder = ConfigurationBuilderFactory.newConfigurationBuilder()
    builder.setStatusLevel(Level.ERROR)

    val console = builder.newAppender("Stdout", "Console")
    builder.add(console)

    val standard = builder.newLayout("PatternLayout")
    standard.addAttribute("pattern", PatternLayout.TTCC_CONVERSION_PATTERN)
    console.add(standard)

    val logger = builder.newLogger("net.maizegenetics", Level.INFO)
    logger.add(builder.newAppenderRef("Stdout"))
    logger.addAttribute("additivity", false)
    builder.add(logger)

    Configurator.reconfigure(builder.build())
}

private fun sendDebugLog4jToStdout() {
    val builder = ConfigurationBuilderFactory.newConfigurationBuilder()
    builder.setStatusLevel(Level.ERROR)

    val console = builder.newAppender("Stdout", "Console")
    builder.add(console)

    val standard = builder.newLayout("PatternLayout")
    standard.addAttribute("pattern", PatternLayout.TTCC_CONVERSION_PATTERN)
    console.add(standard)

    val logger = builder.newLogger("net.maizegenetics", Level.DEBUG)
    logger.add(builder.newAppenderRef("Stdout"))
    logger.addAttribute("additivity", false)
    builder.add(logger)

    Configurator.reconfigure(builder.build())
}