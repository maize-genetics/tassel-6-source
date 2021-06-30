package net.maizegenetics.gui

import javafx.application.Platform
import javafx.geometry.Insets
import javafx.geometry.Pos
import javafx.scene.Scene
import javafx.scene.control.Button
import javafx.scene.control.CheckBox
import javafx.scene.control.TextArea
import javafx.scene.control.Tooltip
import javafx.scene.layout.BorderPane
import javafx.scene.layout.HBox
import javafx.scene.layout.Priority
import javafx.scene.layout.VBox
import javafx.stage.FileChooser
import javafx.stage.Modality
import javafx.stage.Stage
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.prefs.TasselPrefs
import net.maizegenetics.tassel.TASSELMainApp
import net.maizegenetics.util.LoggingUtils
import net.maizegenetics.util.Utils
import org.apache.log4j.Logger
import java.io.IOException
import java.io.OutputStream
import java.io.PrintStream


/**
 * @author Terry Casstevens
 * Created November 07, 2018
 */

private val myLogger = Logger.getLogger(TasselLogging::class.java)

class TasselLogging private constructor() : AbstractPlugin() {

    private val dialog = Stage()
    private val text = TextArea()
    private val textAreaOutputStream = TextAreaStream(text)
    private val printStream = PrintStream(textAreaOutputStream, true)

    init {
        dialog.initModality(Modality.NONE)
        dialog.isResizable = true
        dialog.title = "TASSEL Logging"

        createDialog()
        LoggingUtils.setupLogging(printStream)
        basicLoggingInfo()

        if (TasselPrefs.getLogSendToConsole()) {
            LoggingUtils.setupLogging()
        }
    }

    override fun performFunction(input: DataSet?): DataSet? {
        LoggingUtils.setupLogging(printStream)
        Platform.runLater { dialog.show() }
        fireProgress(100)
        return null
    }

    private fun createDialog() {

        val main = BorderPane()
        val scene = Scene(main)

        val x = TasselPrefs.getLogXDim()
        val y = TasselPrefs.getLogYDim()
        if (x < 50 || y < 50) {
            dialog.setWidth(500.0)
            dialog.setHeight(400.0)
        } else {
            dialog.setWidth(x.toDouble())
            dialog.setHeight(y.toDouble())
        }

        text.isWrapText = true
        text.padding = Insets(10.0)
        text.isEditable = false
        HBox.setHgrow(text, Priority.ALWAYS)
        VBox.setVgrow(text, Priority.ALWAYS)

        val isDebug = CheckBox("Debug Level")
        isDebug.isSelected = TasselPrefs.getLogDebug()
        isDebug.tooltip = Tooltip("Set to show Debug Logging Messages")
        isDebug.setOnAction { event ->
            val debugMode = isDebug.isSelected
            isDebug.isSelected = debugMode
            TasselPrefs.putLogDebug(debugMode)
            LoggingUtils.setupLogging(printStream)
        }

        val closeButton = Button("Close")
        closeButton.setOnAction { close() }

        val clearButton = Button("Clear")
        clearButton.setOnAction {
            text.clear()
            basicLoggingInfo()
        }

        val saveButton = Button("Save")
        saveButton.setOnAction {
            val chooser = FileChooser()
            val theFile = chooser.showOpenDialog(dialog)
            if (theFile != null) {
                try {
                    Utils.getBufferedWriter(theFile).use { writer -> writer.write(text.getText()) }
                } catch (ex: Exception) {
                    showError(ex.message)
                }

            }
        }

        val pnlButtons = HBox()
        pnlButtons.alignment = Pos.CENTER
        pnlButtons.padding = Insets(10.0)
        pnlButtons.spacing = 20.0
        pnlButtons.children.add(closeButton)
        pnlButtons.children.add(clearButton)
        pnlButtons.children.add(saveButton)
        pnlButtons.children.add(isDebug)

        main.center = text
        main.bottom = pnlButtons

        dialog.setScene(scene)

    }

    private fun updateLogging() {
        if (TasselPrefs.getLogSendToConsole()) {
            LoggingUtils.setupLogging()
        } else {
            LoggingUtils.setupDebugLogging(printStream)
        }
    }

    private fun close() {
        TasselPrefs.putLogXDim(dialog.width.toInt())
        TasselPrefs.putLogYDim(dialog.height.toInt())
        dialog.close()
        if (TasselPrefs.getLogSendToConsole()) {
            LoggingUtils.setupLogging()
        }
    }

    override fun icon(): String {
        return "/images/log.gif"
    }

    override fun getButtonName(): String {
        return "Logging"
    }

    override fun getToolTipText(): String {
        return "Logging"
    }

    companion object {
        val instance = TasselLogging()

        fun updateLoggingLocation() {
            instance.updateLogging()
        }
    }
}

fun basicLoggingInfo() {
    myLogger.info("Tassel Version: " + TASSELMainApp.version + "  Date: " + TASSELMainApp.versionDate)
    myLogger.info("Max Available Memory Reported by JVM: " + Utils.getMaxHeapSizeMB() + " MB")
    myLogger.info("Java Version: " + System.getProperty("java.version"))
    myLogger.info("OS: " + System.getProperty("os.name"))
    myLogger.info("Number of Processors: " + Runtime.getRuntime().availableProcessors())
}

class TextAreaStream(private val textArea: TextArea) : OutputStream() {

    override fun write(i: Int) {
        Platform.runLater { textArea.appendText(i.toChar().toString()) }
    }

}