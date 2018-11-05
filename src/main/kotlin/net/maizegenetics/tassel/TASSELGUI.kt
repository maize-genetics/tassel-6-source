package net.maizegenetics.tassel

import javafx.application.Application
import javafx.application.Platform
import javafx.event.ActionEvent
import javafx.event.EventHandler
import javafx.scene.Scene
import javafx.scene.control.*
import javafx.scene.image.ImageView
import javafx.scene.input.KeyCombination
import javafx.scene.layout.*
import javafx.stage.Stage
import javafx.stage.WindowEvent
import net.maizegenetics.analysis.data.FileLoadPlugin
import net.maizegenetics.analysis.data.GenotypeSummaryPlugin
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.Plugin
import org.apache.log4j.Logger
import java.lang.Exception
import java.util.*


/**
 * @author Terry Casstevens
 * Created November 01, 2018
 */

private val myLogger = Logger.getLogger(TASSELGUI::class.java)

class TASSELGUI : Application() {

    lateinit var primaryStage: Stage
    private val myDataTree = DataTree()
    private val myMenuItemHash = HashMap<MenuItem, Plugin>()

    override fun start(stage: Stage) {

        instance = this

        primaryStage = stage

        stage.onCloseRequest = EventHandler<WindowEvent> { stop() }

        val info = TextArea()

        val left = VBox()
        left.children += myDataTree.dataTree
        left.children += info

        val progress = HBox()

        val root = BorderPane()

        val spacer = Region()
        spacer.styleClass.add("menu-bar")
        HBox.setHgrow(spacer, Priority.SOMETIMES)
        root.top = HBox(leftMenuBar(), spacer, rightMenuBar())
        updatePluginsWithGlobalConfigParameters()

        root.left = left
        root.bottom = progress

        stage.title = "TASSEL 6"
        val scene = Scene(root, 800.0, 600.0)
        scene.stylesheets += "/javafx/AppStyle.css"
        stage.scene = scene
        stage.show()

    }

    private fun leftMenuBar(): MenuBar {

        val result = MenuBar()

        result.menus += fileMenu()
        result.menus += dataMenu()

        return result

    }

    private fun rightMenuBar(): MenuBar {
        return MenuBar(helpMenu())
    }

    private fun fileMenu(): Menu {

        val result = Menu("File")

        val autoGuessPlugin = FileLoadPlugin(true, true)
        result.items += createMenuItem(autoGuessPlugin, "O", name = "Open", action = EventHandler {
            autoGuessPlugin.processData(null)
        })

        result.items += createMenuItem(FileLoadPlugin(true))

        result.items += createMenuItem(PreferencesDialog(true))

        result.items += SeparatorMenuItem()

        val exit = MenuItem("Exit")
        exit.onAction = EventHandler { stop() }
        result.items += exit

        return result

    }

    private fun dataMenu(): Menu {

        val result = Menu("Data")

        result.items += createMenuItem(GenotypeSummaryPlugin(true))

        return result
    }

    private fun helpMenu(): Menu {

        val result = Menu("Help")

        result.items += createMenuItem(TasselLogging.getInstance())

        return result

    }

    private fun createMenuItem(plugin: Plugin, mnemonic: String? = null, name: String? = null, action: EventHandler<ActionEvent>? = null): MenuItem {

        val menuItem = MenuItem(name ?: plugin.buttonName)
        try {
            menuItem.graphic = ImageView(plugin.icon())
        } catch (e: Exception) {
            myLogger.warn("Problem loading icon: ${plugin.icon()} for plugin: ${plugin.javaClass}")
        }

        if (mnemonic != null) {
            menuItem.accelerator = KeyCombination.keyCombination(mnemonic)
        }

        if (action == null) {
            menuItem.onAction = EventHandler {
                plugin.performFunction(myDataTree.selectedData())
                //val event = PluginEvent(myDataTree.selectedData())
                //val progressPanel = getProgressPanel()
                //progressPanel.addPlugin(plugin)
                //val thread = ThreadedPluginListener(plugin, event)
                //thread.start()
            }
        } else {
            menuItem.onAction = action
        }

        plugin.addListener(myDataTree)
        myMenuItemHash[menuItem] = plugin
        return menuItem

    }

    fun updatePluginsWithGlobalConfigParameters() {
        myMenuItemHash.values.forEach { (it as AbstractPlugin).setConfigParameters() }
    }

    override fun stop() {
        Platform.exit()
    }

    companion object {
        lateinit var instance: TASSELGUI
    }
}