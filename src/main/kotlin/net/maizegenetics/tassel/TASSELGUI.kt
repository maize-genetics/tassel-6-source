package net.maizegenetics.tassel

import javafx.application.Application
import javafx.application.Platform
import javafx.event.ActionEvent
import javafx.event.EventHandler
import javafx.geometry.Orientation
import javafx.scene.Scene
import javafx.scene.control.*
import javafx.scene.image.ImageView
import javafx.scene.input.KeyCombination
import javafx.scene.layout.*
import javafx.scene.paint.Paint
import javafx.stage.Stage
import javafx.stage.WindowEvent
import net.maizegenetics.analysis.data.FileLoadPlugin
import net.maizegenetics.analysis.data.GenotypeSummaryPlugin
import net.maizegenetics.dna.factor.FeatureTable
import net.maizegenetics.gui.*
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.plugindef.Plugin
import net.maizegenetics.util.TableReport
import org.apache.log4j.Logger
import java.util.*


/**
 * @author Terry Casstevens
 * Created November 01, 2018
 */

private val myLogger = Logger.getLogger(TASSELGUI::class.java)

private val EMPTY_NODE = StackPane().also {
    it.background = Background(BackgroundFill(Paint.valueOf("grey"), null, null))
    it.children += ImageView("/images/Tassel_Logo.png")
}

class TASSELGUI : Application() {

    lateinit var primaryStage: Stage
    private val myDataTree = DataTree()
    private val myInfo = InfoViewer()
    private val myMenuItemHash = HashMap<MenuItem, Plugin>()
    private val myMainPane = BorderPane()
    private val myProgressViewer = ProgressViewer()
    private val myMainControls = MainControls()
    private val myMainView = StackPane()

    override fun start(stage: Stage) {

        instance = this

        primaryStage = stage

        stage.onCloseRequest = EventHandler<WindowEvent> { stop() }

        val left = SplitPane()
        left.orientation = Orientation.VERTICAL
        left.setDividerPositions(0.6)
        left.items += myDataTree.view
        left.items += myInfo.view
        left.prefWidth = 300.0

        val spacer = Region()
        spacer.styleClass.add("menu-bar")
        HBox.setHgrow(spacer, Priority.SOMETIMES)
        myMainPane.top = HBox(leftMenuBar(), spacer, rightMenuBar())
        updatePluginsWithGlobalConfigParameters()

        myMainPane.center = myMainView
        myMainView.children += EMPTY_NODE

        myMainPane.left = left

        val bottomSpacer = Region()
        HBox.setHgrow(bottomSpacer, Priority.ALWAYS)
        val bottom = HBox(myProgressViewer.view, bottomSpacer, myMainControls.view)
        myMainPane.bottom = bottom

        stage.title = "TASSEL 6"
        val scene = Scene(myMainPane, 1000.0, 750.0)
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
        result.items += createMenuItem(autoGuessPlugin, name = "Open", action = EventHandler {
            myProgressViewer.showProgress(autoGuessPlugin)
            Thread { autoGuessPlugin.processData(null) }.start()
        })

        result.items += createMenuItem(FileLoadPlugin(true))

        result.items += MenuItem("Delete", ImageView("/images/trash.gif")).also {
            it.onAction = EventHandler { myDataTree.deleteSelectedNodes() }
        }

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

        result.items += createMenuItem(TasselLogging.instance)

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
                myProgressViewer.showProgress(plugin)
                Thread(Runnable { plugin.performFunction(myDataTree.selectedData()) }).start()
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

    fun changeView(datum: Datum) = Platform.runLater {

        myProgressViewer.showSelected(datum)

        myInfo.show(datum)

        myMainControls.clear()

        myMainView.children.clear()

        when (val data = datum.data) {
            is TableReport -> myMainView.children += TableReportViewer.instance(data).view
            is FeatureTable -> {
                val viewer = FactorTableViewer.instance(data, myMainView.widthProperty())
                myMainView.children += viewer.view
                myMainControls.add(viewer.controls)
            }
            else -> myMainView.children += EMPTY_NODE
        }

    }

    override fun stop() {
        Platform.exit()
    }

    companion object {
        @JvmStatic
        lateinit var instance: TASSELGUI
    }
}