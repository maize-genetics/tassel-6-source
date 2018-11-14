package net.maizegenetics.gui

import javafx.application.Platform
import javafx.geometry.Insets
import javafx.geometry.Pos
import javafx.scene.control.Label
import javafx.scene.control.ProgressBar
import javafx.scene.layout.HBox
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Plugin
import net.maizegenetics.plugindef.PluginEvent
import net.maizegenetics.plugindef.PluginListener

/**
 * @author Terry Casstevens
 * Created November 05, 2018
 */

class ProgressViewer : PluginListener {

    val view = HBox(15.0)
    private val progress = ProgressBar(0.0)
    private lateinit var plugin: Plugin

    init {
        view.alignment = Pos.CENTER_LEFT
        view.padding = Insets(10.0)
        progress.prefWidth = 400.0
    }

    fun showProgress(plugin: Plugin) {
        this.plugin = plugin
        plugin.addListener(this)
        view.children.clear()
        progress.progress = 0.0
        view.children += Label(plugin.buttonName)
        view.children += progress
    }

    override fun dataSetReturned(event: PluginEvent?) {
        finished()
    }

    override fun progress(event: PluginEvent?) {
        val dataSet = event?.getSource() as DataSet
        val value = DataSet.data(dataSet, Int::class.javaObjectType).y
        progress.progress = value.toDouble() / 100.0
        if (progress.progress >= 1.0) finished()
    }

    private fun finished() {
        Platform.runLater {
            view.children.clear()
            view.children += Label("${plugin.buttonName} Finished")
        }
    }

}