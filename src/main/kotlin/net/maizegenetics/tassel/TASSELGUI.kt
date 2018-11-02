package net.maizegenetics.tassel

import javafx.application.Application
import javafx.application.Platform
import javafx.event.ActionEvent
import javafx.event.EventHandler
import javafx.scene.Scene
import javafx.scene.control.*
import javafx.scene.layout.BorderPane
import javafx.scene.layout.HBox
import javafx.scene.layout.VBox
import javafx.stage.Stage
import javafx.stage.WindowEvent

/**
 * @author Terry Casstevens
 * Created November 01, 2018
 */


class TASSELGUI : Application() {

    override fun start(stage: Stage) {

        stage.onCloseRequest = EventHandler<WindowEvent> { stop() }

        val menuBar = MenuBar()
        menuBar.menus += fileMenu()

        val dataTree = TreeView<Any>()

        val info = TextArea()

        val left = VBox()
        left.children += dataTree
        left.children += info

        val progress = HBox()

        val root = BorderPane()
        root.top = menuBar
        root.left = left
        root.bottom = progress

        stage.title = "TASSEL 6"
        stage.scene = Scene(root, 800.0, 600.0)
        stage.show()

    }

    private fun fileMenu(): Menu {

        val result = Menu("File")

        val exit = MenuItem("Exit")
        exit.onAction = EventHandler { stop() }
        result.items += exit

        return result

    }

    override fun stop() {
        Platform.exit()
    }
}