package net.maizegenetics.gui

import javafx.geometry.Insets
import javafx.geometry.Pos
import javafx.scene.Node
import javafx.scene.layout.StackPane

/**
 * @author Terry Casstevens
 * Created November 26, 2018
 */

class MainControls {

    val view = StackPane()

    init {
        view.alignment = Pos.CENTER_RIGHT
        view.padding = Insets(10.0)
    }

    fun clear() {
        view.children.clear()
    }

    fun add(vararg children: Node) {
        view.children.addAll(children)
    }

}