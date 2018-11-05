@file:JvmName("AlertUtils")

package net.maizegenetics.gui

import javafx.scene.control.Alert
import javafx.scene.control.ButtonType

/**
 * @author Terry Casstevens
 * Created November 03, 2018
 */

fun showError(message: String) {
    val alert = Alert(Alert.AlertType.ERROR, message, ButtonType.OK)
    alert.showAndWait()
}

fun showWarn(message: String) {
    val alert = Alert(Alert.AlertType.WARNING, message, ButtonType.OK)
    alert.showAndWait()
}