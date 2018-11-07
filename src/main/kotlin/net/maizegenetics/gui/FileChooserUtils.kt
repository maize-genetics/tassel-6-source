@file:JvmName("FileChooserUtils")

package net.maizegenetics.gui

import javafx.application.Platform
import javafx.concurrent.Task
import javafx.stage.FileChooser
import net.maizegenetics.prefs.TasselPrefs
import net.maizegenetics.tassel.TASSELGUI
import java.io.File

/**
 * @author Terry Casstevens
 * Created November 05, 2018
 */

private val open by lazy {
    FileChooser().also { it.initialDirectory = File(TasselPrefs.getOpenDir()) }
}

fun singeFile(): File? {

    val task = object : Task<File?>() {
        override fun call(): File? {
            return open.showOpenDialog(TASSELGUI.instance.primaryStage)
        }
    }
    Platform.runLater(task)
    val result = task.get()

    if (result != null) {
        TasselPrefs.putOpenDir(result.absolutePath.substringBeforeLast('/'))
    }

    return result

}

fun multipleFiles(): List<File>? {

    val task = object : Task<List<File>?>() {
        override fun call(): List<File>? {
            return open.showOpenMultipleDialog(TASSELGUI.instance.primaryStage)
        }
    }
    Platform.runLater(task)
    val result = task.get()

    if (result != null) {
        TasselPrefs.putOpenDir(result[0].absolutePath.substringBeforeLast('/'))
    }

    return result

}

