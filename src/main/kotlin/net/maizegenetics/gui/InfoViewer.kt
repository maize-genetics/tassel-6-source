package net.maizegenetics.gui

import javafx.scene.control.ScrollPane
import javafx.scene.control.TextArea
import javafx.scene.layout.HBox
import javafx.scene.layout.Priority
import javafx.scene.layout.VBox
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.util.TableReport

/**
 * @author Terry Casstevens
 * Created November 18, 2018
 */

class InfoViewer {

    private val text = TextArea()
    val view = ScrollPane(text)

    init {
        text.isEditable = false
        HBox.setHgrow(view, Priority.ALWAYS)
        VBox.setVgrow(view, Priority.ALWAYS)
        view.isFitToHeight = true;
        view.isFitToWidth = true;
    }

    fun show(datum: Datum) {

        val builder = StringBuilder()
        when (val data = datum.data) {

            is TableReport -> {
                val numColumns = data.columnCount.toLong()
                val numRows = data.rowCount
                builder.append("Table Title: ")
                builder.append(data.tableTitle)
                builder.append("\n")
                builder.append("Number of columns: ")
                builder.append(numColumns)
                builder.append("\n")
                builder.append("Number of rows: ")
                builder.append(numRows)
                builder.append("\n")
                builder.append("Matrix size (excludes row headers): ")
                builder.append((numColumns - 1) * numRows)
                builder.append("\n")
            }

        }

        text.clear()
        text.appendText(builder.toString())

    }

}