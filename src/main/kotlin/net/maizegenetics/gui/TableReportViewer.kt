package net.maizegenetics.gui

import javafx.beans.property.ReadOnlyObjectWrapper
import javafx.collections.FXCollections
import javafx.scene.control.TableColumn
import javafx.scene.control.TableView
import net.maizegenetics.util.TableReport
import java.util.*


/**
 * @author Terry Casstevens
 * Created November 05, 2018
 */

class TableReportViewer(report: TableReport) {

    val view = TableView<Array<Any>>()

    init {

        view.isEditable = false

        val columnNames = report.tableColumnNames
        for (i in columnNames.indices) {
            val column = TableColumn<Array<Any>, Any>(columnNames[i].toString())
            column.setCellValueFactory { ReadOnlyObjectWrapper<Any>(it.value[i]) }
            view.columns.add(column)
        }

        view.items = FXCollections.observableArrayList<Array<Any>>(TableReportList(report, report.rowCount.toInt()))

    }

}

class TableReportList(private val report: TableReport, override val size: Int) : AbstractList<Array<Any>>() {

    override fun get(index: Int): Array<Any> {
        return report.getRow(index.toLong())
    }

}