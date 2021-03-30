package net.maizegenetics.gui

import javafx.beans.property.ReadOnlyObjectWrapper
import javafx.collections.ObservableListBase
import javafx.collections.transformation.SortedList
import javafx.scene.control.TableCell
import javafx.scene.control.TableColumn
import javafx.scene.control.TableView
import net.maizegenetics.dna.map.Chromosome
import net.maizegenetics.util.TableReport
import java.lang.ref.WeakReference
import kotlin.collections.set


/**
 * @author Terry Casstevens
 * Created November 05, 2018
 */

class TableReportViewer private constructor(report: TableReport) {

    val view = TableView<Array<Any>>()

    init {

        view.isEditable = false

        val firstRow = report.getRow(0)

        val columnNames = report.tableColumnNames
        for (i in columnNames.indices) {

            when (firstRow[i]) {

                is String -> {
                    when (columnNames[i].toString().toLowerCase()) {
                        "chromosome", "chr" -> {
                            val column = TableColumn<Array<Any>, Chromosome>(columnNames[i].toString())
                            column.setCellValueFactory { ReadOnlyObjectWrapper<Chromosome>(if (it.value[i] == null) null else Chromosome.instance(it.value[i].toString())) }
                            view.columns.add(column)
                        }
                        "taxa", "taxa name", "taxon", "taxon name" -> {
                            val column = TableColumn<Array<Any>, String>(columnNames[i].toString())
                            column.setCellValueFactory { ReadOnlyObjectWrapper<String>(if (it.value[i] == null) null else it.value[i].toString()) }
                            column.setCellFactory {
                                val result = object : TableCell<Array<Any>, String>() {
                                    override fun updateItem(item: String?, empty: Boolean) {
                                        super.updateItem(item, empty)
                                        text = item
                                        if (item != null && item.matches(Regex("""^[AaBbCcDdEeFf].*"""))) {
                                            style = "-fx-background-color: yellow"
                                        } else {
                                            style = "-fx-background-color: #ab4642"
                                        }
                                    }
                                }
                                result
                            }
                            view.columns.add(column)
                        }
                        else -> {
                            val column = TableColumn<Array<Any>, String>(columnNames[i].toString())
                            column.setCellValueFactory { ReadOnlyObjectWrapper<String>(if (it.value[i] == null) null else it.value[i].toString()) }
                            view.columns.add(column)
                        }
                    }
                }
                is Int -> {
                    val column = TableColumn<Array<Any>, Int>(columnNames[i].toString())
                    column.setCellValueFactory { ReadOnlyObjectWrapper<Int>(if (it.value[i] == null) null else it.value[i] as Int) }
                    view.columns.add(column)
                }
                is Double -> {
                    val column = TableColumn<Array<Any>, Double>(columnNames[i].toString())
                    column.setCellValueFactory { ReadOnlyObjectWrapper<Double>(if (it.value[i] == null) null else it.value[i] as Double) }
                    view.columns.add(column)
                }
                is Float -> {
                    val column = TableColumn<Array<Any>, Float>(columnNames[i].toString())
                    column.setCellValueFactory { ReadOnlyObjectWrapper<Float>(if (it.value[i] == null) null else it.value[i] as Float) }
                    view.columns.add(column)
                }
                is Chromosome -> {
                    val column = TableColumn<Array<Any>, Chromosome>(columnNames[i].toString())
                    column.setCellValueFactory { ReadOnlyObjectWrapper<Chromosome>(if (it.value[i] == null) null else it.value[i] as Chromosome) }
                    view.columns.add(column)
                }
                else -> {
                    val column = TableColumn<Array<Any>, Any>(columnNames[i].toString())
                    column.setCellValueFactory { ReadOnlyObjectWrapper<Any>(it.value[i]) }
                    view.columns.add(column)
                }

            }

        }

        val sorted = SortedList<Array<Any>>(TableReportList(report))
        view.items = sorted
        sorted.comparatorProperty().bind(view.comparatorProperty());

    }

    companion object {

        private val INSTANCES = HashMap<TableReport, WeakReference<TableReportViewer>>()

        @JvmStatic
        fun instance(report: TableReport): TableReportViewer {
            var result = INSTANCES[report]?.get()
            if (result == null) {
                result = TableReportViewer(report)
                INSTANCES[report] = WeakReference(result)
            }
            return result
        }

    }

}

private class TableReportList(val report: TableReport) : ObservableListBase<Array<Any>>() {

    private val CACHE = HashMap<Int, WeakReference<Array<Any>>>()

    override fun get(index: Int): Array<Any> {
        var result = CACHE[index]?.get()
        if (result == null) {
            result = report.getRow(index.toLong())
            CACHE[index] = WeakReference(result)
        }
        return result!!
    }

    override val size: Int
        get() = report.rowCount.toInt()

}