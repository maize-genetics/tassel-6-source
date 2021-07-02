package net.maizegenetics.tassel

import javafx.scene.control.ScrollPane
import javafx.scene.control.SelectionMode
import javafx.scene.control.TreeItem
import javafx.scene.control.TreeView
import net.maizegenetics.analysis.data.GenotypeSummaryPlugin
import net.maizegenetics.dna.map.PositionList
import net.maizegenetics.dna.snp.FilterList
import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.phenotype.Phenotype
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.plugindef.PluginEvent
import net.maizegenetics.plugindef.PluginListener
import net.maizegenetics.taxa.IdentifierSynonymizer
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.taxa.distance.DistanceMatrix
import net.maizegenetics.taxa.tree.Tree
import net.maizegenetics.util.TableReport
import org.apache.logging.log4j.LogManager
import java.util.ArrayList
import kotlin.collections.HashMap

/**
 * @author Terry Casstevens
 * Created November 04, 2018
 */

const val NODE_TYPE_DATA = "Data"
const val NODE_TYPE_RESULT = "Result"
const val NODE_TYPE_SEQUENCE = "Genotype Table"
const val NODE_TYPE_POLYMORPHISMS = "Polymorphisms"
const val NODE_TYPE_NUMERICAL = "Numerical"
const val NODE_TYPE_MATRIX = "Matrix"
const val NODE_TYPE_TREE = "Tree"
const val NODE_TYPE_LISTS = "Lists"
const val NODE_TYPE_FILTERS = "Filters"
const val NODE_TYPE_FUSIONS = "Fusions"
const val NODE_TYPE_SYNONYMS = "Synonyms"
const val NODE_TYPE_DIVERSITY = "Diversity"
const val NODE_TYPE_SNP_ASSAYS = "SNP Assays"
const val NODE_TYPE_LD = "LD"
const val NODE_TYPE_ASSOCIATIONS = "Association"
const val NODE_TYPE_VARIANCES = "Variances"
const val NODE_TYPE_SYNONYMIZER = "Synonymizer"
const val NODE_TYPE_STEPWISE = "Stepwise"
const val NODE_TYPE_GENO_SUMMARY = "Genotype Summary"
const val NODE_TYPE_DEFAULT = NODE_TYPE_DATA

private val myLogger = LogManager.getLogger(DataTree::class.java)

class DataTree : PluginListener {

    private val treeRoot = TreeItem<Any>()
    private val dataTree = TreeView<Any>(treeRoot)
    val view = ScrollPane(dataTree)
    private val treeItems = HashMap<Any, TreeItem<Any>>()

    init {

        dataTree.isShowRoot = false
        dataTree.selectionModel.selectionMode = SelectionMode.MULTIPLE
        treeRoot.isExpanded = true

        dataTree.selectionModel.selectedItemProperty().addListener { observable, oldValue, newValue ->
            val value = newValue?.value
            if (value is Datum) {
                TASSELGUI.instance.changeView(value)
            }
        }

        view.isFitToHeight = true;
        view.isFitToWidth = true;

    }

    fun add(data: DataSet) {

        val creator = data.creator
        for (i in 0 until data.size) {

            val d = data.getData(i)

            val child = when (creator) {
                is GenotypeSummaryPlugin -> addDatum(NODE_TYPE_GENO_SUMMARY, d)
                // is MLMPlugin, is FixedEffectLMPlugin, is EqtlAssociationPlugin -> addDatum(NODE_TYPE_ASSOCIATIONS, d)
                // is SequenceDiversityPlugin -> addDatum(NODE_TYPE_DIVERSITY, d)
                // is LinkageDisequilibriumPlugin -> addDatum(NODE_TYPE_LD, d)
                // is GenotypeTableMask -> addDatum(d)

                else -> when (d.data) {
                    is GenotypeTable -> addDatum(NODE_TYPE_SEQUENCE, d)
                    is IdentifierSynonymizer -> addDatum(NODE_TYPE_SYNONYMIZER, d)
                    is Phenotype -> addDatum(NODE_TYPE_NUMERICAL, d)
                    is DistanceMatrix -> addDatum(NODE_TYPE_MATRIX, d)
                    is TableReport -> addDatum(NODE_TYPE_NUMERICAL, d)
                    is FilterList -> addDatum(NODE_TYPE_FILTERS, d)
                    is Tree -> addDatum(NODE_TYPE_TREE, d)
                    is TaxaList -> addDatum(NODE_TYPE_LISTS, d)
                    is PositionList -> addDatum(NODE_TYPE_LISTS, d)
                    else -> addDatum(NODE_TYPE_DEFAULT, d)
                }

            }

            if (i == 0) {
                dataTree.selectionModel.clearSelection()
                dataTree.selectionModel.select(child)
            }

        }

    }

    private fun addDatum(parent: String, datum: Datum): TreeItem<Any> {
        val node = treeItem(parent)
        val child = TreeItem<Any>(datum)
        child.isExpanded = true
        node.children += child
        return child
    }

    fun treeItem(parent: String): TreeItem<Any> {
        var result = treeItems[parent]
        if (result != null) return result
        result = TreeItem(parent)
        treeItems[parent] = result
        treeRoot.children += result
        result.isExpanded = true
        return result
    }

    fun selectedData(): DataSet {
        val data = ArrayList<Datum>()
        dataTree.selectionModel.selectedItems
                .map { it.value }
                .filter { it is Datum }
                .forEach { data.add(it as Datum) }
        return DataSet(data, null)
    }

    fun deleteSelectedNodes() {
        val items = dataTree.selectionModel.selectedItems
        items.filter { it.value is Datum }
                .forEach { it.parent.children.remove(it) }
    }

    override fun dataSetReturned(event: PluginEvent?) {
        val tds = event?.source as DataSet
        add(tds)
    }

    override fun progress(event: PluginEvent?) {
        // Do nothing
    }

}