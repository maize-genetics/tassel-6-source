package net.maizegenetics.gui

import javafx.scene.control.Label
import javafx.scene.control.ScrollPane
import javafx.scene.layout.GridPane
import net.maizegenetics.dna.factor.FactorTable
import net.maizegenetics.dna.factor.site.SNPSite
import java.lang.ref.WeakReference

/**
 * @author Terry Casstevens
 * Created November 16, 2018
 */

class FactorTableViewer private constructor(val table: FactorTable) {

    private val grid = GridPane()
    val view = ScrollPane(grid)

    init {
        grid.isGridLinesVisible = true
        addTaxa()
        addFactors()
        addValues()
    }


    private fun addTaxa() {

        table.taxa.forEachIndexed { index, taxon ->
            grid.add(Label(taxon.name), 0, index + 1)
        }
    }

    private fun addFactors() {

        table.factors.forEachIndexed { index, factor ->
            val label = Label("${factor.startChr}:${factor.startPos}")
            label.style = "-fx-rotate: -90;"
            grid.add(label, index + 1, 0)
        }

    }

    private fun addValues() {

        table.sites().forEach {
            val site = it.index
            when (it) {
                is SNPSite -> {
                    //for (taxon in 0 until it.taxa.numberOfTaxa()) {
                    for (taxon in 0 until 10) {
                        val cell = Label(it.genotypeAsString(taxon))
                        grid.add(cell, site + 1, taxon + 1)
                    }
                }
            }
        }

    }

    companion object {

        private val INSTANCES = HashMap<FactorTable, WeakReference<FactorTableViewer>>()

        @JvmStatic
        fun instance(table: FactorTable): FactorTableViewer {
            var result = INSTANCES[table]?.get()
            if (result == null) {
                result = FactorTableViewer(table)
                INSTANCES[table] = WeakReference(result)
            }
            return result
        }

    }

}