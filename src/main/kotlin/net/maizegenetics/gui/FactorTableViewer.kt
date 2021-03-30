package net.maizegenetics.gui

import javafx.beans.property.DoubleProperty
import javafx.beans.property.ReadOnlyDoubleProperty
import javafx.beans.property.SimpleDoubleProperty
import javafx.event.EventHandler
import javafx.geometry.Insets
import javafx.geometry.Pos
import javafx.scene.Group
import javafx.scene.Node
import javafx.scene.control.Button
import javafx.scene.control.ScrollPane
import javafx.scene.input.MouseEvent
import javafx.scene.input.ScrollEvent
import javafx.scene.layout.GridPane
import javafx.scene.layout.HBox
import javafx.scene.layout.Pane
import javafx.scene.layout.StackPane
import javafx.scene.paint.Color
import javafx.scene.shape.Rectangle
import javafx.scene.text.Font
import javafx.scene.text.Text
import net.maizegenetics.dna.factor.FactorTable
import net.maizegenetics.dna.factor.site.SNPSite

/**
 * @author Terry Casstevens
 * Created November 16, 2018
 */

private const val RECTANGLE_SIZE = 15.0
private const val LABEL_LENGTH = 10.0

private const val MAX_SCALE = 10.0
private const val MIN_SCALE = 0.1
private const val DELTA = 1.2

private val fonts = HashMap<Double, Font>()

private fun font(size: Double): Font {
    var result = fonts[size]
    if (result == null) {
        result = Font.font("Monospaced", size * 0.8)
        fonts[size] = result
        return result
    } else {
        return result
    }
}

private fun addText(parent: Pane, text: String, size: Double) {
    if (size >= 7.0) {
        val text = Text(text)
        text.font = font(size)
        parent.children += text
    }
}

private fun getVerticalText(text: String, size: Double): Node {
    val result = StackPane(Rectangle(size, LABEL_LENGTH, Color.TRANSPARENT))
    result.alignment = Pos.CENTER
    result.padding = Insets(10.0, 0.0, 10.0, 0.0)
    val text = Text(text)
    text.font = font(size)
    text.rotate = 90.0
    result.children += Group(text)
    return result
}

private fun NUCLEOTIDE_A_RECTANGLE(size: Double): Node {
    val result = StackPane(Rectangle(size, size, Color.RED))
    addText(result, "A", size)
    return result
}

private fun NUCLEOTIDE_C_RECTANGLE(size: Double): Node {
    val result = StackPane(Rectangle(size, size, Color.GREEN))
    addText(result, "C", size)
    return result
}

private fun NUCLEOTIDE_G_RECTANGLE(size: Double): Node {
    val result = StackPane(Rectangle(size, size, Color.YELLOW))
    addText(result, "G", size)
    return result
}

private fun NUCLEOTIDE_T_RECTANGLE(size: Double): Node {
    val result = StackPane(Rectangle(size, size, Color.BLUE))
    addText(result, "T", size)
    return result
}

private fun NUCLEOTIDE_UNKNOWN_RECTANGLE(size: Double): Node {
    val result = StackPane(Rectangle(size, size, Color.GREY))
    return result
}

private fun TAXA_RECTANGLE(name: String, size: Double): Node {
    val result = StackPane(Rectangle(LABEL_LENGTH, size, Color.TRANSPARENT))
    result.alignment = Pos.CENTER
    result.padding = Insets(0.0, 10.0, 0.0, 10.0)
    val text = Text(name)
    text.font = font(size)
    result.children += text
    return result
}

class FactorTableViewer private constructor(val table: FactorTable) {

    private val grid = GridPane()
    private val taxa = GridPane()
    private val factors = GridPane()
    val view = GridPane()
    val scale: DoubleProperty = SimpleDoubleProperty(1.0)
    val width = SimpleDoubleProperty(1000.0)
    private var lastWidth = 0.0

    private val currentSite = 0

    init {

        scale.addListener { observable, oldValue, newValue ->
            println("scale changed: $newValue")
            if (oldValue != newValue) {
                update()
            }
        }

        width.addListener { observable, oldValue, newValue ->
            println("width: ${newValue.toDouble()}")
            update()
        }

        grid.widthProperty().addListener { observable, oldValue, newValue ->
            println("observable: ${observable.value}  oldValue: $oldValue  newValue: $newValue")
            //update(newValue.toDouble())
        }

        taxa.isGridLinesVisible = true
        factors.isGridLinesVisible = true

        val corner = Rectangle(LABEL_LENGTH, LABEL_LENGTH, Color.TRANSPARENT)
        corner.widthProperty().bind(taxa.widthProperty())
        view.add(corner, 0, 0)

        view.add(factors, 1, 0)

        val scrollContent = HBox(taxa, grid)
        val scroll = ScrollPane(scrollContent)
        scroll.hbarPolicy = ScrollPane.ScrollBarPolicy.NEVER

        view.add(scroll, 0, 1, 2, 1)

    }

    fun update() {

        // TODO("Check is this works")
        if (!grid.isVisible) return

        val width = this.width.value

        if (Math.abs(lastWidth - width) < 30) return

        val squareSize = RECTANGLE_SIZE * scale.get()

        taxa.children.clear()
        addTaxa(squareSize)

        val gridWidth = width - (taxa.width + 30.0)

        val maxColumns = Math.floor(gridWidth / squareSize).toInt()

        println("width: $width  taxa width: ${taxa.width}  max columns: $maxColumns")

        grid.children.clear()
        addValues(squareSize, maxColumns)

        factors.children.clear()
        addFactors(maxColumns, squareSize)

    }

    val controls: HBox by lazy {

        val zoomIn = Button("+")
        zoomIn.onAction = EventHandler {
            scale.value = Math.min(MAX_SCALE, scale.value * DELTA)
        }

        val zoomOut = Button("-")
        zoomOut.onAction = EventHandler {
            scale.value = Math.max(MIN_SCALE, scale.value / DELTA)
        }

        val result = HBox(zoomOut, zoomIn)
        result.spacing = 10.0
        result.padding = Insets(3.0)
        result.alignment = Pos.CENTER_RIGHT

        result

    }


    private fun addTaxa(height: Double) {

        taxa.isGridLinesVisible = true

        table.taxa.forEachIndexed { index, taxon ->
            taxa.add(TAXA_RECTANGLE(taxon.name, height), 0, index)
        }

    }

    private fun addFactors(maxColumns: Int, width: Double) {

        factors.isGridLinesVisible = true

        var count = 0
        for (i in currentSite until Math.min(table.factors.size, currentSite + maxColumns)) {
            val factor = table.factors[i]
            val text = "${factor.startChr}:${factor.startPos}"
            val result = getVerticalText(text, width)
            factors.add(result, i - currentSite, 0)
            count++
        }
        println("number factors: $count")

    }

    private fun addValues(squareSize: Double, maxColumns: Int) {

        table.sites().limit(maxColumns.toLong()).forEach {
            val site = it.index
            when (it) {
                is SNPSite -> {
                    for (taxon in 0 until it.taxa.numberOfTaxa()) {
                        val nucleotide = it.genotypeAsString(taxon)
                        when (nucleotide) {
                            "A" -> grid.add(NUCLEOTIDE_A_RECTANGLE(squareSize), site, taxon)
                            "C" -> grid.add(NUCLEOTIDE_C_RECTANGLE(squareSize), site, taxon)
                            "G" -> grid.add(NUCLEOTIDE_G_RECTANGLE(squareSize), site, taxon)
                            "T" -> grid.add(NUCLEOTIDE_T_RECTANGLE(squareSize), site, taxon)
                            else -> grid.add(NUCLEOTIDE_UNKNOWN_RECTANGLE(squareSize), site, taxon)
                        }
                    }
                }
            }
        }

    }

    fun setPivot(x: Double, y: Double) {
        println("scale: ${scale.get()}")
        view.translateX = view.translateX - x
        view.translateY = view.translateY - y
    }

    companion object {

        private val INSTANCES = object : LinkedHashMap<FactorTable, FactorTableViewer>() {
            override fun removeEldestEntry(eldest: MutableMap.MutableEntry<FactorTable, FactorTableViewer>?): Boolean {
                return size > 5;
            }
        }

        @JvmStatic
        fun instance(table: FactorTable, width: ReadOnlyDoubleProperty): FactorTableViewer {
            var result = INSTANCES[table]
            if (result == null) {
                result = FactorTableViewer(table)
                result.width.bind(width)
                INSTANCES[table] = result
            }
            return result
        }

    }

}

internal class DragContext {

    var mouseAnchorX: Double = 0.0
    var mouseAnchorY: Double = 0.0

    var translateAnchorX: Double = 0.0
    var translateAnchorY: Double = 0.0

}

internal class SceneGestures(val viewer: FactorTableViewer) {

    private val canvas = viewer.view
    private val sceneDragContext = DragContext()

    // right mouse button => panning
    val onMousePressedEventHandler: EventHandler<MouseEvent> = EventHandler { event ->
        if (!event.isSecondaryButtonDown)
            return@EventHandler

        sceneDragContext.mouseAnchorX = event.sceneX
        sceneDragContext.mouseAnchorY = event.sceneY

        sceneDragContext.translateAnchorX = canvas.getTranslateX()
        sceneDragContext.translateAnchorY = canvas.getTranslateY()
    }

    // right mouse button => panning
    val onMouseDraggedEventHandler: EventHandler<MouseEvent> = EventHandler { event ->
        if (!event.isSecondaryButtonDown)
            return@EventHandler

        canvas.setTranslateX(sceneDragContext.translateAnchorX + event.sceneX - sceneDragContext.mouseAnchorX)
        canvas.setTranslateY(sceneDragContext.translateAnchorY + event.sceneY - sceneDragContext.mouseAnchorY)

        event.consume()
    }

    /**
     * Mouse wheel handler: zoom to pivot point
     */
    // currently we only use Y, same value is used for X
    // note: pivot value must be untransformed, i. e. without scaling
    val onScrollEventHandler: EventHandler<ScrollEvent> = EventHandler { event ->

        println("onScrollEventHandler: deltaY: ${event.deltaY}")

        var scale = viewer.scale.get()
        val oldScale = scale

        if (event.deltaY < 0)
            scale = reduce(scale)
        else
            scale = increase(scale)

        val f = scale / oldScale - 1

        val dx = event.sceneX - (canvas.boundsInParent.width / 2 + canvas.boundsInParent.minX)
        val dy = event.sceneY - (canvas.boundsInParent.height / 2 + canvas.boundsInParent.minY)

        viewer.scale.set(scale)
        //viewer.setPivot(f * dx, f * dy)

        event.consume()
    }

    companion object {

        val MAX_SCALE = 10.0
        val MIN_SCALE = 0.1
        private val DELTA = 1.2

        fun reduce(scale: Double) = Math.max(MIN_SCALE, scale / DELTA)

        fun increase(scale: Double) = Math.min(MAX_SCALE, scale * DELTA)

    }
}