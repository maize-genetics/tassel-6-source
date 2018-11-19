package net.maizegenetics.dna.factor

import net.maizegenetics.dna.factor.site.FactorSite
import net.maizegenetics.dna.map.GenomicFactorList
import net.maizegenetics.taxa.TaxaList
import java.util.*
import java.util.function.Consumer
import java.util.stream.Stream
import java.util.stream.StreamSupport

/**
 * @author Terry Casstevens
 * Created November 13, 2018
 */

abstract class FactorTable(val taxa: TaxaList, val factors: GenomicFactorList) {

    fun sites(): Stream<FactorSite> {
        return StreamSupport.stream(SiteSpliterator(this), false)
    }

    internal class SiteSpliterator(val table: FactorTable) : Spliterator<FactorSite> {

        private var index = 0
        private val numFactors = table.factors.size

        override fun estimateSize(): Long {
            return numFactors.toLong()
        }

        override fun characteristics(): Int {
            return Spliterator.SIZED
        }

        override fun tryAdvance(action: Consumer<in FactorSite>): Boolean {

            return if (index >= numFactors) {
                false
            } else {
                action.accept(table.site(index++))
                true
            }

        }

        override fun forEachRemaining(action: Consumer<in FactorSite>) {

            while (index < numFactors) {
                action.accept(table.site(index++))
            }

        }

        override fun trySplit(): Spliterator<FactorSite>? {
            return null
        }

    }

    abstract protected fun site(index: Int): FactorSite

}