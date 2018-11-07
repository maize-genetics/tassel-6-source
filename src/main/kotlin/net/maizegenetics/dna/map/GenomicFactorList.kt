package net.maizegenetics.dna.map

import com.google.common.collect.ImmutableSortedSet
import java.util.*

/**
 * @author Terry Casstevens
 * Created November 07, 2018
 *
 * This is a list of genomic factors paired with it's weight (0.0 to 1.0)
 */

class GenomicFactorList private constructor(val factors: ImmutableSortedSet<Pair<GenomicFactor, Float>>) : SortedSet<Pair<GenomicFactor, Float>> by factors {

    class Builder() {

        val builder = ImmutableSortedSet.Builder<Pair<GenomicFactor, Float>>(kotlin.Comparator { o1, o2 ->
            o1.first.compareTo(o2.first)
        })

        fun add(factor: GenomicFactor) {
            builder.add(Pair(factor, 1.0f))
        }

        fun add(factor: GenomicFactor, weight: Float) {
            builder.add(Pair(factor, weight))
        }

        fun build(): GenomicFactorList {
            return GenomicFactorList(builder.build())
        }

    }

}
