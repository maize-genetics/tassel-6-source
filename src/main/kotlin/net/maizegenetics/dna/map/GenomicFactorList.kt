package net.maizegenetics.dna.map

import com.google.common.collect.ImmutableList
import com.google.common.collect.ImmutableSortedSet

/**
 * @author Terry Casstevens
 * Created November 07, 2018
 *
 * This is a list of genomic factors paired with it's weight
 */

class GenomicFactorList private constructor(private val factors: ImmutableList<GenomicFactor>, val weights: ImmutableList<Double>? = null) : List<GenomicFactor> by factors {

    // TODO("Sort")
    // TODO("Validate Order")

    class Builder() {

        val builder = ImmutableList.builder<GenomicFactor>()

        fun add(factor: GenomicFactor) {
            builder.add(factor)
        }

        fun addAll(factors: List<GenomicFactor>) {
            factors.forEach { add(it) }
        }

        fun build(): GenomicFactorList {
            return GenomicFactorList(builder.build())
        }

    }

    class BuilderWithWeights() {

        val builder = ImmutableSortedSet.Builder<Pair<GenomicFactor, Double>>(kotlin.Comparator { o1, o2 ->
            o1.first.compareTo(o2.first)
        })

        fun add(factor: GenomicFactor, weight: Double) {
            builder.add(Pair(factor, weight))
        }

        fun build(): GenomicFactorList {
            val sortedSet = builder.build()
            val factors = ImmutableList.builder<GenomicFactor>()
            val weights = ImmutableList.builder<Double>()
            sortedSet.forEach {
                factors.add(it.first)
                weights.add(it.second)
            }
            return GenomicFactorList(factors.build(), weights.build())
        }

    }

}
