package net.maizegenetics.dna.map

import com.google.common.collect.ImmutableList

/**
 * @author Terry Casstevens
 * Created November 07, 2018
 *
 * This is a list of genomic factors paired with it's weight
 */

class GenomicFactorList private constructor(private val factors: ImmutableList<GenomicFactor>) : List<GenomicFactor> by factors {

    // TODO("Sort")
    // TODO("Validate Order")

    class Builder {

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

}
