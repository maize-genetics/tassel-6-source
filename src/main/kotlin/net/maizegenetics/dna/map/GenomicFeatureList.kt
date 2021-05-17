package net.maizegenetics.dna.map

import com.google.common.collect.ImmutableList

/**
 * @author Terry Casstevens
 * Created November 07, 2018
 *
 * This is a list of genomic factors paired with it's weight
 */

class GenomicFeatureList private constructor(private val factors: ImmutableList<GenomicFeature>) : List<GenomicFeature> by factors {

    // TODO("Sort")
    // TODO("Validate Order")

    class Builder {

        val builder = ImmutableList.builder<GenomicFeature>()

        fun add(factor: GenomicFeature) {
            builder.add(factor)
        }

        fun addAll(factors: List<GenomicFeature>) {
            factors.forEach { add(it) }
        }

        fun build(): GenomicFeatureList {
            return GenomicFeatureList(builder.build())
        }

    }

}
