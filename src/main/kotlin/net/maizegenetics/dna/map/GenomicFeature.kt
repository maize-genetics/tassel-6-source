package net.maizegenetics.dna.map

import com.google.common.collect.Range

/**
 * @author Terry Casstevens
 * Created November 07, 2018
 */

enum class Feature {
    SNP,
    Haplotype,
    Binary,
    Categorical,
    Expression,
    Depth,
    Dosage
}

open class GenomicFeature(val startChr: Chromosome, val startPos: Int, val endChr: Chromosome = startChr, val endPos: Int = startPos) : Comparable<GenomicFeature> {

    override fun compareTo(other: GenomicFeature): Int {
        val result = startChr.compareTo(other.startChr)
        if (result != 0) return result
        return startPos.compareTo(other.startPos)
    }

    val range: Range<ChrPos> by lazy { Range.closed(ChrPos(startChr, startPos), ChrPos(endChr, endPos)) }

}

data class ChrPos(val chromosome: Chromosome, val position: Int) : Comparable<ChrPos> {

    override fun compareTo(other: ChrPos): Int {
        val result = chromosome.compareTo(other.chromosome)
        if (result != 0) return result
        return position.compareTo(other.position)
    }

}
