package net.maizegenetics.dna.map

import com.google.common.collect.Range

/**
 * @author Terry Casstevens
 * Created November 07, 2018
 */

enum class Factor {
    SNP,
    Haplotype,
    Binary,
    Categorical,
    Expression,
    Depth,
    Dosage
}

class GenomicFactor(val factor: Factor, val chromosome: Chromosome, val startPos: Int, val endPos: Int = startPos) : Comparable<GenomicFactor> {

    override fun compareTo(other: GenomicFactor): Int {
        val result = chromosome.compareTo(other.chromosome)
        if (result != 0) return result
        return startPos.compareTo(other.startPos)
    }

    val range: Range<ChrPos> by lazy { Range.closed(ChrPos(chromosome, startPos), ChrPos(chromosome, endPos)) }

}

data class ChrPos(val chromosome: Chromosome, val position: Int) : Comparable<ChrPos> {

    override fun compareTo(other: ChrPos): Int {
        val result = chromosome.compareTo(other.chromosome)
        if (result != 0) return result
        return position.compareTo(other.position)
    }

}
