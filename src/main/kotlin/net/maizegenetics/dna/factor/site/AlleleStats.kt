package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.factor.UNKNOWN_ALLELE
import java.util.*

class AlleleStats(site: FeatureSite) {

    val alleleCounts: List<AlleleCount>
    val numAlleles: Int
    val totalNonMissingAlleles: Int

    init {

        val numTaxa = site.taxa.numberOfTaxa()
        val alleleFreq = IntArray(256)
        for (taxon in 0 until numTaxa) {
            val alleles = site.genotype(taxon)
            alleles.forEach { allele ->
                alleleFreq[allele.toInt() and 0xFF]++
            }
        }

        alleleCounts = sort(alleleFreq)

        numAlleles = alleleCounts.size

        totalNonMissingAlleles = alleleCounts
                .map { it.count }
                .sum()

    }

    fun majorAllele(): Byte {
        return if (numAlleles > 0) {
            alleleCounts.first().allele
        } else {
            UNKNOWN_ALLELE
        }
    }

    fun majorAlleleFrequency(): Double {
        return if (numAlleles > 0) {
            alleleCounts.first().count.toDouble() / totalNonMissingAlleles.toDouble()
        } else {
            0.0
        }
    }

    fun minorAllele(): Byte {
        return if (numAlleles > 1) {
            alleleCounts[1].allele
        } else {
            UNKNOWN_ALLELE
        }
    }

    fun minorAlleleFrequency(): Double {
        return if (numAlleles > 1) {
            alleleCounts[1].count.toDouble() / totalNonMissingAlleles.toDouble()
        } else {
            0.0
        }
    }

    fun totalGametesNonMissingForSite(): Int {
        return alleleCounts
                .map { it.count }
                .sum()
    }

}

private fun sort(data: IntArray): List<AlleleCount> {

    return data
            // don't include last index
            // that represents FactorTableKt.UNKNOWN_ALLELE
            .dropLast(1)
            .mapIndexed { allele, count -> AlleleCount(allele.toByte(), count) }
            .filter { it.count != 0 }
            .sortedWith(Comparator { t, t2 -> t2.count.compareTo(t.count) })

}

data class AlleleCount(val allele: Byte, val count: Int)