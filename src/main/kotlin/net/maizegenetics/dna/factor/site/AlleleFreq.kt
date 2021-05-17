@file:JvmName("AlleleFreq")

package net.maizegenetics.dna.factor.site

import java.util.*
import kotlin.Comparator

fun alleleFreq(site: FeatureSite, maxNumAlleles: Int): Array<IntArray>? {

    val numTaxa = site.taxa.numberOfTaxa()
    val alleleFreq = IntArray(maxNumAlleles)
    for (taxon in 0 until numTaxa) {
        val alleles = site.genotype(taxon)
        alleles
                //  and 0xFF converts byte to unsigned integer
                .filter { allele -> allele.toInt() and 0xFF < maxNumAlleles }
                .forEach { allele ->
                    //println("allele: $allele   allele int: ${allele.toInt()}")
                    alleleFreq[allele.toInt() and 0xFF]++
                }
    }

    val sortedAlleles = sort(alleleFreq, maxNumAlleles)
    val numAlleles = sortedAlleles.size
    //println("alleleFreq: numAlleles: $numAlleles")
    val alleleCounts = Array(2) { IntArray(numAlleles) }
    sortedAlleles.forEachIndexed { index, pair ->
        alleleCounts[0][index] = pair.first
        alleleCounts[1][index] = pair.second
    }

    return alleleCounts

}

private fun sort(data: IntArray, maxNumAlleles: Int): SortedSet<Pair<Int, Int>> {

    return data
            .mapIndexed { allele, count -> Pair(allele, count) }
            .filter { it.second != 0 }
            .take(maxNumAlleles)
            .toSortedSet(Comparator { t, t2 -> t2.second.compareTo(t.second) })

}
