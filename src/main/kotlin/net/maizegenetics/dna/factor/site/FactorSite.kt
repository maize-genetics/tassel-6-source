package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 14, 2018
 */

abstract class FactorSite(val factor: GenomicFactor, val taxa: TaxaList, val weight: Double? = null, val isPhased: Boolean = false) : Comparable<FactorSite> {

    val stats = lazy { AlleleFreqCache.allelesSortedByFrequencyAndCountsNucleotide(0, TODO()) }

    abstract fun ploidy(): Int

    abstract fun genotype(taxon: Int): ByteArray

    abstract fun genotypeAsString(taxon: Int): String

    override fun compareTo(other: FactorSite): Int {
        return factor.compareTo(other.factor)
    }

}