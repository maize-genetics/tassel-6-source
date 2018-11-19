package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 14, 2018
 */

open class FactorSite (val factor: GenomicFactor, val index: Int, val taxa: TaxaList) {

    val stats = lazy { AlleleFreqCache.allelesSortedByFrequencyAndCountsNucleotide(index, TODO()) }

}