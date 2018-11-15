package net.maizegenetics.dna.factor.factorvalues

import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache
import net.maizegenetics.dna.snp.genotypecall.Stats
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 14, 2018
 */

class FactorSite (val factor: GenomicFactor, val index: Int, val taxa: TaxaList, val values: Array<Array<Byte>>) {

    val stats = lazy { AlleleFreqCache.allelesSortedByFrequencyAndCountsNucleotide(index, TODO()) }

}