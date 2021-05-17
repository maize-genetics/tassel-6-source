package net.maizegenetics.dna.factor

import net.maizegenetics.dna.factor.site.FactorSite
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 13, 2018
 */

const val UNKNOWN_ALLELE = 0xFF.toByte()

class FactorTable(val taxa: TaxaList, private var sites: List<FactorSite>) : List<FactorSite> by sites {

    init {
        sites = sites.sorted()
    }

    fun numTaxa() = taxa.size

    fun numFactors() = sites.size

    fun taxa() = taxa

    fun site(index: Int): FactorSite = sites[index]

}