package net.maizegenetics.dna.factor

import net.maizegenetics.dna.factor.site.FeatureSite
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 13, 2018
 */

const val UNKNOWN_ALLELE = 0xFF.toByte()

class FeatureTable(val taxa: TaxaList, private var sites: List<FeatureSite>) : List<FeatureSite> by sites {

    init {
        sites = sites.sorted()
    }

    fun numTaxa() = taxa.size

    fun numFeatures() = sites.size

    fun taxa() = taxa

    fun site(index: Int): FeatureSite = sites[index]

}