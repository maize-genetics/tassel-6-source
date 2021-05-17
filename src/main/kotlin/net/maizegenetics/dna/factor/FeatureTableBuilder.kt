package net.maizegenetics.dna.factor

import net.maizegenetics.dna.factor.site.FeatureSite
import net.maizegenetics.taxa.TaxaList

class FeatureTableBuilder constructor(val taxa: TaxaList) {

    private val sites = mutableListOf<FeatureSite>()

    fun add(site: FeatureSite) {
        sites.add(site)
    }

    fun build(): FeatureTable {
        sites.sort()
        return FeatureTable(taxa, sites)
    }

}