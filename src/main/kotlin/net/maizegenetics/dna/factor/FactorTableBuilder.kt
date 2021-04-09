package net.maizegenetics.dna.factor

import net.maizegenetics.dna.factor.site.FactorSite
import net.maizegenetics.taxa.TaxaList

class FactorTableBuilder constructor(val taxa: TaxaList) {

    private val sites = mutableListOf<FactorSite>()

    fun add(site: FactorSite) {
        sites.add(site)
    }

    fun build(): FactorTable {
        sites.sort()
        return FactorTable(taxa, sites)
    }

}