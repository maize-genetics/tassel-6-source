package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.GenomicFeature
import net.maizegenetics.taxa.TaxaList

class SNPSiteBuilder constructor(val factor: GenomicFeature, val taxa: TaxaList) {

    private val genotypes = ByteArray(taxa.size)

    var isPhased = false

    fun set(taxon: Int, value: Byte): SNPSiteBuilder {
        genotypes[taxon] = value
        return this
    }

    fun build() = SNPSite(factor, taxa, genotypes, isPhased = isPhased)

}