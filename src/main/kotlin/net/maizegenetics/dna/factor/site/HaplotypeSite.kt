package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.util.SuperByteMatrix

/**
 * @author Terry Casstevens
 * Created April 6, 2021
 */

class HaplotypeSite(factor: GenomicFactor, taxa: TaxaList, private val strStates: Array<String>, private val genotypes: SuperByteMatrix, val ploidy: Int = 2, weight: Double? = null, isPhased: Boolean = false) : FactorSite(factor, taxa, weight, isPhased) {

    init {
        require(taxa.size == genotypes.numRows) { "Number of taxa: ${taxa.size} should match number of genotypes: ${genotypes.numRows}." }
    }

    override fun ploidy() = ploidy

    override fun genotype(taxon: Int): ByteArray {
        return genotypes.getAllColumns(taxon)
    }

    override fun genotypeAsString(taxon: Int): String {
        return genotypes.getAllRows(taxon).joinToString(if (isPhased) "|" else "/") { strStates[it.toInt()] }
    }

}