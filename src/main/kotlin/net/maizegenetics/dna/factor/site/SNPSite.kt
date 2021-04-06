package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.dna.snp.GenotypeTableUtils
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 16, 2018
 */

class SNPSite(factor: GenomicFactor, taxa: TaxaList, private val values: ByteArray, weight: Double? = null, isPhased: Boolean = false) : FactorSite(factor, taxa, weight, isPhased) {

    override fun ploidy() = 2

    override fun genotype(taxon: Int): ByteArray = GenotypeTableUtils.getDiploidValues(values[taxon])

    override fun genotypeAsString(taxon: Int): String = NucleotideAlignmentConstants.getNucleotideIUPAC(values[taxon])

}