package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.SNP
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 16, 2018
 */

class SNPSite(factor: SNP, index: Int, taxa: TaxaList, private val values: ByteArray) : FactorSite(factor, index, taxa) {

    fun genotype(taxon: Int) = values[taxon]

    fun genotypeAsString(taxon: Int) = NucleotideAlignmentConstants.getNucleotideIUPAC(values[taxon])

}